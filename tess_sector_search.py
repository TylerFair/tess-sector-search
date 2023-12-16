import os
import argparse
from datetime import datetime
import pandas as pd
import requests
import re
import subprocess
import time
from pbs_script_generator import PBSScript

class SectorSetup:
    """
    The SectorSetup class handles the setup and preparation of a TESS sector for data analysis.
    This includes downloading cadence files, retrieving spacecraft pointings, creating necessary directories,
    configuring light curve (LC) files, and optionally generating a catalog.

    Attributes:
        sector (int): The TESS sector number.
        orbit (int): The TESS orbit number within the sector.
        path (str): The base path for file and directory creation.
        pbs_script (PBSScript): An instance of the PBSScript class for generating PBS scripts.
        cadence_url_pattern (str): The URL pattern for downloading cadence files.
        pointings_url (str): The URL for retrieving spacecraft pointings data.
    """
    def __init__(self, sector, orbit, path="30daytemp/"):
        self.sector = sector
        self.orbit = orbit
        self.path = path
        self.pbs_script = PBSScript(sector, orbit, self.path)  # Creating an instance of PBSScript
        self.cadence_url_pattern = f'https://archive.stsci.edu/hlsps/tica/bundles/s00{sector}/hlsp_tica_tess_ffi_s00{sector}-o{{orbit}}{{orbit_type}}-cam1-ccd1_tess_v01_ffis.sh'
        self.pointings_url = f'https://tess.mit.edu/observations/sector-{sector}/'


    def download_file(self, url, filename):
        response = requests.get(url)
        if response.status_code == 200:
            with open(filename, 'wb') as file:
                file.write(response.content)
                print("File downloaded: %s" % filename)
        else:
            print("Failed to download file. Status code: %d" % response.status_code)

    def extract_cadence_number(self, file_content, orbit_type):
        cadence_match = re.findall(r'--output.*?-o\d+-(\d+)-cam1-ccd1_tess_v01_img.fits', file_content)
        if orbit_type == 'a':
            return cadence_match[0].lstrip('0') if cadence_match else None
        elif orbit_type == 'b':
            return cadence_match[-1].lstrip('0') if cadence_match else None

    def get_cadence_limits(self):
        cadences_csv = f"S{self.sector}_cadences.csv"
        cadences_csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), cadences_csv)

        if not os.path.exists(cadences_csv_path):
            data = {'orbit': [], 'start': [], 'end': []}
            for orbit in ['1', '2']:
                for orbit_type in ['a', 'b']:
                    url = self.cadence_url_pattern.format(sector=self.sector, orbit=orbit, orbit_type=orbit_type)
                    filename = f"downloaded_file_{self.sector}_{orbit}{orbit_type}.sh"
                    self.download_file(url, filename)
                    with open(filename, 'r') as file:
                        file_content = file.read()
                    cadence_number = self.extract_cadence_number(file_content, orbit_type)
                    key = 'start' if orbit_type == 'a' else 'end'
                    data[key].append(int(cadence_number) if cadence_number else None)
                    os.remove(filename)
                data['orbit'].append(int(orbit))
            
            cadences = pd.DataFrame(data)
            cadences.to_csv(cadences_csv_path, index=False)
            print(f"CSV file '{cadences_csv}' created.")
        else:
            print(f"{cadences_csv} already exists, skipping creation.")
            cadences = pd.read_csv(cadences_csv_path)
        orbit_index = 0 if self.orbit % 2 != 0 else 1
        return cadences.iloc[orbit_index]['start'], cadences.iloc[orbit_index]['end']

    def get_pointings(self):
        """
        Retrieves the spacecraft pointing information for a given TESS sector.
        Downloads the pointing data from the TESS website and processes it.

        Returns:
            tuple: A tuple containing the right ascension (ra), declination (dec), and roll of the spacecraft.
        """
        pointings_csv = "S%d.csv" % self.sector
        pointings_csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), pointings_csv)

        # Check if the pointings file already exists
        if not os.path.exists(pointings_csv_path):
            response = requests.get(self.pointings_url)
            if response.status_code == 200:
                dfs = pd.read_html(response.text, flavor='bs4', header=0)
                sector_df = dfs[0]
                spacecraft_ra, spacecraft_dec, spacecraft_roll = sector_df.iloc[0, [1, 2, 3]]
                modified_df = sector_df.iloc[1:, 1:-2].apply(pd.to_numeric, errors='coerce').dropna()
                modified_df['camera'] = range(1, len(modified_df) + 1)
                modified_df = modified_df.reset_index(drop=True)
                modified_df.columns = ['ra', 'dec', 'roll', 'camera']
                self.save_to_csv(modified_df, "S%d.csv" % self.sector)
                return spacecraft_ra, spacecraft_dec, spacecraft_roll
            else:
                print("Failed to retrieve data. Status code: %d" % response.status_code)
                return None, None, None
        else:
            print("Already have pointing file, only retrieving spacecraft.")
            response = requests.get(self.pointings_url)
            if response.status_code == 200:
                dfs = pd.read_html(response.text, flavor='bs4', header=0)
                sector_df = dfs[0]
                spacecraft_ra, spacecraft_dec, spacecraft_roll = sector_df.iloc[0, [1, 2, 3]]
                return spacecraft_ra, spacecraft_dec, spacecraft_roll
            
    def save_to_csv(self, data, filename):
        """
        Saves the given data to a CSV file.

        Parameters:
            data (DataFrame): The pandas DataFrame to be saved.
            filename (str): The name of the CSV file where the data will be saved.
        """
        data.to_csv(filename, index=False)
        print("Data saved to %s" % filename)

    def create_lc_config_files(self, run_path):
        """
        Creates light curve (LC) configuration files for each camera.

        Parameters:
            run_path (str): The directory path where the LC configuration files will be created.
        """
        for cam in range(1, 5):
            config_file_path = os.path.join(run_path, "example-lc-cam%d.cfg" % cam)
            with open(config_file_path, "w") as file:
                file.write("[IOSettings]\n")
                file.write("indir = %s/\n" % os.path.join(self.path, "orbit-%d" % self.orbit))
                file.write("orbit_id = %d\n" % self.orbit)
                file.write("cadence_type = 10\n")
                file.write("sector = %d\n" % self.sector)
                file.write("inlist = lc.ls\n")
                file.write("\n[BLS]\n")
                file.write("# IO settings\n")
                file.write("coljd = 1\n")
                file.write("colmag= 2\n")
                file.write("f0 = 0.0714\n")
                file.write("f1 = 10.\n")
                file.write("fn = 500000\n")
                file.write("\n[LC]\n")
                file.write("# IO settings\n")
                file.write("method = ksp\n")
                file.write("bkspace_min=0.3\n")
                file.write("quadfile = \"cam%d_quat.txt\"\n" % cam)
                file.write("\n[Phot]\n")
                file.write("# IO settings\n")
                file.write("aps='1.75:4.0:3.0,2.5:4.0:3.0,3.0:4.0:3.0,3.5:5.0:4.0,8.0:10.0:5.0'\n")

    def create_directory(self, spacecraft_ra, spacecraft_dec, spacecraft_roll):
        """
        Creates the directory structure necessary for TESS data analysis. This includes directories
        for each camera and CCD, as well as configuration files.

        Parameters:
            spacecraft_ra (float): Right ascension of the spacecraft.
            spacecraft_dec (float): Declination of the spacecraft.
            spacecraft_roll (float): Roll angle of the spacecraft.
        """
        current_date = datetime.now().strftime("%Y%m%d")
        start_cadence, end_cadence = self.get_cadence_limits()
        orbit_path = os.path.join(self.path, "orbit-%d" % self.orbit)
        os.makedirs(orbit_path, exist_ok=True)
        ffi_path = os.path.join(orbit_path, "ffi")
        os.makedirs(ffi_path, exist_ok=True)
        run_path = os.path.join(ffi_path, "run")
        os.makedirs(run_path, exist_ok=True)
        # Create directories for each camera and CCD
        for cam in range(1, 5):
            for ccd in range(1, 5):
                for data_type in ['astrom', 'FITS', 'LC', 'phot', 'sub', 'subphot']:
                    dir_path = os.path.join(ffi_path, "cam%d" % cam, "ccd%d" % ccd, data_type)
                    os.makedirs(dir_path, exist_ok=True)
        # Create the configuration file for FFI
        config_file_path = os.path.join(run_path, "example-ffi.cfg")
        with open(config_file_path, "w") as file:
            file.write("[Setup]\n")
            file.write("indir= %s/ffi/\n" % os.path.join(self.path, "orbit-%d" % self.orbit))
            file.write("orbit_id = %d\n" % self.orbit)
            file.write("cadence_type = 10\n")
            file.write("cadence_limit = %d, %d\n" % (start_cadence, end_cadence))
            file.write("basename = tess%s\n" % current_date)
            file.write("ra0=%.3f\n" % float(spacecraft_ra))
            file.write("dec0=%.3f\n" % float(spacecraft_dec))
            file.write("roll=%.3f\n" % float(spacecraft_roll))
            file.write("[Catalog]\n")
            file.write("catfile='catalog.txt'\n")
            file.write("method='tic'\n")
            file.write("width=15.\n")
            file.write("ratiolim = 31.5\n")
            file.write("leading_cols= 45\n")
            file.write("buffer_rows= 1\n")
            file.write("maglim = 12\n")
            file.write("\n[Fistar]\n")
            file.write("threshold=83333\n")
            file.write("fiign=True\n")
            file.write("\n[Fiphot]\n")
            file.write("magtoflux=18\n")
            file.write("gain=1.0\n")
            file.write("skyfit_sigma=3\n")
            file.write("skyfit_niter=2\n")
            file.write("disjoint_radius=2\n")
            file.write("apertures='1.75:4.0:3.0,2.5:4.0:3.0,3.0:4.0:3.0,3.5:5.0:4.0,8.0:10.0:5.0'\n")
            file.write("\n[Grmatch]\n")
            file.write("order=4\n")
            file.write("unitarity=0.05\n")
            file.write("maxdistance=2\n")
            file.write("maxp=200\n")
 
        self.create_lc_config_files(run_path)

    def generate_catalog(self):
        """
        Generates a catalog of celestial objects based on the TESS sector's data.
        The catalog is used for further analysis in the data processing pipeline.
        """
        log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logfile')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        log_file = os.path.join(log_dir, 'catalog.out')

        # Assuming start_cadence is obtained from get_cadence_limits
        start_cadence, _ = self.get_cadence_limits()

        #subprocess.run(['python', 'catalog_query_tess.py', '--logfile', log_file, '--pointings', f"S{self.sector}.csv", '--orbitid', str(self.orbit), '--cadence', str(start_cadence), '--path', '.','--numsquares', '4', '--sector', str(self.sector)], check=True)     
        subprocess.run(['python', 'catalog_query_tess.py', '--logfile', log_file, '--pointings', f"S{self.sector}.csv", '--orbitid', str(self.orbit), '--cadence', str(start_cadence),'--numsquares', '4', '--sector', str(self.sector), '--path', str(self.path), '--maglim', 12], check=True)
   
    def create_sector_files(self):
        """
        Creates necessary files and directories for TESS sector processing.
        This includes cadence files, pointings, directories, and LC config files.
        """
        # Handle cadences
        start_cadence, end_cadence = self.get_cadence_limits()

        # Handle pointings
        spacecraft_ra, spacecraft_dec, spacecraft_roll = self.get_pointings()

        # Create directory structure
        self.create_directory(spacecraft_ra, spacecraft_dec, spacecraft_roll)

        # Generate LC config files
        run_path = os.path.join(self.path, f"orbit-{self.orbit}", "ffi", "run")
        self.create_lc_config_files(run_path)

        # Optional: Generate catalog if needed
        if self.catalog_required:
            self.generate_catalog()

    def generate_pbs_scripts(self):
            """
            Generates PBS scripts for downloading data, photometry processing,
            and light curve analysis using the PBSScript instance.
            """
            # Generate the download script
            self.pbs_script.create_download_script()

            # Generate the photometry script
            self.pbs_script.create_photometry_script()

            # Generate the light curve script
            self.pbs_script.create_light_curve_script()

class PBSRun:
    """
    The PBSRun class is responsible for executing PBS scripts generated by the PBSScript class.
    It allows for running download, photometry, and light curve scripts independently.

    Attributes:
        pbs_script (PBSScript): An instance of the PBSScript class.
    """
    def __init__(self, pbs_script):
        """
        Initializes the PBSRun with a PBSScript instance.

        Parameters:
            pbs_script (PBSScript): An instance of the PBSScript class.
        """
        self.pbs_script = pbs_script

    def download(self):
        """
        Executes the download script using the PBSScript instance.
        """
        if hasattr(self.pbs_script, 'create_download_script'):
            print("Running download script.")
            self.pbs_script.run_download_script()
        else:
            print("Download script not found or not generated.")

    def photometry(self):
        """
        Executes the photometry script using the PBSScript instance.
        """
        if hasattr(self.pbs_script, 'create_photometry_script'):
            print("Running photometry script.")
            self.pbs_script.run_photometry_script()
        else:
            print("Photometry script not found or not generated.")

    def lc(self):
        """
        Executes the light curve script using the PBSScript instance.
        """
        if hasattr(self.pbs_script, 'create_light_curve_script'):
            print("Running light curve script.")
            self.pbs_script.run_light_curve_script()
        else:
            print("Light curve script not found or not generated.")

def main():
    parser = argparse.ArgumentParser(description="Create directory structure for orbit data.")
    parser.add_argument("sector", type=int, help="Sector number")
    parser.add_argument("orbit", type=int, help="Orbit ID number")
    parser.add_argument("--path", default="30daytemp", help="Base path for the directories. Default is '30daytemp'")
    parser.add_argument("--catalog", action='store_true', help="Generate catalogs.")
    parser.add_argument("--download", action='store_true', help='Run Download PBS Script')
    parser.add_argument("--photometry", action='store_true', help='Run photometry PBS Script')
    parser.add_argument("--lc", action='store_true', help='Run light curve PBS Script')

    args = parser.parse_args()

    tess_processor = SectorSetup(args.sector, args.orbit, args.path)
    tess_processor.catalog_required = args.catalog  # Set the flag based on the user's choice

    # Create all sector files
    tess_processor.create_sector_files()

    # Generate PBS Scripts
    tess_processor.generate_pbs_scripts()
    # Initialize PBSRun with the created PBSScript instance
    pbs_run = PBSRun(tess_processor.pbs_script)

    # Count the number of scripts to run
    scripts_to_run = sum([args.download, args.photometry, args.lc])

    if args.download:
        pbs_run.download()
        if scripts_to_run > 1:
            print("Waiting 6 hours before running the next script")
            time.sleep(21600)  # 6 hours

    if args.photometry:
        pbs_run.photometry()
        # Check if light curve script is also scheduled to run
        if args.lc:
            print("Waiting 12 hours before running the light curve script")
            time.sleep(43200)  # 12 hours

    if args.lc:
        pbs_run.lc()

if __name__ == "__main__":
    main()

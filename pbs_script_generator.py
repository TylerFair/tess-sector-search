import os
import subprocess

class PBSScript:
    """
    A class to handle the creation and submission of various PBS scripts for tasks 
    such as downloading, photometry, and light curve processing for TESS (Transiting 
    Exoplanet Survey Satellite) data.

    Attributes:
        sector (int): The sector number for TESS data processing.
        orbit (int): The orbit number for TESS data processing.
        path (str): Base file path for script creation and data processing.
    """
    def __init__(self, sector, orbit, path='~/30daytemp/'):
        """
        Constructs all the necessary attributes for the PBSScript object.

        Parameters:
            sector (int): The sector number for TESS data processing.
            orbit (int): The orbit number for TESS data processing.
            path (str): Base file path for script creation and data processing.
        """
        self.sector = sector
        self.orbit = orbit
        self.path = os.path.expanduser(path)
        self.common_header = self.create_common_header()

    def create_common_header(self, user_email="default@usq.edu.au"):
        """
        Creates a common header for PBS scripts.

        Parameters:
            user_email (str): Email address for PBS notifications.

        Returns:
            str: A string containing the common header for PBS scripts.
        """
        header = "#!/bin/bash\n\n"
        header += "#PBS -l nodes=1:ppn=10\n"
        header += "#PBS -q default\n"
        header += "#PBS -j oe\n"
        header += "#PBS -P 210101893\n"
        header += "#PBS -m abe\n"
        header += f"#PBS -M {user_email}\n"
        return header

    def create_download_script(self):
        """
        Creates a script for downloading TESS data.

        The script is saved to a file whose path is constructed based on the sector and orbit.

        Returns:
            str: The path to the created download script.
        """
        output_file_name = f"downloadfits_S{self.sector}_o{self.orbit}.pbs"
        output_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), output_file_name)
        start_cadence, end_cadence = self.get_cadence_limits()

        orbit_binary = self.orbit%2
        if orbit_binary == 0:
            orbit_binary += 2
        
        with open(output_file_path, "w") as f:
            f.write(self.common_header)
            f.write("#PBS -N download_tica\n")
            f.write("#PBS -l walltime=06:00:00\n")
            f.write("#PBS -o logfile/curl.out\n")
            f.write("module load curl/7.80.0-gcc-6vy\n")
            for cam in range(1, 5):
                    for ccd in range(1, 5):
                        f.write("curl --proxy http://webproxy.usq.edu.au:8080 -f --parallel --parallel-max 10 --create-dirs --output '{self.path}/S%d/orbit-%d-tica/cam%d-ccd%d/hlsp_tica_tess_ffi_s%.4d-o%d-00#1-cam%d-ccd%d_tess_v01_img.fits' https://archive.stsci.edu/hlsps/tica/s%.4d/cam%d-ccd%d/hlsp_tica_tess_ffi_s%.4d-o%d-00[%d-%d]-cam%d-ccd%d_tess_v01_img.fits 2>&1 | tee curl.log\n" % (int(self.sector),int(self.orbit),cam,ccd,int(self.sector),orbit_binary,cam,ccd,int(self.sector),cam,ccd,int(self.sector),orbit_binary,int(start_cadence),int(end_cadence),cam,ccd))
        try:
            subprocess.run(['dos2unix', output_file_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Failed to convert {output_file_name} to Unix format: {e}")
        
        return output_file_path

    def create_photometry_script(self):
        """
        Creates a script for running photometry processing on TESS data.

        The script is saved to a file whose path is constructed based on the sector and orbit.

        Returns:
            str: The path to the created photometry script.
        """
        output_file_name = f"runphot_S{self.sector}_o{self.orbit}.pbs"
        output_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), output_file_name)

        with open(output_file_path, "w") as f:
            f.write(self.common_header)
            f.write("#PBS -N qlp_phot\n")
            f.write("#PBS -l walltime=12:00:00\n")
            f.write("#PBS -o logfile/phot.out\n")
            f.write("module unload python\n")
            f.write("module load python/3.9.13-gcc-jhn\n")
            f.write("export PATH=$PATH:\"/home/u8015661/FFITools-0.6.1/PATools/bin/\":\"/home/u8015661/FFITools-0.6.1/LCTools/bin/\"\n")
            f.write("export PYTHONPATH=$PYTHONPATH:\"/home/u8015661/FFITools-0.6.1/PATools/\":\"/home/u8015661/FFITools-0.6.1/LCTools/\"\n")
            f.write("export PYTHONPATH=$PYTHONPATH:\"/home/u8015661/tess-point/\"\n")
            f.write(f"runpath=\"{self.path}/orbit-{self.orbit}/ffi/run/\"\n")
            for cam in range(1, 5):
                for ccd in range(1, 5):
                    f.write(f"patools-rename-tica -c $runpath\"example-ffi.cfg\" -I \"{self.path}/S{self.sector}/orbit-{self.orbit}-tica\" -n 60 -a {cam},{ccd} --bgsub\n")
                    f.write(f"patools-readheader -c $runpath\"example-ffi.cfg\" -a {cam},{ccd} -o $runpath\"orbit{self.orbit}_header_cam{cam}ccd{ccd}.csv\"\n")
                    f.write(f"patools-mediansub -c $runpath\"example-ffi.cfg\" -n 60 -a {cam},{ccd} -q $runpath\n")
            f.write(f"patools-reduce fieldbgsub -c $runpath\"example-ffi.cfg\" -n 60 --debug\n")
        
        # Convert script to Unix format to remove ^M characters
        try:
            subprocess.run(['dos2unix', output_file_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Failed to convert {output_file_name} to Unix format: {e}")
        
        return output_file_path

    def create_light_curve_script(self):
        """
        Creates a script for generating light curves from TESS data.

        The script is saved to a file whose path is constructed based on the sector and orbit.

        Returns:
            str: The path to the created light curve script.
        """
        output_file_name = f"runlc_S{self.sector}_o{self.orbit}.pbs"
        output_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), output_file_name)

        with open(output_file_path, "w") as f:
            f.write(self.common_header)
            f.write("#PBS -N qlp_lc\n")
            f.write("#PBS -l walltime=12:00:00\n")
            f.write("#PBS -o logfile/lc.out\n")
            f.write("module unload python\n")
            f.write("module load python/3.9.13-gcc-jhn\n")
            f.write("export PATH=$PATH:\"/home/u8015661/FFITools-0.6.1/PATools/bin/\":\"/home/u8015661/FFITools-0.6.1/LCTools/bin/\"\n")
            f.write("export PYTHONPATH=$PYTHONPATH:\"/home/u8015661/FFITools-0.6.1/PATools/\":\"/home/u8015661/FFITools-0.6.1/LCTools/\"\n")
            f.write(f"runpath=\"{self.path}/orbit-{self.orbit}/ffi/run/\"\n")
            for cam in range(1, 5):
                for ccd in range(1, 5):
                    f.write(f"lctools-phottoh5-bgsub -c $runpath\"example-lc-cam{cam}.cfg\" -n 60 --debug --logfile - -a {cam},{ccd} 2>&1 | tee logfile/lcgen.log\n")

        # Convert script to Unix format to remove ^M characters
        try:
            subprocess.run(['dos2unix', output_file_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Failed to convert {output_file_name} to Unix format: {e}")

    def _run_script(self, output_file_path):
        """
        Submits a given script to the PBS queue.

        Parameters:
            output_file_path (str): The path to the script to be submitted.
        """
        try:
            subprocess.run(['qsub', output_file_path], check=True)
            print(f"Submitted PBS script: {output_file_path}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to submit PBS script: {e}")

    def run_download_script(self):
        """
        Submits the download script to the PBS queue.
        """
        output_file_path = f"downloadfits_S{self.sector}_o{self.orbit}.pbs"
        self._run_script(output_file_path)

    def run_photometry_script(self):
        """
        Submits the photometry script to the PBS queue.
        """
        output_file_path = f"runphot_S{self.sector}_o{self.orbit}.pbs"
        self._run_script(output_file_path)

    def run_light_curve_script(self):
        """
        Submits the light curve script to the PBS queue.
        """
        output_file_path = f"runlc_S{self.sector}_o{self.orbit}.pbs"
        self._run_script(output_file_path)

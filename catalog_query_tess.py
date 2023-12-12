
import optparse
import logging
import os
import urllib.request
import sys
import numpy as np
import pandas as pd
from astropy.io import fits
import astropy.wcs as wcs
import astropy.coordinates as coord
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord 
from astropy.time import Time
from astroquery.vizier import Vizier
import astroquery
from urllib.error import URLError
from astroquery.mast import Catalogs

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '..'))

import patools 
from patools.frame import Frame
from patools.field import Field
from patools.util.configurable import Configuration

logging.basicConfig(level=logging.DEBUG, filename='default_log_file.log', filemode='w', format='%(asctime)s - %(levelname)s - %(message)s')

def cmd_parse():
    p = optparse.OptionParser()
    p.add_option('--version', action="store_true", help="dislay the version ")
    p.add_option('--logfile', default='example.log', 
                 help='The destination of the log file, - outputs the logs to stdout')
    p.add_option('--pointings', help='the pointing file that corresponding to the four cameras')
    p.add_option('-O', '--outdir', help='The output directory for the catalog files')
    p.add_option('--debug', action='store_true', help='output debug level message')
    p.add_option('--maglim', default=12, help='the magnitude limit of the catalog')
    p.add_option('--width', default=0.3, help='the radius of the region to query')
    p.add_option('--orbitid', help='the orbitid of the catalogfile')
    p.add_option('--cadence', help='cadence number of one of the images (preferably start cadence')
    p.add_option('--numsquares', default=4, help='the number of subsquares to break the catalog image into')
    p.add_option('--sector', help='sector number')
    options, arguments = p.parse_args()
    
    return [options, arguments]


class Catalog(object):
    def __init__(self, catfile, ra, dec, width, colid=1, colra=2, coldec=3,
                 colmag=4, colx=5, coly=6, ra0=None, dec0=None, method=None, maglim=12):
        self.name = catfile
        self.ra = float(ra)
        self.dec = float(dec)
        if ra0 is None:
            self.ra0 = self.ra
        else:
            self.ra0 = float(ra0)
        if dec0 is None:
            self.dec0 = self.dec
        else:
            self.dec0 = float(dec0)
        self.width = float(width)
        self.colid = colid
        self.colra = colra
        self.coldec = coldec
        self.colmag = colmag
        self.colx = colx
        self.coly = coly
        self.method = method
        self.maglim = float(maglim)
        print("Config Catalog: name=%s, ra=%f, dec=%f, ra0=%f, dec0=%f, width=%f, method=%s, maglim=%f" % (self.name, self.ra, self.dec, self.ra0, self.dec0, self.width, self.method, self.maglim))
    def query(self, maglim=None):
        if maglim is None:
            maglim = self.maglim
        if self.method is None:
            # query the catalog with Vizier 
            v = Vizier(columns=['UCAC4', '_RAJ2000', 'e_RAJ2000', '_DEJ2000',
                                'e_DEJ2000', 'Jmag', 'Kmag', 'Vmag', 'pmRA', 'pmDE',
                                'e_pmRA', 'e_pmDE', 'rmag', 'imag'],
                       row_limit="unlimited", column_filters={"Vmag": "<%f" % (maglim)})
            # print self.ra,self.dec
            result = v.query_region(SkyCoord(ra=self.ra*u.degree,
                                                   dec=self.dec*u.degree,
                                                   frame='icrs'),
                                    width=self.width*u.degree, catalog=["UCAC4"])

            ids = result[0]['UCAC4']
            ras = result[0]['_RAJ2000'][:]
            decs = result[0]['_DEJ2000'][:]
            e_ra = result[0]['e_RAJ2000'][:]
            e_de = result[0]['e_DEJ2000'][:]
            jmag = result[0]['Jmag'][:]
            kmag = result[0]['Kmag'][:]
            vmag = result[0]['Vmag'][:]
            with open(self.name, mode='w') as fout:
                for i in range(len(ids)):
                    # ID[i]=''.join(ID[i].split('-'))
                    if jmag[i] == "--":
                        jmag[i] = np.nan
                    if kmag[i] == "--":
                        kmag[i] = np.nan
                    if vmag[i] == "--":
                        vmag[i] = np.nan
                    fout.write("%s %12.7f %12.7f %5.3f %d %d %5.3f %5.3f %5.3f\n" %
                               (ids[i], ras[i], decs[i], vmag[i], e_ra[i], e_de[i],
                                jmag[i], kmag[i], vmag[i]))
        elif self.method == 'TIC':
            # query the catalog from tic
            #TBD: reimplement a MAST TIC query
            c = SkyCoord(self.ra0*u.degree, self.dec0*u.degree, frame='icrs')
            print(self.maglim)
            catalog_data = Catalogs.query_criteria(coordinates=c, catalog="TIC", radius= "%s deg" % (self.width), Tmag=[-3,self.maglim]).to_pandas()
            #catalog_data = Catalogs.query_region("%f %f" % (self.ra0, self.dec0), catalog="TIC", radius= self.width, Tmag=[-3, self.maglim]).to_pandas()
            #print(catalog_data.columns) 
            epoch0 = 2015.5
            epoch_now = 2021.5
            catalog_data = catalog_data[catalog_data["plx"]>0]
            print(catalog_data.head())
            c = SkyCoord(np.array(catalog_data["RA_orig"])*u.degree, np.array(catalog_data["Dec_orig"])*u.degree, pm_ra_cosdec = np.array(catalog_data["pmRA"])*np.cos(np.array(catalog_data["Dec_orig"])/180.*np.pi)*u.mas/u.yr, pm_dec=np.array(catalog_data["pmDEC"])*u.mas/u.yr, distance=(1000./np.array(catalog_data["plx"]))*u.pc, obstime=Time('2015-06-30'))
            c1 = c.apply_space_motion(new_obstime=Time('2021-06-30'))

            #newra = np.array(catalog_data["RA_orig"]) + (2021.5-2015.5) * catalog_data["pmRA"] * 1.e-3 /3600. * np.cos(np.array(catalog_data["Dec_orig"])/180.*np.pi) 
            #newdec = np.array(catalog_data["Dec_orig"]) + (2021.5-2015.5) * catalog_data["pmDEC"] * 1.e-3/3600.
            catalog_data["ra"] = c1.ra
            catalog_data["dec"] = c1.dec
            catalog_data.to_csv(self.name, columns=["ID","ra","dec","Tmag","pmRA","pmDEC","Jmag","Hmag","Kmag"], header=False, index=False)
        elif self.method == 'Gaia':
            c = SkyCoord(self.ra0*u.degree, self.dec0*u.degree, frame='icrs')
            catalog_data = Catalogs.query_region(c, catalog="Gaia", radius= "%s deg" % (self.width), version=2).to_pandas()
            newra = np.array(catalog_data["ra"]) + (2021.5-2015.5) * catalog_data["pmRA"] * 1.e-3 /3600. * np.cos(np.array(catalog_data["dec"])/180.*np.pi) 
            newdec = np.array(catalog_data["dec"]) + (2021.5-2015.5) * catalog_data["pmDEC"] * 1.e-3/3600.
            catalog_data["ra"] = newra
            catalog_data["dec"] = newdec
            print(self.ra0, self.dec0, self.width)
            catalog_data.to_csv(self.name, columns=["source_id","ra","dec","phot_g_mean_mag","pmRA","pmdDEC","phot_bp_mean_mag","phot_gp_mean_mag","radius"], header=False, index=False)
        else:
            raise AttributeError("not implemented yet for catalog " \
                                  "query method %s" % self.method) 
        return
    
def process_catalogs(sector, orbit, numsquares, base_path):
    """
    Processes catalogs for a given sector orbit.

    Args:
        sector (int): TESS sector number.
        orbit (int): TESS orbit number.
        numsquares (int): Number of subsquares per CCD.
        base_path (str): Base path for data storage.

    Returns:
        str: Error messages, if any.
    """
    downloadfits(sector, orbit, numsquares, base_path)
    all_pointing_info = get_subsquare_centers_and_radius(numsquares, orbit, base_path)

    error_messages = ''
    for cam in range(1, 5):
        for ccd in range(1, 5):
            ccd_catalog, error_messages = process_cam_ccd(cam, ccd, all_pointing_info, orbit, base_path, error_messages)
    return error_messages

def process_cam_ccd(cam, ccd, all_pointing_info, orbit, base_path, error_messages):
    """
    Processes a single camera and CCD.

    Args:
        cam (int): Camera number.
        ccd (int): CCD number.
        all_pointing_info (dict): Pointing information.
        orbit (int): TESS orbit number.
        base_path (str): Base path for data storage.
        error_messages (str): Accumulated error messages.

    Returns:
        tuple: CCD catalog dataframe and updated error messages.
    """
    ccd_catalog = []
    pointing_info = all_pointing_info[f"pointing_info_cam{cam}_ccd{ccd}"]

    for i, row in pointing_info.iterrows():
        ra, dec, width = row['Center RA'], row['Center Dec'], row['Circle Radius']
        catfile = os.path.join(base_path, f'orbit-{orbit}/ffi/run/orbit{orbit}_header_cam{cam}ccd{ccd}_subsquare{i}.txt')
        ccd_catalog, error_messages = query_catalog(catfile, ra, dec, width, orbit, cam, ccd, i, error_messages)

    ccd_df = pd.concat(ccd_catalog, axis=0, ignore_index=False)
    ccd_df = ccd_df.drop_duplicates()

    save_catalog(ccd_df, orbit, cam, ccd, base_path)
    return ccd_catalog, error_messages

def query_catalog(catfile, ra, dec, width, orbit, cam, ccd, subsquare_index, error_messages):
    """
    Queries catalog for a given subsquare.

    Args:
        catfile (str): Catalog file path.
        ra (float): Right ascension.
        dec (float): Declination.
        width (float): Width of the search area.
        orbit (int): TESS orbit number.
        cam (int): Camera number.
        ccd (int): CCD number.
        subsquare_index (int): Index of the subsquare.
        error_messages (str): Accumulated error messages.

    Returns:
        tuple: Subsquare catalog dataframe and updated error messages.
    """
    retry_count = 1
    cat = Catalog(catfile, ra, dec, width, colid=1, colra=2, coldec=3, colmag=4, colx=5, coly=6, ra0=ra, dec0=dec, method="TIC", maglim=args.maglim)
    for attempt in range(retry_count + 1):
        try:
            cat.query()
            subsquare_catalog = cat
            tempdf = pd.read_csv(catfile, names=["ID", "ra", "dec", "Tmag", "pmRA", "pmDEC", "Jmag", "Hmag", "Kmag"], delimiter=',')
            os.remove(catfile)
            return tempdf, error_messages
        except RemoteServiceError as e:
            if attempt < retry_count:
                print(f"Error occurred: {e}. Retrying once for subsquare {subsquare_index} of cam{cam}-ccd{ccd}...")
                time.sleep(5)
            else:
                error_messages += f"No stars found in subsquare {subsquare_index} of cam{cam}-ccd{ccd}: {e}\n"
    return None, error_messages

def save_catalog(ccd_df, orbit, cam, ccd, base_path):
    """
    Saves the catalog dataframes to files.

    Args:
        ccd_df (DataFrame): The CCD catalog dataframe.
        orbit (int): TESS orbit number.
        cam (int): Camera number.
        ccd (int): CCD number.
        base_path (str): Base path for data storage.
    """
    bright_df = ccd_df[ccd_df['Tmag'] <= 10]
    bright_df.to_csv(os.path.join(base_path, f'orbit-{orbit}/ffi/run/catalog_{orbit}_{cam}_{ccd}_bright.txt'), header=False, index=False, sep=' ', na_rep='')

    full_df = ccd_df[ccd_df['Tmag'] <= 12]
    full_df.to_csv(os.path.join(base_path, f'orbit-{orbit}/ffi/run/catalog_{orbit}_{cam}_{ccd}_full.txt'), header=False, index=False, sep=' ', na_rep='')

def main():
    parser = argparse.ArgumentParser(description="Generate catalogs for TESS data.")
    parser.add_argument("--sector", type=int, required=True, help="TESS sector number.")
    parser.add_argument("--orbit", type=int, required=True, help="TESS orbit number.")
    parser.add_argument("--cadence", type=int, required=True, help="Start cadence number.")
    parser.add_argument("--numsquares", type=int, default=4, help="Number of subsquares per CCD.")
    parser.add_argument("--path", type=str, default="30daytemp", help="Base path for data storage.")

    args = parser.parse_args()
    base_path = os.path.expanduser(f'~/{args.path}')

    if not os.path.exists(base_path):
        os.makedirs(base_path, exist_ok=True)

    error_messages = process_catalogs(args.sector, args.orbit, args.numsquares, base_path)
    print(error_messages, end='')

if __name__ == '__main__':
    main()

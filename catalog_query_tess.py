#!/usr/bin/env python
# 
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT) 
# 
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 



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
from astroquery.exceptions import RemoteServiceError
import time

#from astroquery.mast import Conf 
#Conf.timeout = 3600
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
    p.add_option('--maglim', default=13.5, help='the magnitude limit of the catalog')
    p.add_option('--width', default=0.3, help='the radius of the region to query')
    p.add_option('--orbitid', help='the orbitid of the catalogfile')
    p.add_option('--cadence', help='cadence number of one of the images (preferably start cadence')
    p.add_option('--numsquares', help='the number of subsquares to break the catalog image into')
    p.add_option('--sector', help='sector number')
    p.add_option('--path', default='30daytemp', help='The path to the orbit directory')
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
    
class FitsDownloader:
    def __init__(self, sector, orbitid, cadence, path):
        self.sector = sector
        self.cadence = cadence
        self.path = path
        self.orbitid = orbitid
        self.orbit = 2 if self.orbitid % 2 == 0 else 1
        self.url_pattern = f'https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:HLSP/tica/s00{self.sector}/cam{{cam}}-ccd{{ccd}}/hlsp_tica_tess_ffi_s00{self.sector}-o{self.orbit}-00{self.cadence}-cam{{cam}}-ccd{{ccd}}_tess_v01_img.fits'

    def download_fits(self):
        for cam in range(1, 5):
            for ccd in range(1, 5):
                self._download_single_fits(cam, ccd)

    def _download_single_fits(self, cam, ccd):
        url = self.url_pattern.format(cam=cam, ccd=ccd)
        output_file = f'{self.path}/orbit-{self.orbitid}/ffi/run/cam{cam}-ccd{ccd}.fits'
        try:
            urllib.request.urlretrieve(url, output_file)
        except URLError as e:
            logging.error(f'Error downloading {url}: {e}')
 
class SubsquareInfoGenerator:
    def __init__(self, num_squares, path, orbitid):
        self.num_squares = num_squares
        self.path = path
        self.orbitid = orbitid

    def get_subsquare_centers_and_radius(self):
        all_pointing_info = {}
        for cam in range(1, 5):
            for ccd in range(1, 5):
                filename = f"{self.path}/orbit-{self.orbitid}/ffi/run/cam{cam}-ccd{ccd}.fits"
                with fits.open(filename) as file:
                    data = file[0].data
                    world = WCS(file[0].header)

                subsquare_info = self._calculate_subsquare_info(data, world)
                df_name = f"pointing_info_cam{cam}_ccd{ccd}"
                all_pointing_info[df_name] = pd.DataFrame(subsquare_info)
        return all_pointing_info

    def _calculate_subsquare_info(self, data, world):
        image_shape = data.shape
        squares = (image_shape[0] // self.num_squares, image_shape[1] // self.num_squares)

        subsquare_info = []
        for i in range(self.num_squares):
            for j in range(self.num_squares):
                x_min, x_max = i*squares[0], (i+1)*squares[0]
                y_min, y_max = j*squares[1], (j+1)*squares[1]

                coords = world.pixel_to_world_values([x_min, x_max], [y_min, y_max])
                center_ra, center_dec = np.median(coords[0]), np.median(coords[1])

                x_sep_deg = np.abs(world.pixel_to_world(x_min, y_min).separation(world.pixel_to_world(x_max, y_min)).deg)
                y_sep_deg = np.abs(world.pixel_to_world(x_min, y_min).separation(world.pixel_to_world(x_min, y_max)).deg)

                square_cross_deg = np.sqrt(x_sep_deg**2 + y_sep_deg**2)
                circle_radius = 0.5 * square_cross_deg * 1.2

                subsquare_info.append({'Center RA': center_ra, 'Center Dec': center_dec, 'Circle Radius': circle_radius})
        return subsquare_info
    
class SectorProcessor:
    def __init__(self, all_pointing_info, path, orbitid, maglim):
        self.all_pointing_info = all_pointing_info
        self.path = path
        self.orbitid = orbitid
        self.maglim = maglim

    def process_sector(self):
        error_messages = ''

        for cam in range(1, 5):
            for ccd in range(1, 5):
                ccd_catalog, error_messages = self._process_ccd(cam, ccd, error_messages)

                ccd_df = pd.concat(ccd_catalog, axis=0, ignore_index=False)
                ccd_df = ccd_df.drop_duplicates()

                bright_df = ccd_df[ccd_df['Tmag'] <= 10]
                bright_df.to_csv(f'{self.path}/orbit-{self.orbitid}/ffi/run/catalog_{self.orbitid}_{cam}_{ccd}_bright.txt', header=False, index=False, sep=' ', na_rep='')

                full_df = ccd_df[ccd_df['Tmag'] <= 12]
                full_df.to_csv(f'{self.path}/orbit-{self.orbitid}/ffi/run/catalog_{self.orbitid}_{cam}_{ccd}_full.txt', header=False, index=False, sep=' ', na_rep='')

        return error_messages

    def _process_ccd(self, cam, ccd, error_messages):
        ccd_catalog = []
        catfiles_to_remove = []  # List to store catfile paths
        pointing_info = self.all_pointing_info[f"pointing_info_cam{cam}_ccd{ccd}"]

        for i, row in pointing_info.iterrows():
            ra = row['Center RA']
            dec = row['Center Dec']
            width = row['Circle Radius']
            catfile = f"{self.path}/orbit-{self.orbitid}/ffi/run/orbit{self.orbitid}_header_cam{cam}ccd{ccd}_subsquare{i}.txt"
            catfiles_to_remove.append(catfile)  # Add catfile path to the list

            cat = Catalog(catfile, ra, dec, width, colid=1, colra=2, coldec=3, colmag=4, colx=5, coly=6, ra0=ra, dec0=dec, method="TIC", maglim=self.maglim)
            ccd_catalog, error_messages = self._query_catalog(cat, catfile, i, cam, ccd, ccd_catalog, error_messages)

        # Remove subsquare files after all have been processed
        for catfile in catfiles_to_remove:
            os.remove(catfile)

        return ccd_catalog, error_messages


    def _query_catalog(self, cat, catfile, i, cam, ccd, ccd_catalog, error_messages):
        retry_count = 1
        for attempt in range(retry_count + 1):
            try:
                cat.query()
                tempdf = pd.read_csv(catfile, names=["ID", "ra", "dec", "Tmag", "pmRA", "pmDEC", "Jmag", "Hmag", "Kmag"], delimiter=',')
                ccd_catalog.append(tempdf)
                break
            except RemoteServiceError as e:
                if attempt < retry_count:
                    print(f"Error occurred: {e}. Retrying once for subsquare {i} of cam{cam}-ccd{ccd}...")
                    time.sleep(5)
                else:
                    error_messages += f"No stars found in subsquare {i} of cam{cam}-ccd{ccd}: {e}\n"
        return ccd_catalog, error_messages



if __name__ == '__main__':
    options, args = cmd_parse()
    if options.logfile:
        logging.basicConfig(level=logging.DEBUG, filename=options.logfile, filemode='w', format='%(asctime)s - %(levelname)s - %(message)s')
    patools.setup_logging(debug=options.debug, filename=options.logfile, multi=False)
    logger = logging.getLogger(__name__)

    #pointings = options.pointings
    #out_dir = options.outdir
    path = options.path

    maglim = float(options.maglim)
    orbitid = int(options.orbitid)
    numsquares = int(options.numsquares)
    cadence = int(options.cadence)
    sector = int(options.sector)
  
    
    downloader = FitsDownloader(sector, orbitid, cadence, path)
    downloader.download_fits()
    
    generator = SubsquareInfoGenerator(numsquares, path, orbitid)
    all_pointing_info = generator.get_subsquare_centers_and_radius()
    
    processor = SectorProcessor(all_pointing_info, path, orbitid, maglim)
    error_messages = processor.process_sector() 

    print(error_messages)

if __name__ == "__main__":
    main()   

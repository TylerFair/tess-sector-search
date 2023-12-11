# TESS Sector Search and Analysis Tools

This project provides a suite of tools for analyzing TESS (Transiting Exoplanet Survey Satellite) sectors, handling the extraction of cadences, spacecraft and camera pointings, directory creation for sector searches, and catalog generation. It also includes functionality to generate and execute PBS scripts for high-performance computing (HPC) environments.

- Determine start and end cadences for each orbit in a sector.
- Retrieve and utilize spacecraft and camera pointing data.
- Create necessary directories for a sector search.
- Generate star catalogs for specific camera/CCD images.
- Generate and execute PBS scripts for data processing on HPC systems.

# Modules

- tess_sector_search.py: Main script for sector analysis, including cadence handling, pointing information retrieval, directory setup, and catalog creation.
- pbs_script_generator.py: Module for generating and executing PBS scripts for HPC usage.

# Usage

# Command-Line Execution

Run the tess_sector_search.py script with the desired sector and orbit numbers. Optional arguments include specifying a path, generating catalogs, and controlling the execution of download, photometry, and light curve scripts.

- python tess_sector_search.py sector orbit [--path=default] [--catalog=False] [--download=False] [--photometry=False] [--lc=False]

# Module Import
The script can also be used as a module, allowing the execution of specific functions as needed.

    from tess_sector_search import SectorSetup, PBSRun

Initialize the SectorSetup with desired parameters

    tess_processor = SectorSetup(sector, orbit, path)
    tess_processor.catalog_required = catalog

Create all necessary files for the sector

    tess_processor.create_sector_files()

Generate PBS scripts for data processing

    tess_processor.generate_pbs_scripts()

Initialize PBSRun with the created PBSScript instance

    pbs_run = PBSRun(tess_processor.pbs_script)

Execute specific scripts

    if download:
        pbs_run.download()

    if photometry:
        pbs_run.photometry()

    if lc:
        pbs_run.lc()
        
# Dependencies

- Python 3.x
- Requests
- Pandas
- Subprocess

Ensure all dependencies are installed to use the scripts effectively.

# Authors

Tyler R. Fairnington, building off of public QLP by Chelsea X. Huang. 

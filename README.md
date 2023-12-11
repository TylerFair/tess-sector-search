This code builds off of the public QLP FFI code by C. X. Huang.

tess_sector_search.py can be used to find the start and end cadences of each orbit in a sector, save and use the spacecraft and camera pointings and create directories for a sector search.
It can also generate catalogs for all stars in a given cam/ccd image. 

pbs_script_generator.py is a module invoked in tess_sector_search.py and can be used to generate PBS scripts for HPC use, as well as executing the same scripts. 

### Example Usage ###
if executing on shell through main: 
python tess_sector_search.py <sector> <orbit> <--path=default> <--catalog=False> <--download=False> <--photometry=False> <--lc=False>

else:
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

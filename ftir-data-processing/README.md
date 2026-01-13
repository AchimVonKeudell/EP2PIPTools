
# FTIR-data-processing

## Synopsis

This project can simulate IR spectra and fit measured spectra, to characterise the corresponding parameters. See sub-readme file
[README.md](/ftir_data_processing/vibrationalspectrum/README.md).

The project is used in the following publications by our group:
- [Yu et al 2024 J. Phys. D: Appl. Phys. 57 245203](https://doi.org/10.1088/1361-6463/ad2ef6)
- [Vervloedt and von Keudell 2024 Plasma Sources Sci. Technol. 33 045005](https://doi.org/10.1088/1361-6595/ad38d6)
- [Vervloedt and von Keudell Plasma Chem Plasma Process 45, 1551–1565 (2025)](https://doi.org/10.1007/s11090-025-10584-x) 
- Seferoglu et al (to be written)


## Installation 
1. install python following your preferred method. Personally, I made use of [miniforger](https://conda-forge.org/miniforge/), 
which allows you to use conda to control your python environment for different projects.  
2. Install the required packages by


    pip install -r ./requirments.txt

3. Configure SQLite database for the molecule of your choosing, see [First Steps.pdf](docs/First%20Steps.pdf)


## Workings:
The project can be operated using command prompt (Windows), or similar methods when using another operating software. 
Spectra can be simulated and fitted after configuring the Python environment and installing the required packages, see 
`Installation` section. 

### Simulate spectra
Simulate spectra: see scripts in [simulate_spectra](ftir_data_processing\simulate_spectra). Mulitple scripts are written 
to simulate the vibrational spectrum of different molecules. Please adjust the parameter `hitran_database_file` according 
to the location of the SQlite database on your computer.

### Fit spectra
Call [run.py](run.py) in command prompt to fit the data following a configured `.yaml` file. The template for these files 
can be found under [templates](templates).

Try calling 

    python run.py --help

It provides information on which functions are available to the user. 


## Database Configuration 
The project makes use of an SQLite database structure to efficiently store and index data on molecules, isotopologues, absorption structures and data
on line-by-line parameters. These parameters are stored files with a `.sqlite` extension. This database is filled using the
HITRAN or HITEMP database. This database is configured as follows

This HITRAN input file directory should contain `molparam.txt` (contains data on molecules and 
isotopologues), line parameters files `*.par` (from this folder will be imported as long as they start with 
the 2-digit molecular code defined by HITRAN, followed by an underscore and some additional text), and `parsum.dat`
(contains the partition sums of the various isotopologues to be included).

An example of a filename would be `02_2000-2125_HITEMP2010.par`, which contains HITEMP data of the CO2 molecule 
(with molecule id 2).

Running [generate_database.py](generate_database.py) will read all available HITRAN/HITEMP data from 
the specified HITRAN directory and store it in a sqlite database, which is used by this script to select certain 
absorption lines to do simulations.
 
## API Reference 
 
No API is available yet. 


## Author
Steijn Vervloedt (Ruhr Universty Bochum, 2022-2026)

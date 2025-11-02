# SDL-HST-Assignment-Daragh-Haughton ReadMe

## HST Photometry pipeline for NGC 1261

This script performs photometry on Hubble space telescope imagery in 2 filters (F555W and F336W) of the clobular cluster NGC 1261. It removes cosmic rays, identifies stellar sources, performs photometry on the sources to extract the magnitudes, and produces a Hertzprung Russel (HR) Diagram (more specifically, a colour (F336 - F555) vs magnitude (F555) diagram.

## Repository Structure
SDL-HST-Assignment-Daragh-Haughton/
├── data/
│ ├── F336W/
│ │ ├── file1.fits
│ │ ├── file2.fits
│ │ └── file3.fits
│ └── F555W/
│ ├── file1.fits
│ ├── file2.fits
│ └── file3.fits
├── myscript.py
└── README.md

## Requirements
- astropy
- numpy
- matplotlib
- scipy
- pandas

  These can be installed with:
  pip install astropy numpy matplotlib scipy pandas

  ## Usage
  This code can be ran by:
  1. Simplest case:
      Defaults (percentile=93, saves CSV, no plot)
      $ python myscript.py

  2. With a custom data folder
      $ python myscript.py /path/to/data
      $ python myscript.py ../hst_data

  3. Custom data folder AND percentile value
      $ python myscript.py data 95
      $ python myscript.py /path/to/data 99
      $ python myscript.py data 90

  4. No CSV saving
      $ python myscript.py nocsv
      $ python myscript.py data nocsv
      $ python myscript.py data 95 nocsv

  5. Generate overlay plot
      $ python myscript.py plot
      $ python myscript.py data plot
      $ python myscript.py data 95 plot

  6. Combine flags (the order does not matter for the flags)
      $ python myscript.py data 95 nocsv plot
      $ python myscript.py data plot nocsv
      $ python myscript.py nocsv plot
      $ python myscript.py data 99 plot nocsv

##  Outputs
The script will generate the following:
1. hrdiagramofNGC1261.png - A HR diagram
2. SuccessfulCandPos.png - OPTIONAL - A plot showing a transformed image and the coordinates of the final candidates overlaid - Only produced if plot flag used
3. 5 CSV Files - OPTIONAL:
   i) Accepted (Quality) - filter: F336W
   ii) Accepted (Quality) - filter: F555W
   iii) Rejected (Quality) - filter: F336W
   iv) Rejected (Quality) - filter: F555W
   v) Finalcatalogue.csv

   The accepeted CSV files contain the sources accepeted based on PSF filters being passed
   The rejected CSV files contain the sources rejected based on PSF filters not being passed, as    well as the reason for their not passing
   The final candidates CSV file contains the Accepeted sources which also passed the cross         filter position match check.

   These files are only produced if the 'nocsv' flag is not used.


## Overview of process
1. File read in: FITS files are loaded from relevant filter directories
2. Processing: Cosmic rays are removed via median combination of the files, the background is subtracted
3. Detection: The sources over the threshold variable are identified
4. Photometry: Circular masking is performed, to obtain the magnitude of the source
5. Candidate Selection: Based on SigNoise ratio, ellipticity, spatial extent and edge filtration, successful candidates from each filter are selected.
6. Cross matching: Match sources between filters based on their positions in both filters
7. Visualisation: HR diagram is generated for the successful candidates


## Notes:
# NB: The percentile threshold value significantly impacts runtime:
At default value (93), which produces the best plot and best number of sources, the current runtime is 29 minutes! Many attempts were made to reduce this, although these resulted in poorer diagrams and lower overall source populations. 

95 is a good threshold value for verifying the code works and produces a decent diagram, although this also has a long runtime of around 12 minutes. 

# GenAI was used during the making of this, please see comments in code to see where, and what it was used for.
Specific model used: OpenAI's ChatGPT-4-Turbo. 

# I have inlcuded a folder of previously generated results, with default parameters used

## Name: Daragh Haughton
## Student No. 25237942

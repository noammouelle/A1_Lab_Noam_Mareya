INSTRUCTIONS FOR RUNNING OUR CODE



1. CONTENTS ###################################################################
 
Both versions analyse the image contained in "mosaic.fits".

This folder contains two working versions of our code to obtain a galaxy number
count.

The version 1 (variable aperture shape) is composed of 3 files:
    1. "A1_tools.py", which is a package containing all the routines needed to 
    analyse the image.
    
    2. "cleaning_data.py", which removes unwanted features from the image and
    creates a file, "modified_img.txt", containing the improved data.
    
    2. "img_analysis.py", which detects sources in the image and produces a txt
    file, "run2_calibrated_stats.txt", containing the calibrated source list.
    
    3. "analysing_calibrated_data.py", which analyses the file created by the 
    first program.
    
    4. "run1_calibrated_stats.txt", which contains data already collected to
    avoid having to fully run the programs.
    
The version 2 (circular aperture) is composed of 3 files:
    1. "aperture_tools.py", which contains various functions required for 
    circular aperture photometry.
    
    2. "variable_aperture.py", which contains the function to run variable
    circular radius aperture photometry. It tests the function on a couple
    small patches and plots these for easy viewing.
    
    3. "number-count-plot.py" contains the analysis needed to produce a 
    number count plot.
    
2. INSTRUCTIONS ###############################################################

For version 1:
    a. Files must be executed in the following order: "img_analysis.py", 
    "analysing_calibrated_data.py". 
    b. Files were divided in Spyder cells which run one after the other.
    
For version 2:
    a. First execute "aperture_tools.py". This is split up into cells
    as well which can be run one after the other.
    b. To test out the algorithm, run cells one after the other in the 
    "variable_aperture.py" file.
    c. A number count plot can be produced from the pre-set catalogue data
    "catalogue_final_3465_3_9.txt", which is also in this GitHub repo.
    
    

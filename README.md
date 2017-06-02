# plastic-networks Readme
====================

This repository will store code associated with our plastic network model for tidewater glaciers.  

As of 1 June 2017: 

-Major model functions are housed in "plastic_utilities.py" and the network approach is handled in a separate script for each glacier simulated.  

-Utilities under development are in plastic_utilities_v2.py" and "plastic-glacier-class.py".  We are working to change how the model runs, so that it will be based on a class that contains model functions and can read in and store data for each glacier.

-A working example, covering three branches of Columbia Glacier and making an aerial difference plot, is in "columbia_aerialplot-example.py".  The bed topography and thickness data necessary to initialize it was provided by Bob McNabb (then of UAF, now Oslo) and is stored as "Data/McNabb_11J249_S1.nc".  The coordinates of the hand-selected flowlines are in "Data/mainline.txt", "easttribline.txt", and "westtribline.txt"

-Additional functions used to process the Alaska tidewater glacier data provided by Bob McNabb, Christian Kienholz, Matthias Huss, and Shad O'Neel can be found in "alaska_processing.py", but we do not yet have an example up using that code.  We are archiving it because these are the functions we used for Hubbard Glacier, as described in our Frontiers manuscript.

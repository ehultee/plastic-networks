# plastic-networks Readme
====================

This repository will store code associated with our plastic network model for calving glaciers.  

As of 21 June 2018:
-All functions associated with the network approach are implemented in a set of classes in "flowline_class_hierarchy.py".

    -Class "Ice" stores useful constants e.g. density of ice, characteristic lengths for nondimensionalization
    
    -Class "Branch" stores some functions applicable to portions of glacier networks that do not necessarily reach the terminus.  Branch allows networks to be initialized with partial flowlines that are then connected by KDTree search.
    
    -Class "Flowline" stores many of the functions that must be computed separately for each flowline in a network, including the creation of snapshot profiles and the calculation of net volume change between successive profiles
    
    -Class "Network" stores functions that apply to collections of interacting flowlines.  In particular, time evolution is always done from the network class, even if the network consists of a single flowline.

-"Greenland-flux_projections.py" stores code used to generate and plot simulations associated with EHU's dissertation and presentation at POLAR2018.

-"Greenland-summary_plotting.py" stores function to quickly read in model output from saved simulations (saved as pickle files with Network.save_network) and display it together for multi-scenario analysis.




As of 1 June 2017: 

-Major model functions are housed in "plastic_utilities.py" and the network approach is handled in a separate script for each glacier simulated.  

-Utilities under development are in "plastic_utilities_v2.py" and "plastic-glacier-class.py".  We are working to change how the model runs, so that it will be based on a class that contains model functions and can read in and store data for each glacier.

-A working example, covering three branches of Columbia Glacier and making an aerial difference plot, is in "columbia_aerialplot-example.py".  The bed topography and thickness data necessary to initialize it was provided by Bob McNabb and is stored as "Data/McNabb_11J249_S1.nc".  The coordinates of the hand-selected flowlines are in "Data/mainline.txt", "easttribline.txt", and "westtribline.txt"

-Additional functions used to process the Alaska tidewater glacier data provided by Bob McNabb, Christian Kienholz, Matthias Huss, and Shad O'Neel can be found in "alaska_processing.py", but we do not yet have an example up using that code.  We are archiving it because these are the functions we used for Hubbard Glacier, as described in our manuscript:
    L. Ultee & J. N. Bassis (2017), "A plastic network approach to model calving glacier advance and retreat", Frontiers in Earth Sciences.

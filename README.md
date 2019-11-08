# plastic-networks Readme
====================

Welcome!  This is a working repository associated with the SERMeQ model for marine-terminating glaciers.  You'll find core model code in the SERMeQ folder; with working examples separated by their science application.

## The model
SERMeQ - the Simple Estimator of Retreat Magnitude and Ice Flux (Q) - is a vertically-integrated, width-averaged model that tracks the loss of glacier mass to iceberg calving.
It self-consistently determines a rate of terminus advance/retreat, and associated upstream thickening/thinning, based on surface mass balance and glacier geometry.

The repository is called plastic-networks to highlight two key features: the glacier rheology is similar to the _plastic_ approximation of Nye, such that ice can break in a boundary layer along the glacier bed and at the front,
and the model formulation handles interacting _networks_ of glacier flowlines.
For more on the model physics, please see

	* Ultee and Bassis (2016). The future is Nye: an extension of the perfect plastic approximation to tidewater glaciers. _Journal of Glaciology_ 62(236), 1143–1152. doi:10.1017/jog.2016.108. [PDF](http://ehultee.github.io/files/jog-2016_108.pdf)
	* Ultee and Bassis (2017). A plastic network approach to model calving glacier advance and retreat. _Frontiers in Earth Sciences_ 5(24). doi:10.3389/feart.2017.00024. [PDF](http://ehultee.github.io/files/feart-2017_00024.pdf)
	* Bassis and Ultee (2019). A thin film viscoplastic model for calving glaciers: an upper bound on calving retreat. _Journal of Geophysical Research: Earth Surface_. doi:10.1029/2019JF005160. [PDF](http://ehultee.github.io/files/BassisUltee-thin_film_viscoplastic-2019.pdf)


## Authors
*Lizz Ultee* (repo owner, "EHU" author attribution in code) developed this model at the University of Michigan in conversation with *Jeremy Bassis* ("JNB" attribution in code).


## Major updates
As of 3 Sept 2019:
* Core model code is in the SERMeQ subfolder, which has an \__init\__.py file so that Python 2.x treats it as a package to import.  
All scripts in subfolders Hindcasts-MS, Projections, and Visualization have an updated import command (using sys.path) that should preserve all functionality.

As of 29 Aug 2019:
* Core model code has been moved into SERMeQ subfolder for import.  If you are trying to run a script in another subfolder and you run into errors importing modules, this is probably why.


As of 21 June 2018:
* All functions associated with the network approach are implemented in a set of classes in "flowline_class_hierarchy.py".

    * Class "Ice" stores useful constants e.g. density of ice, characteristic lengths for nondimensionalization
    
    * Class "Branch" stores some functions applicable to portions of glacier networks that do not necessarily reach the terminus.  Branch allows networks to be initialized with partial flowlines that are then connected by KDTree search.
    
    * Class "Flowline" stores many of the functions that must be computed separately for each flowline in a network, including the creation of snapshot profiles and the calculation of net volume change between successive profiles
    
    * Class "Network" stores functions that apply to collections of interacting flowlines.  In particular, time evolution is always done from the network class, even if the network consists of a single flowline.

* "Greenland-flux_projections.py" stores code used to generate and plot simulations associated with EHU's dissertation and presentation at POLAR2018.

* "Greenland-summary\_plotting.py" stores function to quickly read in model output from saved simulations (saved as pickle files with Network.save_network) and display it together for multi-scenario analysis.




As of 1 June 2017: 

* Major model functions are housed in "plastic_utilities.py" and the network approach is handled in a separate script for each glacier simulated.  

* Utilities under development are in "plastic\_utilities_v2.py" and "plastic-glacier-class.py".  We are working to change how the model runs, so that it will be based on a class that contains model functions and can read in and store data for each glacier.

* A working example, covering three branches of Columbia Glacier and making an aerial difference plot, is in "columbia\_aerialplot-example.py".  The bed topography and thickness data necessary to initialize it was provided by Bob McNabb and is stored as "Data/McNabb\_11J249_S1.nc".  The coordinates of the hand-selected flowlines are in "Data/mainline.txt", "easttribline.txt", and "westtribline.txt"

* Additional functions used to process the Alaska tidewater glacier data provided by Bob McNabb, Christian Kienholz, Matthias Huss, and Shad O'Neel can be found in "alaska_processing.py", but we do not yet have an example up using that code.  We are archiving it because these are the functions we used for Hubbard Glacier, as described in our manuscript:
    L. Ultee & J. N. Bassis (2017), "A plastic network approach to model calving glacier advance and retreat", Frontiers in Earth Sciences.

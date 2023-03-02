# SuperNEMO timestamps analysis
Code for and results from anaysis of timestamps from SuperNEMO tracker
This code can be found also found in folder /sps/nemo/scratch/amendl/SuperNEMO-timestamps-analysis on CCLyon.
If you have any questions or remarks feel free to contact me at [adam.mendl@cvut.cz](mailto:adam.mendl@cvut.cz)
## .root files
 - **correlation_plots.root** - correlations between different ways how to obtain z value
 - **output.root** - Please not that event is used only when all timestamps are valid!
 - **without_propagation_time_correction.root** - z distribution of hits (absolute value of distance from the middle) without correction for variable propagation time
 - **with_propagation_time_correction.root** - z distribution of hits (absolute value of distance from the middle of the detector) abs(r5-r6)/(r5+r6-2r0)
## Code
All code relies on RED format. It is build by placing it in folder with RED CMakeLists.txt file a callign ``cmake CMakeLists.txt`` and ``make``.
 - **analysis.cxx** - produces ``output.root``
 - **correlations.cxx** - produces ``correlation_plots.root``

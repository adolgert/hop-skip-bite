
# Create the streets and blocks and have a look at them before intersection.
houses.txt streetsxy.txt streetspts.txt nocross.pdf: generate_landscape.R try.R
	R --no-save --args 1000 30 streets < generate_landscape.R

# Slow. Find intersections among streets and blocks.
landscape.pdf crossings.txt: houses.txt streetsxy.txt try.R
	R --no-save --args houses.txt streets < calculate_crossings.R

# Run the sir on the landscape.
sir.h5: try.R sir_run.R
	R --no-save --args houses.txt crossings.txt sir.h5 2 < sir_run.R

# Run a verhulst on the landscape.
verhulst.h5: verhulst_run.R try.R sir_run.R houses.txt crossings.txt
	R --no-save --args houses.txt crossings.txt verhulst.h5 1 -1 < verhulst_run.R

# Run verhulst on a single house.
verhulst_single.h5: verhulst_run.R try.R sir_run.R one_house.txt one_house_cross.txt
	R --no-save --args one_house.txt one_house_cross.txt verhulst_single.h5 20 5.0 < verhulst_run.R

# Get a plot of verhulst on a single house.
bug_trajectory.png: try.R bugplot.R verhulst_single.h5
	R --no-save --args verhulst_single.h5 bug_trajectory.png < bugplot.R

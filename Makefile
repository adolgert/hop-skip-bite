
houses.txt streetsxy.txt streetspts.txt nocross.pdf: generate_landscape.R try.R
	R --no-save --args 100 20 streets < generate_landscape.R

landscape.pdf crossings.txt: houses.txt streetsxy.txt try.R
	R --no-save --args houses.txt streets < calculate_crossings.R

sir.h5: try.R sir_run.R
	R --no-save --args houses.txt crossings.txt sir.h5 2 < sir_run.R

verhulst.h5: try.R sir_run.R houses.txt crossings.txt
	R --no-save --args houses.txt crossings.txt verhulst.h5 2 < verhulst_run.R

verhulst_single.h5:  try.R sir_run.R one_house.txt one_house_cross.txt
	R --no-save --args one_house.txt one_house_cross.txt verhulst_single.h5 2 < verhulst_run.R

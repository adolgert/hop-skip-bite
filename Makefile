
houses.txt streetsxy.txt streetspts.txt nocross.pdf: generate_landscape.R try.R
	R --no-save --args 100 20 streets < generate_landscape.R

landscape.pdf crossings.txt: houses.txt streetsxy.txt try.R
	R --no-save --args houses.txt streets < calculate_crossings.R

z.h5: try.R sir_run.R
	R --no-save --args houses.txt crossings.txt z.h5 2 < sir_run.R

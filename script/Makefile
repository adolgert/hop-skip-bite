REXEC=R --no-save
BOOST=/home/ajd27/Documents/boost_1_57_0

.PHONY: hopskip

hopskip:
	{ \
	cd .. ;\
	R CMD build hopskip ;\
	for i in $$(seq 1 500); do echo; done ;\
	R CMD INSTALL --configure-args="--with-boost=$(BOOST)" hopskip_1.0.tar.gz ;\
	}

# Create the streets and blocks and have a look at them before intersection.
10/houses.txt streetsxy.txt streetspts.txt nocross.pdf: generate_landscape.R try.R
	$(REXEC) --args 20 5 streets 33334 < generate_landscape.R
	$(REXEC) --args houses.txt streets < calculate_crossings.R
	DDIR=10 && mkdir -p $$DDIR  && cp *.txt $$DDIR && cp landscape.pdf $$DDIR

# Create the streets and blocks and have a look at them before intersection.
100/houses.txt streetsxy.txt streetspts.txt nocross.pdf: generate_landscape.R try.R
	$(REXEC) --args 100 10 streets 33333 < generate_landscape.R
	$(REXEC) --args houses.txt streets < calculate_crossings.R
	DDIR=100 && mkdir -p $$DDIR  && cp *.txt $$DDIR && cp landscape.pdf $$DDIR

# Create the streets and blocks and have a look at them before intersection.
1000/houses.txt streetsxy.txt streetspts.txt nocross.pdf: generate_landscape.R try.R
	$(REXEC) --args 1000 30 streets 33333 < generate_landscape.R
	$(REXEC) --args houses.txt streets < calculate_crossings.R
	DDIR=1000 && mkdir -p $$DDIR  && cp *.txt $$DDIR && cp landscape.pdf $$DDIR

# Run the sir on the landscape.
sir.h5: try.R sir_run.R
	$(REXEC) --args houses.txt crossings.txt sir.h5 2 < sir_run.R

# Run a verhulst on the landscape.
verhulst.h5: verhulst_run.R try.R sir_run.R houses.txt crossings.txt
	$(REXEC) --args 10/houses.txt 10/crossings.txt verhulst.h5 1 10 < verhulst_run.R

# Run verhulst on a single house.
verhulst_single.h5: verhulst_run.R try.R sir_run.R one_house.txt one_house_cross.txt
	$(REXEC) --args one_house.txt one_house_cross.txt verhulst_single.h5 20 5.0 < verhulst_run.R

# Get a plot of verhulst on a single house.
bug_trajectory.png: try.R bugplot.R verhulst_single.h5
	$(REXEC) --args verhulst_single.h5 bug_trajectory.png < bugplot.R

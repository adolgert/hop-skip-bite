# Requires:
#   BOOST libraries: boost.org
#   Semi-Markov library: https://github.com/afidd/Semi-Markov
#

BOOST=/home/ajd27/Documents/boost_1_57_0
# Different Boost installations have different suffixes.
# If there is no suffix, use "BOOSTVARIANT=".
BOOSTVARIANT=
SEMIMARKOV=/usr/local/include/semimarkov-0.1

# Use local.mk to change defaults for the variables above.
-include local.mk
# -DSMVHIDELOG -pg
OPT=-std=c++11 -DHAVE_CONFIG_H -fPIC -g -O2 -DSMVHIDELOG
INCLUDES=-I$(SEMIMARKOV) -I. -I$(BOOST)/include -I$(HDF5)/include
LIBS=-L$(BOOST)/lib -L$(HDF5)/lib64  \
	-lboost_unit_test_framework$(BOOSTVARIANT) \
	-lboost_log_setup$(BOOSTVARIANT) -lboost_log$(BOOSTVARIANT) \
	-lboost_chrono$(BOOSTVARIANT) -lboost_thread$(BOOSTVARIANT) \
	-lboost_date_time$(BOOSTVARIANT) -lboost_filesystem$(BOOSTVARIANT) \
	-lboost_program_options$(BOOSTVARIANT) -lboost_random$(BOOSTVARIANT) \
	-lboost_system$(BOOSTVARIANT) -lboost_timer$(BOOSTVARIANT) \
	-lpthread

ALLFILES=simple_hazard.cpp simple_hazard.hpp main.cpp spatial_process.hpp Makefile

all: simple_hazard hsb_version.hpp
.PHONY: all

simple_hazard: simple_hazard.o main.o
	$(CXX) $(OPT) -o simple_hazard simple_hazard.o main.o $(LIBS)

main.o: main.cpp simple_hazard.hpp
	$(CXX) -fPIC $(INCLUDES) $(OPT) main.cpp -c -o main.o

simple_hazard.o: simple_hazard.cpp simple_hazard.hpp spatial_process.hpp
	$(CXX) $(OPT) -fPIC $(INCLUDES) -c simple_hazard.cpp -o simple_hazard.o

hsb_version.hpp:
	python getgit.py

clean:
	rm -f *.o simple_hazard hsb_version.hpp

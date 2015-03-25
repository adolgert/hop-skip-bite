constexpr char VERSION[]=R"(git@github.com:adolgert/hop-skip-bite.git:31aa55ebb1ba3ae581dbd16fdfd722906417fc39)";

constexpr char CFG[]=R"(# Requires:
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

all: simple_hazard
.PHONY: all

simple_hazard: simple_hazard.o main.o
    $(CXX) $(OPT) -o $@ $< $(LIBS)

main.o: main.cpp simple_hazard.hpp
    $(CXX) main.cpp -DHAVE_CONFIG_H -std=c++11 -fPIC $(INCLUDES) $(OPT) \
    -c -o main.o

simple_hazard.o: simple_hazard.cpp simple_hazard.hpp spatial_process.hpp
    $(CXX) $(OPT) -fPIC -o simple_hazard.o $(LIBS)

hsb_version.hpp: Makefile
    python getgit.py

clean:
    rm -f *.o sirexp together rider individual ensemble_sum sirdemo_version.hpp
)";

constexpr char COMPILETIME[]=R"(2015-03-18T13:44:27.977963)";


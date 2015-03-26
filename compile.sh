R CMD build hopskip
for i in {1..500}; do echo; done
R CMD INSTALL --configure-args="--with-boost=/home/ajd27/Documents/boost_1_57_0" hopskip_1.0.tar.gz

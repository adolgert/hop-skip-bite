
source("try.R")
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
houses <- args[1]
crossings <- args[2]
output <- args[3]
run.cnt <- as.integer(args[4])
if (length(args)>4) {
  max.time <- as.numeric(args[5])
} else {
  max.time <- -1
}
verhulst.generate(housesfile=houses, crossingsfile=crossings,
  outfile=output, run_cnt=run.cnt, max_time=max.time)

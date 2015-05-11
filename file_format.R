library(rhdf5)

create_file <- function(filename) {
  if (file.exists(filename)) {
    file.remove(filename)
  }
  h5createFile(filename)
  h5createGroup(filename, "trajectory")
  h5createGroup(filename, "other")
}


write_single_dataset <- function(outfile, trajectory_group, dset, name, type) {
  event_dataset=paste(trajectory_group, name, sep="/")
  h5createDataset(file=outfile, dataset=dset, dims=c(length(dset)),
    H5type=type)
  h5write(dset, file=outfile, name=event_dataset)
}


write_events <- function(outfile, event_frame) {
  trajectory_name=next_dataset(outfile)
  trajectory_group=paste("trajectory", trajectory_name, sep="/")
  h5createGroup(outfile, trajectory_group)
  write_single_dataset(outfile, trajectory_group, event_frame$event,
    "Event", "H5T_STD_I64LE")
  write_single_dataset(outfile, trajectory_group, event_frame$who,
    "Who", "H5T_STD_I64LE")
  write_single_dataset(outfile, trajectory_group, event_frame$whom,
    "Whom", "H5T_STD_I64LE")
  write_single_dataset(outfile, trajectory_group, event_frame$when,
    "When", "H5T_STD_F64LE")
}


next_dataset <- function(infile) {
  datasets=h5ls(file=infile, recursive=TRUE)
  names=datasets[[2]]
  max_cnt=0
  for (name in names) {
    print(name)
    if (substr(name, start=1, stop=4)=="dset") {
      numeric_part=substr(name, start=5, stop=100)
      cnt=as.integer(numeric_part)
      print(paste("found", cnt))
      max_cnt<-max(max_cnt, cnt)
    }
  }
  paste("dset", max_cnt+1, sep="")
}

test_create_file <- function() {
  outfile="z.h5"
  create_file(outfile)
  events=data.frame(event=c(1:100), who=c(1:100), whom=c(1:100),
    when=as.double(c(1:100)))
  write_events(outfile, events)
  write_events(outfile, events)
  write_events(outfile, events)
}

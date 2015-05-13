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
  h5createDataset(file=outfile, dataset=event_dataset, dims=c(length(dset)),
    H5type=type)
  h5write(dset, file=outfile, name=event_dataset)
}

# locations are a spatstat dataset
write_locations <- function(outfile, locations) {
  cnt<-length(locations$x)
  print("About to create dataset")
  h5createDataset(file=outfile, dataset="locations",
    dims=c(cnt,2), H5type="H5T_IEEE_F64LE")
  m<-matrix(nr=cnt, nc=2)
  m[1:cnt,1]<-locations$x
  m[1:cnt,2]<-locations$y
  print("about to write to dataset")
  h5write(m, file=outfile, name="locations")
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
    "When", "H5T_IEEE_F64LE")
}


next_dataset <- function(infile) {
  datasets=h5ls(file=infile, recursive=TRUE)
  names=datasets[[2]]
  max_cnt=0
  for (name in names) {
    if (substr(name, start=1, stop=4)=="dset") {
      numeric_part=substr(name, start=5, stop=100)
      cnt=as.integer(numeric_part)
      max_cnt<-max(max_cnt, cnt)
    }
  }
  paste("dset", max_cnt+1, sep="")
}


# The functor arguments are locations and
# a data frame containing (event, who, whom, when)
foreach.trajectory <- function(infile, functor) {
  matrix_locations<-h5read(infile, "locations")
  lcnt<-dim(matrix_locations)[1]
  locations<-ppp(matrix_locations[1:lcnt,1], matrix_locations[1:lcnt, 2])
  datasets=h5ls(file=infile, recursive=TRUE)
  names=datasets[[2]]
  for (name in names) {
    if (substr(name, start=1, stop=4)=="dset") {
      trajectory_group<-paste("trajectory", name, sep="/")
      event<-h5read(infile, paste(trajectory_group, "Event", sep="/"),
        bit64conversion="int")
      who<-h5read(infile, paste(trajectory_group, "Who", sep="/"),
        bit64conversion="int")
      whom<-h5read(infile, paste(trajectory_group, "Whom", sep="/"),
        bit64conversion="int")
      when<-h5read(infile, paste(trajectory_group, "When", sep="/"))
      df<-data.frame(when=when, event=event, who=who, whom=whom)
      functor(locations, df)
    }
  }
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

#######################################################
##
## kt-prep-data.R
##
## prepare the kannisto-thatcher data for use
## in the old-age mortality analysis
##
#######################################################

root.dir <- file.path("..")

#root.dir <- file.path("~", "Dropbox", "oldageff", "paper-code")
code.dir <- file.path(root.dir, "code")
kt.dir <- file.path(root.dir, "data", "kt", "raw_scraped")
data.dir <- file.path(root.dir, "data")

library(mortfit)
library(ggplot2)
library(reshape2)
library(methods)

num.folds <- 5

set.seed(100)

## load the data we'll use
kt.dat <- load.kt.data(kt.desc.file=file.path(data.dir, "kt-data-touse.csv"),
                       kt.dir=kt.dir)

## now partition it into folds for fitting via cross-validation
kt.dat$num.folds <- num.folds

max.age <- 104
min.age <- 80
kt.dat$max.age <- max.age
kt.dat$min.age <- min.age

make.folds <- TRUE

## first go through each dataset and keep only given age range,
## and also only keep times with nonzero exposure
tmp <- plyr::llply(kt.dat$kt.data,
             function(x) {
               x@data <- subset(x@data, ((x@age.offset+age) <= max.age) &
                                        ((x@age.offset+age) >= min.age))
               x@data <- subset(x@data, Nx > 0)
               return(x)
             },
             .progress="text")

## then go through and partition each dataset into folds
tmp2 <- plyr::llply(tmp,
              function(x) {
                partition.into.folds(num.folds, x)
              },
              .progress="text")

kt.dat$kt.data <- tmp2
names(kt.dat$kt.data) <- plyr::laply(kt.dat$kt.data,
                               function(x) { x@name })

with(kt.dat,
     save(list=names(kt.dat),
          file=file.path(data.dir, "kt-data.RData")))





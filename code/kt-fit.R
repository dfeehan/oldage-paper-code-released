#!/usr/bin/env R

###########################################################
## kt-fit.R
##
## fit the mortality models to the kannisto-thatcher
## database
##
## uses the output of kt-prep-data.R
##
###########################################################

root.dir <- file.path("..")

data.dir <- file.path(root.dir, "data")
out.dir <- file.path(root.dir, "out")

library(mortfit)
library(methods)

################################################################
# Load the data
################################################################

parallel <- TRUE

## number of processors to claim if parallel is TRUE...
numproc <- 20

if(parallel) {
  library(doMC)
  registerDoMC(numproc)

}

load(file.path(data.dir, "kt-data.RData"))

## this is the name that the associated directory / files / etc 
## for this run will have
run.name <- "cohort"

## tweak these to produce fits using subsets of the
## models and the data
tofit.dat <- grab.kt.data(kt.data,
                          tags=c(type="cohort"))

out.dir <- path.expand(file.path(out.dir, run.name))

dir.create(out.dir, showWarnings = FALSE)

cat("out.dir is ", out.dir, "\n")

err.file <- file.path(out.dir, paste0("fit-", run.name, "-errors.log"))
cat("Errors for run starting ",
    format(Sys.time(), "%a %b %d %X %Y"),
    "\n",
    file=err.file, append=FALSE)

tofit.models <- binomial.models

this.seed <- 101009

logfile <- file(file.path(out.dir, paste0("fit-", run.name, "-data.log")), open="wt")
sink(logfile, split=TRUE)
sink(logfile, type="message")
cat("TOTAL DATASETS TO FIT: ", length(tofit.dat), "\n")
cat("TOTAL TO FIT: ", length(tofit.dat)*length(tofit.models), "\n")
cat("STARTING SEED: ", this.seed, "\n")
set.seed(this.seed)
system.time(
            fits <- plyr::llply(tofit.dat,
              function(thisdat) {
                cat("\n[DATA] starting data ", thisdat@name, "\n")
                res <- plyr::llply(tofit.models,
                             function(mod) {
                               #cat("\nstarting model", mod@name, " - ", thisdat@name, "\n")
                               thisres <- tryCatch(
                                          mort.fit(mod,
                                                   thisdat,
                                                   ##combinedFit,
                                                   optimMultipleFit,
                                                   M=10,
                                                   verbose=TRUE,
                                                   ##verbose=FALSE,
                                                   ignore.folded=FALSE),
                                          error=function(err) {
                                            cat("problem with: ", mod@name, " - ",
                                                thisdat@name, " - ",
                                                err$message,
                                                "\n========================\n",
                                                file=err.file,
                                                append=TRUE)
                                            return(NULL)
                                          })
                               return(thisres)
                             },
                             .parallel=FALSE)
                ## grab errors and keep track of them
                todrop <- plyr::laply(res, is.null)
                keptfits <- res[!todrop]
                return(create.mortalityFits(keptfits))
              },
              .progress="text",
              .parallel=parallel)
)

cat("\nEND OF RUN\n",
    file=err.file, append=TRUE)

save(tofit.models,
     tofit.dat,
     fits,
     #fits.df,
     file=file.path(out.dir, paste0("kt-", run.name, "-fit.RData")))

sink(NULL, type="message")
sink(NULL)


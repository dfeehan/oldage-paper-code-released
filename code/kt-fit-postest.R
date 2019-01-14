#!/usr/bin/env R

###########################################################
## kt-fit-postest.R
##
## summarize the results of fitting all of the models
## in preparation for actual analysis of results
##
## uses the output of fit-kt-data.R
##
###########################################################


root.dir <- file.path("..")

out.dir <- file.path(root.dir, "out")

library(methods)
library(mortfit)
library(tidyverse)
library(gridExtra)
library(ggplot2)

run.name <- "cohort"

out.dir <- path.expand(file.path(out.dir, run.name))

make.plots <- FALSE

################################################################
# Load the data
################################################################
load(file=file.path(out.dir, paste0("kt-", run.name, "-fit.RData")))

is.null.fit <- plyr::laply(fits, is.null)
fits <- fits[!is.null.fit]

fits.b <- plyr::llply(fits,
                function(x) {
                  fn <- names(x@fits)
                  x@fits <- x@fits[grep("Binomial", fn)]
                  return(x)
                })


fits.b.df <- plyr::llply(fits.b,
                   as.data.frame,
                   order.by="AIC")

######################################################################
## summarize the parameter estimates by functional form,
## to get some idea of plausible values...
##
#thetas <- plyr::llply(fits,
#                function(tf) {
#                  res <- plyr::ldply(tf@fits,
#                               function(x) {
#                                 thetas <- rep(NA,4)
#                                 thetas[1:length(x@theta.hat)] <- x@theta.hat
#                                 thismname <- gsub(".* - .+ - (.+) -.*",
#                                                     "\\1",
#                                                     x@name)
#                                 thislname <- gsub(".* - (.+) - .+ -.*",
#                                                     "\\1",
#                                                     x@name)                                 
#                                 thisdname <- gsub("(.*) - .+ - .+ -.*",
#                                                     "\\1",
#                                                     x@name)
#                                 cname <- gsub("(.*)-\\w-.*", "\\1", thisdname)
#                                 sex <- gsub("\\w+-(\\w)-.*","\\1",thisdname)
#                                 time <- gsub(".*-\\w-\\w+-(\\d+).*", "\\1", thisdname)
#                                 thisname <- x@name
#                                 
#                                 return(data.frame(modelname=thislname,
#                                                   hazardname=thismname,
#                                                   dataname=thisdname,
#                                                   country=cname,
#                                                   sex=sex,
#                                                   time=time,
#                                                   fullname=thisname,
#                                                   theta1=thetas[1],
#                                                   theta2=thetas[2],
#                                                   theta3=thetas[3],
#                                                   theta4=thetas[4]))
#                               })
#                })
#thetas <- do.call("rbind", thetas)
#thetas$.id <- NULL
#thetas.b <- plyr::dlply(subset(thetas, modelname=="Binomial"), plyr::.(hazardname), function(x) { x })
#
## grab the middle 90 pct of each param...
#tidx <- paste("theta",1:4,sep="")
#theta.b.distns <- plyr::llply(thetas.b,
#                      function(ff) {
#                        res <- plyr::aaply(ff[,tidx],
#                                     2,
#                                     quantile,
#                                     probs=c(.05,.5,.95),
#                                     na.rm=TRUE)
#                                     ##na.rm=FALSE)                        
#                      })

## save the list of parameters for each functional form to a file...
## this will be useful in thinking about plausible values later.

#save(thetas, theta.b.distns, thetas.b, 
#     file=file.path(out.dir, 
#		    paste0("kt-", run.name, "-parameters-byfunctionalform.RData")))

#theta.b.dat <- plyr::ldply(names(theta.b.distns),
#                   function(td) {
#                     thisdat <- theta.b.distns[[td]]
#                     return(data.frame(hazard=td,
#                                       param=rownames(thisdat),
#                                       thisdat))
#                   })
#write.csv(theta.b.dat, 
#          file.path(out.dir, 
#                    paste0("kt-", run.name, "-binomial-hazard-parameter-distns.csv")))

######################################################################
## summarize the results by country, time, and sex

## create vars for country, sex, year, likelihood (Poisson/Binomial/etc),
##             and hazard (Gompertz, Kannisto, etc)
fits.b.df <- plyr::llply(fits.b.df,
                 function(x) {
                   x$country <- gsub("(\\w+)\\-.*", "\\1", x$data)
                   x$year <- gsub("\\w+\\-.\\-\\w+\\-(\\d+)", "\\1", x$data)
                   x$sex <- gsub("\\w+\\-(.)\\-\\w+\\-\\d+", "\\1", x$data)
                   x$model.likelihood <- gsub("(\\w+) \\- .*", "\\1", x$model)
                   x$hazard <- gsub("\\w+ \\- (.*)", "\\1", x$model)
                   return(x)
                 })

## compute delta AIC and the AIC rank for each dataset we
## fit the models to
tmp.b <- plyr::llply(fits.b.df,
             function(x) {
               x$AIC.rank <- 1:nrow(x)
               x$delta.AIC <- x$AIC[1]
               x$delta.AIC <- x$AIC - x$AIC[1]
               x$BIC.rank <- rank(x$BIC)
               x$delta.BIC <- x$BIC - min(x$BIC)
               x$CV.rmse.Dx.rank <- rank(x$CV.rmse.Dx)
               x$delta.CV <- x$CV.rmse.Dx - min(x$CV.rmse.Dx)

               return(x)
             })

## make the results into one big, long dataset
uberfits.b <- do.call("rbind", tmp.b)

##########################################################
## relabel countries, hazards, and sexes using readable names
##

uberfits.b <- uberfits.b %>%
       mutate(country=
                 dplyr::recode(country,
                 'belgi' = 'Belgium',
                 'scotl' = 'Scotland',
                 'denma' = 'Denmark',
                 'franc' = 'France',
                 'germw' = 'W. Germany',
                 'italy' = 'Italy',
                 'japan' = 'Japan',
                 'nethe' = 'Netherlands',
                 'sweka' = 'Sweden',
                 'switz' = 'Switzerland'))

uberfits.b <- uberfits.b %>%
  mutate(sex=
           dplyr::recode(sex,
                         'f'='Females',
                         'm'='Males'))

summ.uberfits.b <- plyr::ddply(uberfits.b,
                       plyr::.(model,country,sex),
                       function(model.df) {
                         mda <- mean(model.df$delta.AIC)
                         mdb <- mean(model.df$delta.BIC)
                         mdcv <- mean(model.df$delta.CV)

                         mr <- mean(model.df$AIC.rank)
                         mcv <- mean(model.df$CV.rmse.Dx)

                         return(data.frame(model.likelihood=model.df$model.likelihood[1],
                                           hazard=model.df$hazard[1],
                                           mean.delta.AIC=mda,
                                           mean.delta.BIC=mdb,
                                           mean.delta.CV=mdcv,
                                           mean.AIC.rank=mr,
                                           mean.CV.RMSE.Dx=mcv,
                                           country=model.df$country[1]))
                       })


summ.uberfits.cs.b <- plyr::ddply(summ.uberfits.b,
                       plyr::.(country,sex),
                       function(cs.df) {
                         cs.df$mean.delta.AIC.rank <- rank(cs.df$mean.delta.AIC)
                         cs.df$mean.delta.BIC.rank <- rank(cs.df$mean.delta.BIC)
                         return(cs.df)
                       })

summ.uberfits.cs.b <- arrange(summ.uberfits.cs.b,
                         country, sex, mean.AIC.rank)

save(summ.uberfits.b,
     summ.uberfits.cs.b,
     tofit.dat,
     uberfits.b,
     fits,
     file=file.path(out.dir, 
                    paste0("kt-summary-", run.name, "-fits.RData")))

#########################
## optionally save plots of all the fits

if (make.plots) {
  
  pdf(file.path(out.dir, paste0(run.name, "-fit-plots-binomial.pdf")), height=12, width=12)
  plyr::l_ply(fits,
        function(x) {
          plot(x, Dx=FALSE)
        })
  dev.off()

}


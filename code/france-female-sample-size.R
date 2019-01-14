#!/usr/bin/env R

###########################################################
## france-females-sample-size.R
##
## illustrate changes in model fit as sample size changes
##
## uses the output of kt-prep-data.R
##
###########################################################

root.dir <- file.path("..")

data.dir <- file.path(root.dir, "data")
out.dir <- file.path(root.dir, "out")

library(mortfit)
#library(plyr)
library(methods)
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)

parallel <- TRUE
#parallel <- FALSE

## number of processors to claim if parallel is TRUE...
numproc <- 20

if(parallel) {
  library(doMC)
  registerDoMC(numproc)
  
}

tofit.models <- binomial.models

################################################################
# Load the data
################################################################

load(file.path(data.dir, "kt-data.RData"))

## this is the name that the associated directory / files / etc 
## for this run will have
run.name <- "cohort"

## tweak these to produce fits using subsets of the
## models and the data
french_female_cohorts <- grab.kt.data(kt.data,
                          tags=c(type="cohort",
                                 country='franc',
                                 sex='f',
                                 time='1872'))

# given a mortalityData or mortalityDataFolded object,
# simulate reducing the sample size by a factor of frac
reduce.cohort <- function(base_cohort, frac) {

  sampled_cohort <- base_cohort
  sampled_cohort_data <- sampled_cohort@data
  sampled_cohort_data$Nx <- round(frac * sampled_cohort_data$Nx)
  sampled_cohort_data$Dx <- round(frac * sampled_cohort_data$Dx)
  sampled_cohort@data <- sampled_cohort_data
  sampled_cohort@name <- paste0(sampled_cohort@name, " - frac ", frac)

  return(sampled_cohort)

}

pared.cohorts <- map(french_female_cohorts,
    function(base_cohort) {
      
      these.pared.cohorts <- map(seq(.05, 1, by=.05),
                                 ~ reduce.cohort(base_cohort, .x))
      
      return(these.pared.cohorts)
    })

pared.cohorts <- flatten(pared.cohorts)

fits <- plyr::llply(pared.cohorts,
                    function(thisdat) {
                        
                      res <- plyr::llply(tofit.models,
                                         function(mod) {
                                           thisres <- tryCatch(
                                             mort.fit(mod,
                                                      thisdat,
                                                      ##combinedFit,
                                                      optimMultipleFit,
                                                      M=10,
                                                      verbose=TRUE,
                                                      ignore.folded=TRUE),
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
                                         })
                      return(create.mortalityFits(res))
                    },
                    .parallel=parallel)

fits.df <- plyr::llply(fits,
                       as.data.frame,
                       order.by="AIC")

fits.df <- plyr::llply(fits.df,
                         function(x) {
                           x$country <- gsub("(\\w+)\\-.*", "\\1", x$data)
                           x$year <- gsub("\\w+\\-.\\-\\w+\\-(\\d+) - .*", "\\1", x$data)
                           x$sex <- gsub("\\w+\\-(.)\\-\\w+\\-\\d+", "\\1", x$data)
                           x$model.likelihood <- gsub("(\\w+) \\- .*", "\\1", x$model)
                           x$hazard <- gsub("\\w+ \\- (.*)", "\\1", x$model)
                           x$frac <- gsub("(.*)(frac )(.*)", "\\3", x$data)
                           return(x)
                         })

## compute delta AIC and the AIC rank for each dataset we
## fit the models to
fits.df <- plyr::llply(fits.df,
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
fits.df <- do.call("rbind", fits.df)

fits.df <- fits.df %>% filter(! hazard %in% c("Weibull", "Constant Hazard"))

## france females: 1866 to 1892
cur.year <- 1872

pdf(file=file.path(out.dir, run.name, "for-paper", "sample-size-plot.pdf"))
  toplot <- fits.df %>% 
    filter(year==cur.year) %>%
    mutate(frac = as.numeric(frac))
  fracplot <- ggplot(toplot, aes(x=frac, y=delta.AIC, color=hazard, group=hazard)) +
    #geom_line(aes(x=frac, y=delta.AIC, color=hazard, group=hazard)) +
    coord_cartesian(xlim = c(min(toplot$frac), max(toplot$frac) + .02)) +
    geom_line() +
    geom_text_repel(
      data = subset(toplot, frac == max(frac)),
      aes(label = hazard),
      size = 3.5,
      direction="y",
      nudge_x = 1,
      nudge_y = 1,
      segment.color = NA,
      show.legend=FALSE
    ) +
    xlab("Cohort size (% of original cohort)") +
    ylab(expression(Delta*AIC)) +
    ggtitle(paste0("French Females ", cur.year)) +
    guides(color=guide_legend(title='')) +
    scale_x_continuous(labels = scales::percent,
                       breaks=seq(from=0,to=1,by=.1)) +
    theme_minimal() +
    theme(legend.position=c(.1,.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(fracplot)
dev.off()

save(tofit.models,
     pared.cohorts,
     fits,
     fits.df,
     file=file.path(out.dir, run.name, paste0("kt-samplesize-fit.RData")))

#load(file=file.path(out.dir, paste0("kt-", run.name, "-fit.RData")))





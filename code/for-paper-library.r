
##########################################################
## get fits, plots, etc for a single age/sex/country example
##
do.example <- function(ex.dat, ex.name, all.ff) {

    ex.plot <- plot(ex.dat)
    fits.ex <- csy.fit(ex.dat=ex.dat, ex.name=ex.name)
    ex <- csy.example(fits.ex, all.ff=all.ff)
    tp.ex <- csy.fits.2panelplot(fits.ex, all.ff, line.size=1.5, line.alpha=.8) 
    errp.ex <- csy.fits.errorplots(fits=fits.ex, all.ff=all.ff)

    return(list(data=ex.dat,
                name=ex.name,
                fn=ex.dat@name,
                fits=fits.ex,
                onebyone=ex,
                together=tp.ex,
                errs=errp.ex))
}

##########################################################
## fit models to an example country-year
##
csy.fit <- function(ex.dat, ex.name) {

    tofit.models <- binomial.models

    system.time(
                all.fits <- plyr::llply(list(ex.dat),
                  function(thisdat) {
                    cat("\nstarting data ", thisdat@name, "\n")
                    res <- plyr::llply(tofit.models,
                                 function(mod) {
                                   cat("\nstarting model", mod@name, " - ", thisdat@name, "\n")
                                   thisres <- tryCatch(
                                              mort.fit(mod,
                                                       thisdat,
                                                       ##combinedFit,
                                                       ##M=3,
                                                       cpp=FALSE,
                                                       optimFit,
                                                       ##gridSearchFit,
                                                       keep.fold.fits=TRUE,
                                                       ##verbose=TRUE,
                                                       verbose=FALSE,
                                                       ## for now, don't fit CV versions
                                                       ## (for speed, until we get kinks ironed out)
                                                       ignore.folded=FALSE,
                                                       #ignore.folded=TRUE,
                                                       random.start=FALSE),
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
                  .parallel=FALSE)
    )
    names(all.fits) <- plyr::laply(list(ex.dat), function(x) x@name)

    fits <- all.fits[[1]]

    return(fits)

}

# hack to avoid blue circles in publication
tmpplot <- function (x, y = NULL, Dx = FALSE) 
{
  if (!Dx) {
    #dataplot <- plot(x@data, x@model@hazard, theta = x@theta.hat)
    theta <- x@theta.hat
    y <- x@model@hazard
    x <- x@data
    
    if (is.null(theta)) {
      theta <- y@theta.default
    }
    age <- x@data$age
    age.offset <- x@age.offset
    mu <- y@haz.fn(theta, age)
    pi <- mortfit:::haz.to.prob(y@haz.fn, theta, age)
    haz.dat <- data.frame(age = age + age.offset, mu = mu, pi = pi)
    plot.dat <- x@data
    plot.dat$age <- plot.dat$age + age.offset + 0.5
    dataplot <- ggplot(plot.dat) + 
      geom_line(aes(x = age, y = mu), 
                                          data = haz.dat) + 
      geom_point(aes(x = age, y = (Dx/(Nx - 0.5 * Dx)), size = Nx), 
                 pch = 1) + labs(x = "age", y = "hazard /\ncentral death rate") + 
      ggtitle(paste0(x@name,"\n", y@name, "\n", "(", paste(round(theta, 4), collapse = ", "),")")) + 
      scale_size_area() 
    
  }
  else {
    obsDx <- x@data@data$Dx
    obsNx <- x@data@data$Nx
    fitDx <- x@fitted.values@fitted.Dx
    ages <- x@fitted.values@age + x@fitted.values@age.offset
    dataplot <- ggplot(data = data.frame(obsDx = obsDx, fitDx = fitDx, 
                                         ages = ages, obsNx = obsNx)) + 
      geom_point(aes(x = ages, 
                     y = obsDx, size = obsNx), 
                 #color = "blue", 
                 pch = 1) + 
      geom_line(aes(x = ages, y = fitDx), color = "red") + 
      geom_point(aes(x = ages, y = fitDx), color = "red", 
                 pch = 3) + xlab("age") + ylab("number of deaths") + 
      scale_size_area() + ggtitle(x@name)
  }
  return(dataplot)
}
#########################################################
## example of model fits to one country-sex-cohort
##
## fits - should be the fits to the given country-sex-cohort
csy.example <- function(ex.fit, all.ff) {

  ## make a set of plots for this country-sex-year
  ex.fit.names <- plyr::laply(ex.fit@fits, function(x) { x@model@name })
  tograb <- which(grepl("Binomial", ex.fit.names) & 
                  grepl(paste(all.ff,collapse="|"), ex.fit.names))
  ex.fit@fits <- ex.fit@fits[tograb]

                      
  resplots <- plyr::llply(ex.fit@fits,
                    function(x) {
                      thisn <- gsub(".* - (.*)$", "\\1", x@model@name)
                      # hack to fix quadratic -> log-quadratic 
                      if (thisn == 'Quadratic') {
                        #thisn <- "Log-Quadratic"
                        thisn <- "Log-quadratic"
                      }
                      if (thisn == 'Log-Quadratic') {
                        #thisn <- "Log-Quadratic"
                        thisn <- "Log-quadratic"
                      }
                      
                      #thisp <- plot.mortalityDataWithFit(x,Dx=FALSE)
                      
                      thisp <- tmpplot(x,Dx=FALSE)
                      
                      thisp <- thisp + 
                        theme(legend.position="none",
                              axis.title.y = element_text(face='bold')) + 
                        ggtitle(thisn) +
                        ylab("Hazard/\nCentral Death Rate")
                      return(thisp)
                    })

  resplots.dx <- plyr::llply(ex.fit@fits,
                    function(x) {
                      thisn <- gsub(".* - (.*)$", "\\1", x@model@name)
                      # hack to fix quadratic -> log-quadratic 
                      if (thisn == 'Quadratic') {
                        #thisn <- "Log-Quadratic"
                        thisn <- "Log-quadratic"
                      }
                      #thisp <- plot.mortalityDataWithFit(x,Dx=TRUE)
                      thisp <- tmpplot(x,Dx=TRUE)
                      thisp <- thisp + theme(legend.position="none") + ggtitle(thisn)
                      return(thisp)
                    })

  fits.df <- as.data.frame(ex.fit, order.by="AIC")

  fits.df <- transform(fits.df,
                       country = gsub("(\\w+)\\-.*", "\\1", data),
                       year = gsub("\\w+\\-.\\-\\w+\\-(\\d+)", "\\1", data),
                       sex = gsub("\\w+\\-(.)\\-\\w+\\-\\d+", "\\1", data),
                       model.likelihood = gsub("(\\w+) \\- .*", "\\1", model),
                       hazard = gsub("\\w+ \\- (.*)", "\\1", model))

  fits.df <- transform(fits.df,
                       country = gsub("(\\w+)\\-.*", "\\1", data),
                       year = gsub("\\w+\\-.\\-\\w+\\-(\\d+)", "\\1", data),
                       sex = gsub("\\w+\\-(.)\\-\\w+\\-\\d+", "\\1", data),
                       model.likelihood = gsub("(\\w+) \\- .*", "\\1", model),
                       hazard = gsub("\\w+ \\- (.*)", "\\1", model))

  fits.df$AIC.rank <- 1:nrow(fits.df)

  uberfits <- transform(fits.df,
                        SSE.Dx.rank = rank(SSE.Dx),
                        delta.AIC = AIC - AIC[1],
                        BIC.rank = rank(BIC),
                        delta.BIC = BIC - min(BIC),
                        CV.rmse.Dx.rank = rank(CV.rmse.Dx))

  ## make a table of fit summaries for this country-sex-year
  resdat <- subset(uberfits,
                   hazard %in% all.ff)

  #resdat$RMSE <- (1/(104-80+1))*sqrt(resdat$SSE.Dx)
  num.ages <- nrow(ex.fit@data@data)
  resdat$RMSE <- (1/num.ages)*sqrt(resdat$SSE.Dx)
  
  rownames(resdat) <- resdat$hazard

  resdat <- subset(resdat,
                   select=c(log.likelihood,
                            SSE.Dx,                   
                            SSE.Dx.rank,
                            AIC,
                            AIC.rank,
                            delta.AIC,
                            BIC,
                            delta.BIC,
                            BIC.rank,
                            CV.rmse.Dx,
                            CV.rmse.Dx.rank))

  colnames(resdat) <- c("log-likelihood",
                        "SSE",
                        "SSE rank",
                        "AIC",
                        "AIC rank",
                        "$\\Delta$ AIC",
                        "BIC",
                        "$\\Delta$ BIC",                        
                        "BIC rank",
                        "CV RMSE(Dx)",
                        "CV rank")

  return(list(plots=resplots,
              plots.dx=resplots.dx,
              table=resdat))
}

#### example 2 panel plot (hazards and deaths)
csy.fits.2panelplot <- function(x, all.ff, line.size=1.5, line.alpha=.8, point.size=5) {

    obsdat <- x@data@data
    obsdat$age <- obsdat$age + x@data@age.offset
    obsdat <- obsdat[, c('Nx', 'Dx', 'age')]
    ## use central death rate as observed hazard
    obsdat <- transform(obsdat,
                        mu = Dx/(Nx - 0.5*Dx))
    obsdat$model <- 'observed'

    fitvals <- plyr::llply(x@fits,
                     function(thisfit) {
                         modelName <- thisfit@model@name
                         fitted.Dx <- thisfit@fitted.values@fitted.Dx
                         raw.ages <- thisfit@fitted.values@age
                         fitted.ages <- raw.ages + thisfit@fitted.values@age.offset

                         mu <- thisfit@model@hazard@haz.fn(thisfit@theta.hat,
                                                           raw.ages)

                         return(data.frame(model=modelName,
                                           Dx=fitted.Dx,
                                           # repeat Nx from observed data
                                           Nx=obsdat$Nx,
                                           mu=mu,
                                           age=fitted.ages))
                     })

    plotdat <- do.call("rbind", fitvals)
    rownames(plotdat) <-  NULL

    ## remove prefixes to model description
    ## (TODO -- might eventually want to make this more general)
    plotdat$model <- gsub(".* - (.*)$", "\\1", plotdat$model)

    plotdat <- subset(plotdat, model %in% all.ff)

    plotdat <- rbind(obsdat, plotdat)
    ## plot halfway through interval
    plotdat$age <- plotdat$age + 0.5

    ## TODO -- left off here: trying to figure out why gompertz is strange in this
    ##         combined plot

    # plot for Dx=FALSE
    tp.haz <- ggplot(plotdat) +
          geom_line(data=subset(plotdat, model != 'observed'),
                    aes(x=age, y=mu, color=model), size=line.size, alpha=line.alpha) +
          geom_point(data=subset(plotdat, model == 'observed'),
                    aes(x=age, y=mu, size=Nx), pch=1) +
          #geom_point(data=subset(plotdat, model != 'observed'),
          #          aes(x=age, y=mu, color=model, pch=model)) +
          scale_size_area() + 
          scale_color_hazards_slides + 
          theme(legend.position=c(.15, .65), panel.border=element_blank()) +
          guides(color=guide_legend(order=1), 
                 size="none") +
                 #size=guide_legend(title="", order=2)) +
          ylab("hazard")

    # plot for Dx=TRUE
    tp.dx <- ggplot(plotdat) +
          #geom_line(data=subset(plotdat, model != 'observed'),
          #          aes(x=age, y=Dx, color=model), size=line.size, alpha=line.alpha) +
          geom_point(data=subset(plotdat, model != 'observed'),
                    aes(x=age, y=Dx, color=model), pch='x', size=point.size) +
          geom_point(data=subset(plotdat, model=='observed'),
                     aes(x=age, y=Dx, size=Nx), pch=1) +
          scale_size_area() + 
          scale_color_hazards_slides +
          theme(legend.position=c(.8, .65), panel.border=element_blank()) +
          guides(color=guide_legend(order=1), 
                 size="none") +
                 #size=guide_legend(title="", order=2)) +
          ylab("number of deaths")

    return(list(plots=list(plot.haz=tp.haz, plot.dx=tp.dx),
                data=plotdat))
}


#### make a list of error plots, one for each hazard
csy.fits.errorplots <- function(fits, all.ff) {
    names(all.ff) <- all.ff
    res <- plyr::llply(all.ff,
                 function(this.ff) {
                     fit.idx <- which(str_detect(names(fits@fits), this.ff))
                     csy.fits.errorplot(fits[[fit.idx]], this.ff=this.ff)                   
                 })
}

#### example error illustration plot for an invidual hazard
csy.fits.errorplot <- function(this.fit, this.ff, title=NULL) {

    obsdat <- this.fit@data@data
    obsdat$age <- obsdat$age + this.fit@data@age.offset
    obsdat <- obsdat[, c('Nx', 'Dx', 'age')]
    ## use central death rate as observed hazard
    obsdat <- transform(obsdat,
                        cdr = Dx/(Nx - 0.5*Dx))

    modelName <- this.fit@model@name
    fitted.Dx <- this.fit@fitted.values@fitted.Dx
    raw.ages <- this.fit@fitted.values@age
    fitted.ages <- raw.ages + this.fit@fitted.values@age.offset

    mu <- this.fit@model@hazard@haz.fn(this.fit@theta.hat,
                                    raw.ages)

    plotdat <- data.frame(model=modelName, 
                          Dx.fit=fitted.Dx,
                          mu.fit=mu, 
                          age=fitted.ages)

    ## remove prefixes to model description
    ## (TODO -- might eventually want to make this more general)
    plotdat$model <- gsub(".* - (.*)$", "\\1", plotdat$model)

    plotdat <- subset(plotdat, model %in% this.ff)

    plotdat <- merge(obsdat, plotdat, by=c("age"))
    plotdat <- transform(plotdat,
                         Dx.fit.err=Dx.fit - Dx)

    ## plot halfway through interval
    plotdat$age <- plotdat$age + 0.5

    ## TODO -- left off here: trying to figure out why gompertz is strange in this
    ##         combined plotmethods

    pred.pch <- 4

    # plot for Dx=TRUE
    tp.dx.noerr <- ggplot(plotdat) +
                 geom_point(aes(x=age, y=Dx.fit, color=model), pch=pred.pch, size=2) +
                 ##geom_line(aes(x=age, y=Dx.fit, color=model), size=2) +
                 geom_point(aes(x=age, y=Dx, size=Nx), pch=1) +
                 scale_size_area() + 
                 scale_color_hazards_slides +
                 theme(legend.position=c(.8, .8)) +
                 guides(color='none') +
                 ylab("number of deaths") + ggtitle(title)

    tp.dx.err <- ggplot(plotdat) +
                 geom_point(aes(x=age, y=Dx.fit, color=model), pch=pred.pch, size=2) +
                 ##geom_line(aes(x=age, y=Dx.fit, color=model), size=2) +
                 geom_point(aes(x=age, y=Dx, size=Nx), pch=1) +
                 geom_linerange(aes(x=age, ymin=Dx, ymax=Dx + Dx.fit.err),
                                color='red', size=1.5) +
                 scale_size_area() + 
                 scale_color_hazards_slides +
                 theme(legend.position=c(.8, .8)) +
                 guides(color='none') +
                 ylab("number of deaths") + ggtitle(title)

    return(list(name=this.ff,
                plots=list(without.err=tp.dx.noerr, 
                           with.err=tp.dx.err), 
                data=plotdat))
}


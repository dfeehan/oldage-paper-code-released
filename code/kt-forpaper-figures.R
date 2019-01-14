###########################################################
## kt-forpaper-figures.R
##
## produce the figures for the old-age mortality paper
##
## uses the output of fit-kt-data.R and kt-fit-postest.R
##
###########################################################

root.dir <- ".."

library(tidyverse)
library(mortfit)
library(devtools)
library(gridExtra)
library(xtable)
library(GGally)
library(cluster)
library(forcats)
library(pander)
library(ggrepel)
library(methods)
#library(ggedit)

run.name <- "cohort"

out.dir <- path.expand(file.path(root.dir, "out", run.name))

out.summ.dir <- out.dir
paper.out.dir <- file.path(out.dir, "for-paper")

code.dir <- path.expand(file.path(root.dir, "code"))

dir.create(paper.out.dir, showWarnings = FALSE)

## Load the data
load(file=file.path(out.summ.dir, 
                    paste0("kt-summary-", run.name, "-fits.RData")))
#load(file=file.path(out.summ.dir, 
#                    paste0("kt-", run.name, "-parameters-byfunctionalform.RData")))

## load the ggplot theme...
source(file.path(root.dir, "code", "oldage_theme.R"))
theme_set(oldage_theme)

## we're only considering this subset of the functional forms
## in our analysis
all.ff <- c("Gompertz", "Kannisto", "Makeham",
            "Weibull",
            "Log-Quadratic", "Perks", "Beard",
            "Logistic", "Lynch-Brown")


## has the csy.example function
source(file.path(code.dir, "for-paper-library.r"))

############
## cohort sizes at ages 80 and 95 for descriptive table
cohort.sizes <- plyr::ldply(tofit.dat, 
                            function(x) { 
                              res <- data.frame(country = x@tags$country,
                                                time = as.numeric(x@tags$time),
                                                sex = x@tags$sex,
                                                size80 = x@data$Nx[1],
                                                size95 = x@data$Nx[15])
                              return(res)
                            })

cohort.sizes %>%
  group_by(country) %>%
  summarize(
            cohort.start = min(time),
            cohort.end = max(time),
            number = 2*(cohort.end-cohort.start+1),
            n.80.mean = format(round(mean(size80)),big.mark=','),
            n.95.mean = format(round(mean(size95)),big.mark=',')
    ) %>%
  arrange(-number) -> avg.cohort.sizes

tosave.tab <- avg.cohort.sizes %>%
  ungroup() %>%
  select("Country"=country, 
         "Cohort start"=cohort.start,
         "Cohort end"=cohort.end,
         "Number of cohorts"=number,
         "Avg. size (80)"=n.80.mean, 
         "Avg. size (95)"=n.95.mean) %>%
  mutate(Country = fct_recode(Country, 
                              "Japan"="japan",
                              "West Germany"="germw",
                              "Italy"="italy",
                              "France"="franc",
                              "Netherlands"="nethe",
                              "Belgium"="belgi",
                              "Sweden"="sweka",
                              "Switzerland"="switz",
                              "Denmark"="denma",
                              "Scotland"="scotl"))


cat(
pandoc.table.return(tosave.tab %>% as_data_frame(), 
             emphasize.rownames=FALSE,
             keep.line.breaks=TRUE,
             split.tables=Inf,
             caption="
Cohorts from the Kannisto-Thatcher Database on Old Age Mortality
used in this analysis. Only cohorts with the highest-quality data and
deaths reported at least up to age 105 by single year of age are included.
The final dataset has data for 360 country-sex-cohorts from 10 countries. 
@jdanov_kannisto_2008 has a detailed discussion of data quality and
Appendix [-@sec:ap-data] has more information about how the analysis dataset
was constructed. {#tbl:kt-data}
             "),
file=file.path(paper.out.dir, "analysis-sample.md")
)

panderOptions('big.mark',"")

#########
## save a graph and table of the specific example used in the paper
## (Danish males, 1895)
ex1.data <- "denma-m-cohort-1895"

ex1 <- csy.example(fits[[ex1.data]], all.ff=all.ff)

## manually set y range on plots for consistency when printing on one page
ex1$plots <- plyr::llply(ex1$plots,
                   function(x) { 
                     return(x + ylim(0,0.3))
                   })

pdf(file.path(paper.out.dir, paste0(ex1.data, ".pdf")),
    height=10, width=8)
do.call(grid.arrange, c(ex1$plots, list(ncol=2)))
dev.off()


### save table of specific example for main text 
tosave.tab <- ex1$table[,1:6]

cat(
pandoc.table.return(tosave.tab %>% as_data_frame(), 
             digits=c(6, 6, 1, 7, 1, 2),
             emphasize.rownames=FALSE,
             split.tables=Inf,
             caption="
             Measurements of model fit for the cohort of Danish males born in 1895. {#tbl:denma95} 
             "),
file=file.path(paper.out.dir, paste0(ex1.data, ".md")))

##############
### AIC -- good vs bad plots
delta.goodbad.summ <- uberfits.b %>%
  group_by(sex, hazard) %>%
  summarize(
            mean.delta.AIC.lt2 = mean(delta.AIC <= 2.0),
            mean.delta.AIC.gt10 = mean(delta.AIC >= 10.0)
            ) %>%
  arrange(sex, desc(mean.delta.AIC.lt2))

write_csv(delta.goodbad.summ,
          path=file.path(paper.out.dir, "summ-delta-goodbad.csv"))


toplot <- delta.goodbad.summ %>%
  gather(range, fraction, -hazard, -sex) 
 delta.aic.goodbad.plot <-
  ggplot(delta.goodbad.summ %>%
         filter(hazard != 'Constant Hazard')) +
  geom_point(aes(x=mean.delta.AIC.lt2, y=mean.delta.AIC.gt10, shape=sex)) +
  geom_text_repel(aes(x=mean.delta.AIC.lt2, y=mean.delta.AIC.gt10, label=hazard)) +
  xlim(0,1) +
  ylim(0,1) +
  xlab(expression(bold(Fraction~of~Cohorts~With~Delta*AIC<=2))) +
  ylab(expression(bold(Fraction~of~Cohorts~With~Delta*AIC>10))) +
  scale_shape_manual(name="",
                     values=c('Females'=1, 'Males'=2)) +
  theme_minimal() +
  theme(legend.position="bottom",
        axis.title.y = element_text(face='bold')) +
  guides(shape = guide_legend(override.aes = list(size=3)))

pdf(file.path(paper.out.dir, "delta-aic-goodbad.pdf"),
    width=8, height=8)
print(delta.aic.goodbad.plot)
dev.off()

##############
### AIC -- good vs bad plots by country

delta.goodbad.country.summ <- uberfits.b %>%
  group_by(country, sex, hazard) %>%
  summarize(
            mean.delta.AIC.lt2 = mean(delta.AIC <= 2.0),
            mean.delta.AIC.gt10 = mean(delta.AIC >= 10.0)
            ) %>%
  arrange(sex, desc(mean.delta.AIC.lt2))


set.seed(101)
delta.aic.goodbad.country.plot <-
  ggplot(delta.goodbad.country.summ %>%
         filter(! hazard %in% c('Constant Hazard', 'Weibull')) %>%
         #filter(! hazard %in% c('Constant Hazard')) %>%
         filter(country %in% c('Denmark', 'France', 'Italy', 'Netherlands', 'Sweden')) %>%
         mutate(hazard = dplyr::recode(hazard,
                                        'Log-Quadratic'='Log-quadratic'))
         ) +
  geom_point(aes(x=mean.delta.AIC.lt2, y=mean.delta.AIC.gt10)) +
  geom_text_repel(aes(x=mean.delta.AIC.lt2, 
                      y=mean.delta.AIC.gt10, 
                      label=hazard),
                  max.iter=1000,
                  size=2.5) +
  facet_grid(country ~ sex) +
  #facet_grid(sex ~ country) +
  xlim(0,1) +
  ylim(0,1) +
  xlab(expression(bold(Fraction~of~Cohorts~With~Delta*AIC<=2))) +
  ylab(expression(bold(Fraction~of~Cohorts~With~Delta*AIC>10))) +
  theme_bw(base_size=10) +
  coord_equal() +
  theme(legend.position="bottom",
        strip.background=element_rect(fill=NA),
        strip.text=element_text(face='bold')) 
delta.aic.goodbad.country.plot

ggsave(plot=delta.aic.goodbad.country.plot,
       filename=file.path(paper.out.dir, "delta-aic-goodbad-country.pdf"),
       height=10, width=6.5)

##############
### BIC -- good vs bad plots

delta.bic.goodbad.country.summ <- uberfits.b %>%
  group_by(country, sex, hazard) %>%
  summarize(
            mean.delta.BIC.lt2 = mean(delta.BIC <= 2.0),
            mean.delta.BIC.gt10 = mean(delta.BIC >= 10.0)
            ) %>%
  arrange(sex, desc(mean.delta.BIC.lt2))

set.seed(101)
delta.bic.goodbad.country.plot <-
  ggplot(delta.bic.goodbad.country.summ %>%
         filter(! hazard %in% c('Constant Hazard', 'Weibull')) %>%
         filter(country %in% c('Denmark', 'France', 'Italy', 'Netherlands', 'Sweden'))) +
  geom_point(aes(x=mean.delta.BIC.lt2, y=mean.delta.BIC.gt10)) +
  geom_text_repel(aes(x=mean.delta.BIC.lt2, 
                      y=mean.delta.BIC.gt10, 
                      label=hazard),
                  max.iter=1000,
                  size=2.5) +
  facet_grid(country ~ sex) +
  xlim(0,1) +
  ylim(0,1) +
  xlab(expression(Fraction~of~cohorts~with~Delta*BIC<=2)) +
  ylab(expression(Fraction~of~cohorts~with~Delta*BIC>10)) +
  theme_bw(base_size=10) +
  coord_equal() +
  theme(legend.position="bottom",
        strip.background=element_rect(fill=NA)) 

ggsave(plot=delta.bic.goodbad.country.plot,
       filename=file.path(paper.out.dir, "delta-bic-goodbad-country.pdf"),
       height=10, width=6.5)

### BIC -- good vs bad plots
delta.bic.goodbad.summ <- uberfits.b %>%
  group_by(sex, hazard) %>%
  summarize(
            mean.delta.BIC.lt2 = mean(delta.BIC <= 2.0),
            mean.delta.BIC.gt10 = mean(delta.BIC >= 10.0)
            ) %>%
  arrange(sex, desc(mean.delta.BIC.lt2))

delta.bic.goodbad.plot <-
  ggplot(delta.bic.goodbad.summ %>%
         filter(hazard != 'Constant Hazard')) +
  geom_point(aes(x=mean.delta.BIC.lt2, y=mean.delta.BIC.gt10, shape=sex)) +
  geom_text_repel(aes(x=mean.delta.BIC.lt2, y=mean.delta.BIC.gt10, label=hazard)) +
  xlim(0,1) +
  ylim(0,1) +
  xlab(expression(Fraction~of~cohorts~with~Delta*BIC<=2)) +
  ylab(expression(Fraction~of~cohorts~with~Delta*BIC>10)) +
  scale_shape_manual(name="",
                     values=c('Females'=1, 'Males'=2)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  guides(shape = guide_legend(override.aes = list(size=3)))

pdf(file.path(paper.out.dir, "delta-bic-goodbad.pdf"),
    width=8, height=8)
print(delta.bic.goodbad.plot)
dev.off()


##########################################################
## figure with boxplots of the delta AIC distributions
## over the entire dataset

## we want a log scale; recode 0s to a very small value so that we can still plot them
smallest.nonzero <- function(x) { min(x[x != 0]) }
eps <- smallest.nonzero(uberfits.b$delta.AIC)

# this function is taken from this thread:
# https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4
#
# it puts axis labels in scientific notation
fancy_scientific <- function(l) { 
  # turn in to character string in scientific notation 
  l <- format(l, scientific = TRUE) 
  # quote the part before the exponent to keep all the digits 
  l <- gsub("^(.*)e", "'\\1'e", l) 
  
  # remove the + in the exponent, if there is one
  l <- gsub("\\+", "", l)
  
  # turn the 'e+' into plotmath format 
  l <- gsub("e", "%*%10^", l) 
  # return this as an expression 
  parse(text=l) 
}

########### 
## which models have the lowest median Delta AIC across all
## cohorts and both sexes?
## NB: WE REFER TO THIS IN THE TEXT

tosave.med.deltaaic <- 
  uberfits.b %>% 
  group_by(hazard) %>%
  summarize(median.Delta.AIC = median(delta.AIC),
            iqr.Delta.AIC = IQR(delta.AIC)) %>%
  filter(! hazard %in% c('Constant Hazard')) %>%
  arrange(median.Delta.AIC) %>%
  ungroup() %>%
  select("Model"=hazard, 
         "Median $\\Delta$ AIC"=median.Delta.AIC,
         "IQR $\\Delta$ AIC"=iqr.Delta.AIC)

cat(
#pandoc.table.return(data.frame(tosave.tab), 
pandoc.table.return(tosave.med.deltaaic %>% as_data_frame(), 
             emphasize.rownames=FALSE,
             keep.line.breaks=TRUE,
             split.tables=Inf,
             caption="
Median $\\Delta$ AIC by model, across all cohorts. {#tbl:median-delta-aic}
             "),
file=file.path(paper.out.dir, "median-delta-aic.md")
)

tosave.med.deltaaic.bysex <- 
  uberfits.b %>% 
  group_by(sex, hazard) %>%
  summarize(median.Delta.AIC = median(delta.AIC),
            iqr.Delta.AIC = IQR(delta.AIC)) %>%
  filter(! hazard %in% c('Constant Hazard')) %>%
  arrange(sex, median.Delta.AIC) %>%
  ungroup() %>%
  select("Sex" = sex,
         "Model"=hazard, 
         "Median $\\Delta$ AIC"=median.Delta.AIC,
         "IQR $\\Delta$ AIC"=iqr.Delta.AIC)

cat(
#pandoc.table.return(data.frame(tosave.tab), 
pandoc.table.return(tosave.med.deltaaic.bysex %>% as_data_frame(), 
             emphasize.rownames=FALSE,
             keep.line.breaks=TRUE,
             split.tables=Inf,
             caption="
Median $\\Delta$ AIC by model and sex. {#tbl:median-delta-aic}
             "),
file=file.path(paper.out.dir, "median-delta-aic-bysex.md")
)

########################
## boxplots showing the distn of delta AIC
## for all models except constant hazard and Weibull
## (NB: y axis is logged)

delta.aic.summ.box.plot <- uberfits.b %>% 
  filter(! hazard %in% c('Constant Hazard')) %>%
  #mutate(hazard = forcats::fct_reorder(hazard, sex, delta.AIC)) %>%
  ggplot(.,
         ## NB: this orders by the median delta.AIC across both sexes
         ## (which is why neither panel is perfectly in order)
         aes(x=forcats::fct_reorder(hazard, delta.AIC), 
             y=delta.AIC + eps, 
             group=interaction(hazard, sex))) +
  geom_boxplot() +
  scale_y_log10(labels=fancy_scientific) +
  #coord_trans(y="log10") +
  facet_grid(~ sex) +
  ylab(expression(bold(log[10](Delta*AIC)))) +
  xlab("") +
  theme_minimal(base_size=14) +
  theme(axis.text.x=element_text(angle=90, hjust=1, face='bold'),
        strip.text=element_text(face='bold'))
delta.aic.summ.box.plot

pdf(file.path(paper.out.dir, "delta-aic-boxplot.pdf"),
    width=9, height=4.5)
print(delta.aic.summ.box.plot)
dev.off()


#######################
## boxplots showing the distn of delta BIC
## for all models except constant hazard and Weibull
## (NB: y axis is logged)
delta.bic.summ.box.plot <- uberfits.b %>% 
  filter(! hazard %in% c('Constant Hazard')) %>%
  ggplot(.,
         ## NB: this orders by the median delta.BIC across both sexes
         ## (which is why neither panel is perfectly in order)
         aes(x=forcats::fct_reorder(hazard, delta.BIC), 
             y=delta.BIC + eps, 
             group=interaction(hazard, sex))) +
  geom_boxplot() +
  scale_y_log10(labels=fancy_scientific) +
  facet_grid(~ sex) +
  ylab(expression(log[10](Delta*BIC))) +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1))
delta.bic.summ.box.plot

pdf(file.path(paper.out.dir, "delta-bic-boxplot.pdf"),
    width=9, height=4.5)
print(delta.bic.summ.box.plot)
dev.off()


########################
## CV results

CV.subsetrank.plot.dat <- uberfits.b %>%
  filter(! hazard %in% c('Constant Hazard', 'Weibull')) %>%
  group_by(data) %>%
  mutate(CV.rank.subset = dense_rank(CV.rmse.Dx)) 
  
CV.subsetrank.med <- CV.subsetrank.plot.dat %>%
  group_by(sex, hazard) %>%
  summarize(CV.rank.median = median(CV.rank.subset))

CV.subsetrank.plot.dat <- CV.subsetrank.plot.dat %>% left_join(CV.subsetrank.med)

## order hazards by median rank
CV.subsetrank.plot.dat <- CV.subsetrank.plot.dat %>%
  mutate(hazard = fct_reorder(hazard, CV.rank.median))

CV.subsetrank.med <- CV.subsetrank.med %>%
  ungroup() %>%
  mutate(hazard = fct_relevel(hazard, levels(CV.subsetrank.plot.dat$hazard)))
  
CV.subsetrank.plot <-  ggplot(CV.subsetrank.plot.dat) +
  geom_bar(aes(x=CV.rank.subset)) +
  geom_vline(aes(xintercept=CV.rank.median), 
             linetype=2,
             color='red',
             data=CV.subsetrank.med) +
  #facet_grid(sex ~ hazard)
  facet_grid(hazard ~ sex) +
  xlab("Model rank, according to CV") +
  ylab("Number of cohorts") 
#CV.subsetrank.plot

ggsave(plot=CV.subsetrank.plot,
       width=8, height=11,
       filename=file.path(paper.out.dir, "summ-CV-subsetrank.pdf"))


###################
### CV: good/bad plots

cv.goodbad.summ <- CV.subsetrank.plot.dat %>%
  group_by(sex, hazard) %>%
  summarize(
    frac.cv.top2 = mean(CV.rank.subset <= 2.0),
    frac.cv.bot2 = mean(CV.rank.subset >= 7.0)
  ) %>%
  arrange(sex, desc(frac.cv.top2))

write_csv(cv.goodbad.summ,
          path=file.path(paper.out.dir, "summ-cv-goodbad.csv"))

toplot <- cv.goodbad.summ %>%
  gather(range, fraction, -hazard, -sex) 

delta.cv.goodbad.plot <-
  ggplot(cv.goodbad.summ %>%
           filter(hazard != 'Constant Hazard')) +
  geom_point(aes(x=frac.cv.top2, y=frac.cv.bot2, shape=sex)) +
  geom_text_repel(aes(x=frac.cv.top2, y=frac.cv.bot2, label=hazard)) +
  xlim(0,1) +
  ylim(0,1) +
  xlab(expression(Fraction~of~cohorts~with~CV*top*2)) +
  ylab(expression(Fraction~of~cohorts~with~CV*bottom*2)) +
  scale_shape_manual(name="",
                     values=c('Females'=1, 'Males'=2)) +
  theme_minimal() +
  theme(legend.position="bottom") +
  guides(shape = guide_legend(override.aes = list(size=3)))
#delta.cv.goodbad.plot

pdf(file.path(paper.out.dir, "summ-CV-goodbad.pdf"),
    width=8, height=8)
print(delta.cv.goodbad.plot)
dev.off()

##############
### CV -- good vs bad plots by country

cv.goodbad.country.summ <- CV.subsetrank.plot.dat %>%
  group_by(country, sex, hazard) %>%
  summarize(
    frac.cv.top2 = mean(CV.rank.subset <= 2.0),
    frac.cv.bot2 = mean(CV.rank.subset >= 7.0)
  ) %>%
  arrange(sex, desc(frac.cv.top2))

set.seed(101)
delta.cv.goodbad.country.plot <-
  ggplot(cv.goodbad.country.summ %>%
           filter(hazard != 'Constant Hazard') %>%
           filter(country %in% c('Denmark', 'France', 'Italy', 'Netherlands', 'Sweden'))) +
  geom_point(aes(x=frac.cv.top2, y=frac.cv.bot2, shape=sex)) +
  geom_text_repel(aes(x=frac.cv.top2, y=frac.cv.bot2, label=hazard),
                  max.iter=1000,
                  #max.iter=10000,
                  size=2.5) +
  facet_grid(country ~ sex) +
  xlim(0,1) +
  ylim(0,1) +
  xlab(expression(Fraction~of~cohorts~with~CV*top*2)) +
  ylab(expression(Fraction~of~cohorts~with~CV*bottom*2)) +
  scale_shape_manual(name="",
                     values=c('Females'=1, 'Males'=2)) +
  theme_bw(base_size=10) +
  coord_equal() +
  theme(legend.position="bottom",
        strip.background=element_rect(fill=NA)) 
  

ggsave(plot=delta.cv.goodbad.country.plot,
       filename=file.path(paper.out.dir, "cv-goodbad-country.pdf"),
       height=10, width=6.5)

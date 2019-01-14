###########################################################
##  kt-forpaper-clustering.R
##
## produce the figures for the old-age mortality paper
##
## uses the output of fit-kt-data.R and kt-fit-postest.R
##
###########################################################

root.dir <- ".."

library(mortfit)
library(tidyverse)

library(gridExtra)
library(ggplot2)
library(xtable)
library(GGally)
library(cluster)
library(Rtsne)
library(gridExtra)
library(grid)

run.name <- "cohort"

out.dir <- path.expand(file.path(root.dir, "out", run.name))

out.summ.dir <- out.dir
paper.out.dir <- file.path(out.dir, "for-paper")

code.dir <- path.expand(file.path(root.dir, "code"))

## Load the data
load(file=file.path(out.summ.dir, 
                    paste0("kt-summary-", run.name, "-fits.RData")))

## load the ggplot theme...
source(file.path(root.dir, "code", "oldage_theme.R"))
theme_set(oldage_theme)

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
    mutate(hazard=
              dplyr::recode(hazard,
                            'Quadratic'='Log-Quadratic'))

uberfits.b <- uberfits.b %>%
  mutate(sex=
           dplyr::recode(sex,
                         'f'='Females',
                         'm'='Males'))

## dataset with only the relationship between delta AICs
tp <- uberfits.b %>%
  ungroup() %>%
  mutate(hazard = dplyr::recode(hazard,
                                'Constant Hazard'='ConstantHazard',
                                'Log-Quadratic'='LogQuadratic',
                                'Lynch-Brown'='LynchBrown')) %>%
  select(hazard, country, sex, year, delta.AIC) %>%
  spread(hazard, delta.AIC)

############################################################
## hierarchical clustering

## hierarchical clustering
## see https://github.com/andrie/ggdendro
library(ggdendro)
tp.pairs.hc <- hclust(dist(t(tp %>% 
                             dplyr::rename('Log-quadratic'='LogQuadratic',
                                           'Lynch-Brown'='LynchBrown') %>%
                             select(-country, -sex, -year, -Weibull, -ConstantHazard))),
                      method="average")

model.dend.plot <- ggdendrogram(dendro_data(tp.pairs.hc, type='rectangle'), 
             rotate=TRUE,
             size=4) +
  ggtitle(expression(bold(Hierarchical~clustering~of~model~Delta*AIC))) +
  xlab("Dissimilarity index")

ggsave(plot=model.dend.plot,
       width=5, height=3,
       filename=file.path(paper.out.dir, "delta-aic-model-dendrogram.pdf"))




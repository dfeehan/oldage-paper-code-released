##############################################
## oldage_theme.R
##
## ggplot2 theme for slides used in
## oldage paper
##

#oldage_theme <- function(base_size=12, base_family="",...) {
#  modifyList(theme_bw(base_size=base_size, base_family=base_family),
#             list( panel.grid.major = theme_blank(),
#                   panel.grid.minor = theme_blank()))
#}

oldage_theme <- theme_minimal(base_size=12) +
#oldage_theme <- theme_bw(base_size=12) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())

#oldage_theme <- function(base_size=12, base_family="",...) {
#  modifyList(theme_bw(base_size=base_size, base_family=base_family),
#             list( panel.grid.major = theme_blank(),
#                   panel.grid.minor = theme_blank()))
#}


## see https://github.com/hadley/ggplot2/wiki/Themes
##theme_minimal_light <- function (base_size = 12, base_family = "", ...){
##  modifyList (theme_minimal (base_size = base_size, base_family = base_family), 
##              list (axis.ticks = theme_segment (colour = "grey50"), 
##                    axis.text.x = theme_text (colour = "grey33"), 
##                    axis.text.y = theme_text (colour = "grey33")))
##}

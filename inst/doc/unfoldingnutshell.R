## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----breakfast, message = FALSE-----------------------------------------------
library(smacof)
head(breakfast)

## ----privacy, message = FALSE-------------------------------------------------
library(MPsychoR)
data(Privacy)
Privacy_rev <- 101-Privacy
head(Privacy_rev)

## ----breakunf-----------------------------------------------------------------
un_breakfast <- unfolding(breakfast)
un_breakfast

## ----breakconf, fig.width = 6, fig.height = 6---------------------------------
plot(un_breakfast, main = "Configuration Breakfast Data")

## ----privunf------------------------------------------------------------------
un_privacy <- unfolding(Privacy_rev)
un_privacy

## ----privconf, fig.width = 6, fig.height = 6----------------------------------
plot(un_privacy, main = "Unfolding Configuration Internet Privacy", 
     label.conf.rows = list(label = FALSE))

## ----ggconfu, message = FALSE, fig.width = 6, fig.height = 6------------------
library(ggplot2)
conf_items <- as.data.frame(un_privacy$conf.col)
conf_persons <- as.data.frame(un_privacy$conf.row)
p <- ggplot(conf_persons, aes(x = D1, y = D2)) 
p + geom_point(size = 0.5, colour = "gray") + coord_fixed() + 
  geom_point(aes(x = D1, y = D2), conf_items, colour = "cadetblue") + 
  geom_text(aes(x = D1, y = D2, label = rownames(conf_items)), 
            conf_items, colour = "cadetblue", vjust = -0.8) + 
  ggtitle("Unfolding Configuration Internet Privacy")


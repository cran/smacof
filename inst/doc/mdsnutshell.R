## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----ekman--------------------------------------------------------------------
library(smacof)
ekman

## ----ekmanD, eval = FALSE-----------------------------------------------------
#  ekmanD <- sim2diss(ekman, method = "confusion", to.dist = TRUE)

## ----wenchuan, message = FALSE------------------------------------------------
library(MPsychoR)
data(Wenchuan)
head(Wenchuan)

## ----wendelta-----------------------------------------------------------------
Rmat <- cor(Wenchuan, use = "pairwise.complete.obs")
Delta <- sim2diss(Rmat, method = "corr", to.dist = TRUE)
round(Delta, 2)

## ----wenmds-------------------------------------------------------------------
mds_wen1 <- mds(Delta)
mds_wen1

## ----wenmds2, fig.width = 5, fig.height = 5-----------------------------------
mds_wen2 <- mds(Delta, type = "ordinal")
plot(mds_wen1, plot.type = "Shepard", main = "Shepard Diagram (Ratio Transformation)")
plot(mds_wen2, plot.type = "Shepard", main = "Shepard Diagram (Ordinal Transformation)")

## ----wenord-------------------------------------------------------------------
mds_wen2

## ----confplot, fig.width = 6, fig.height = 6----------------------------------
plot(mds_wen2, main = "PTSD Symptoms Configuration")

## ----ggconf, message = FALSE, fig.width = 6, fig.height = 6-------------------
library(ggplot2)
conf_wen <- as.data.frame(mds_wen2$conf)
p <- ggplot(conf_wen, aes(x = D1, y = D2, label = rownames(conf_wen))) 
p + geom_point(size = 0.9) + coord_fixed() + geom_text(size = 3.5, vjust = -0.8) + 
  ggtitle("PTSD Symptoms Configuration")


## ----preliminaries, echo=FALSE, results='hide'------------------------
library("knitr")
opts_chunk$set(highlight=FALSE, prompt=TRUE, background='#FFFFFF')
options(replace.assign=TRUE, width=72, prompt="R> ")

## ----setup, include=FALSE, cache=FALSE--------------------------------
render_sweave()

## ----intelligence-----------------------------------------------------
library("smacof")
idiss <- sim2diss(intelligence[,paste0("T", 1:8)])  ## convert correlation 
fitrat <- mds(idiss)                                ## ratio MDS
fitint <- mds(idiss, type = "interval")             ## interval MDS
fitord <- mds(idiss, type = "ordinal")              ## ordinal MDS
fitspl <- mds(idiss, type = "mspline")              ## spline MDS
round(c(fitrat$stress, fitint$stress, fitord$stress, fitspl$stress), 3)

## ----shepard-plot, eval=FALSE, echo=FALSE-----------------------------
#  op <- par(mfrow = c(2,2))
#  plot(fitrat, plot.type = "Shepard", main = "Ratio MDS")
#  plot(fitint, plot.type = "Shepard", main = "Interval MDS")
#  plot(fitord, plot.type = "Shepard", main = "Ordinal MDS")
#  plot(fitspl, plot.type = "Shepard", main = "Spline MDS")
#  par(op)

## ----shepard-plot1, echo=FALSE, fig.width=8, fig.height=8-------------
op <- par(mfrow = c(2,2))
plot(fitrat, plot.type = "Shepard", main = "Ratio MDS")
plot(fitint, plot.type = "Shepard", main = "Interval MDS")
plot(fitord, plot.type = "Shepard", main = "Ordinal MDS")
plot(fitspl, plot.type = "Shepard", main = "Spline MDS")
par(op)

## ----lawlerfit--------------------------------------------------------
WishD <- sim2diss(wish, method = 7)  ## convert into dissimilarities
fitWish <- mds(WishD)                ## fit ratio MDS

## ----rstarts, cache=TRUE----------------------------------------------
stressvec <- NULL
fitbest <- NULL
set.seed(123)
for(i in 1:100) {
  fitran <- mds(WishD, init = "random")
  stressvec[i] <- fitran$stress
  if (i == 1) fitbest <- fitran
  if (stressvec[i] < fitbest$stress) fitbest <- fitran
}
round(fitbest$stress, 4)

## ----ic-plot, eval=FALSE----------------------------------------------
#  set.seed(123)
#  icWish <- icExplore(WishD, nrep = 100, returnfit = TRUE)
#  plot(icWish, main = "IC Plot Wish")

## ----ic-plot1, echo=FALSE, cache=TRUE, out.width='0.7\\textwidth'-----
set.seed(123)
icWish <- icExplore(WishD, nrep = 100, returnfit = TRUE)
plot(icWish, main = "IC Plot Wish")

## ----lawler-----------------------------------------------------------
LawlerD <- sim2diss(Lawler, to.dist = TRUE)
fitLaw <- mds(LawlerD)

## ----rstress, cache=TRUE----------------------------------------------
set.seed(123)
rstress <- randomstress(n = attr(LawlerD, "Size"), ndim = 2, nrep = 500, 
                        type = "ratio")

## ----mrstress---------------------------------------------------------
lbound <- mean(rstress) - 2*sd(rstress)
round(lbound, 4)

## ----permlaw, cache=TRUE----------------------------------------------
set.seed(123)
permLaw <- permtest(fitLaw, nrep = 500, verbose = FALSE)
permLaw

## ----wenmds-----------------------------------------------------------
data(Wenchuan, package = "MPsychoR")
Wdelta <- dist(t(Wenchuan))                                        
fitWen <- mds(Wdelta, type = "interval")        
round(fitWen$stress, 4)

## ----wenperm, cache=TRUE----------------------------------------------
set.seed(123)
permWen <- permtest(fitWen, data = Wenchuan, method.dat = "euclidean", 
                    nrep = 1000, verbose = FALSE)
permWen

## ----wenperm-plot, eval=FALSE, echo = FALSE---------------------------
#  op <- par(mfrow = c(1,2))
#  plot(permWen)
#  perm5 <- quantile(permWen$stressvec, probs = 0.05)    ## critical value
#  hist(permWen$stressvec, xlim = c(0.10, 0.40), xlab = "Stress Values",
#       main = "Wenchuan Permutation Histogram")
#  abline(v = perm5, lty = 2, col = "red", lwd = 2)
#  abline(v = fitWen$stress)
#  par(op)

## ----wenperm-plot1, echo=FALSE, fig.width=12, fig.height=6------------
op <- par(mfrow = c(1,2))
plot(permWen)
perm5 <- quantile(permWen$stressvec, probs = 0.05)    ## critical value
hist(permWen$stressvec, xlim = c(0.10, 0.40), xlab = "Stress Values", 
     main = "Wenchuan Permutation Histogram")
abline(v = perm5, lty = 2, col = "red", lwd = 2)             
abline(v = fitWen$stress)       
par(op)

## ----jackmds----------------------------------------------------------
data(Rogers, package = "MPsychoR")
RogersSub <- Rogers[,1:16]
RogersD <- dist(t(RogersSub))
fitRogers <- mds(RogersD, type = "ordinal")
jackRogers <- jackknife(fitRogers)
jackRogers

## ----jack-plot, eval=FALSE--------------------------------------------
#  plot(jackRogers)

## ----jack-plot1, echo=FALSE, out.width='0.7\\textwidth'---------------
plot(jackRogers)

## ----bootmds, cache=TRUE----------------------------------------------
set.seed(123)
bootRogers <- bootmds(fitRogers, RogersSub, method.dat = "euclidean", 
                      nrep = 500)
bootRogers

## ----boot-plot, eval=FALSE--------------------------------------------
#  plot(bootRogers)

## ----boot-plot1, echo=FALSE, out.width='0.7\\textwidth'---------------
plot(bootRogers)

## ----confmds----------------------------------------------------------
fitRogers2 <- mds(RogersD)               
confRogers <- confEllipse(fitRogers2)    

## ----cell-plot, eval=FALSE--------------------------------------------
#  plot(confRogers, eps = 0.01, ylim = c(-0.11, 0.11),
#       ell = list(lty = 1, col = "gray"))

## ----cell-plot1, echo=FALSE, out.width='0.7\\textwidth'---------------
plot(confRogers, eps = 0.01, ylim = c(-0.11, 0.11), 
     ell = list(lty = 1, col = "gray"))

## ----faceexp----------------------------------------------------------
Delta <- FaceExp
attr(Delta, "Labels") <- NULL            ## remove object labels
fitFace <- mds(Delta, type = "ordinal")  ## fit ordinal MDS

## ----facecons---------------------------------------------------------
Z <- FaceScale[, c(1,3)]                 ## select PU and TS as covariates
fitFaceC1 <- smacofConstraint(Delta, type = "ordinal", 
  constraint = "diagonal", external = Z, constraint.type = "ordinal", 
  init = fitFace$conf)

## ----trafo-plot, eval=FALSE, echo=FALSE-------------------------------
#  op <- par(mfrow = c(1,2))
#  ind1 <- order(Z[,1])
#  plot(Z[ind1, 1], fitFaceC1$external[ind1, 1], type = "s", lwd = 2,
#       xlab = "Original Scale", ylab = "Transformed Scale",
#       main = "Pleasant-Unpleasant Transformation")
#  ind2 <- order(Z[,2])
#  plot(Z[ind2, 2], fitFaceC1$external[ind2, 2], type = "s", lwd = 2,
#       xlab = "Original Scale", ylab = "Transformed Scale",
#       main = "Tension-Sleep Transformation")
#  par(op)

## ----trafo-plot1, echo=FALSE, fig.width = 12, fig.height = 5.5--------
op <- par(mfrow = c(1,2))
ind1 <- order(Z[,1])
plot(Z[ind1, 1], fitFaceC1$external[ind1, 1], type = "s", lwd = 2,
     xlab = "Original Scale", ylab = "Transformed Scale", 
     main = "Pleasant-Unpleasant Transformation")
ind2 <- order(Z[,2])
plot(Z[ind2, 2], fitFaceC1$external[ind2, 2], type = "s", lwd = 2,
     xlab = "Original Scale", ylab = "Transformed Scale", 
     main = "Tension-Sleep Transformation")
par(op)

## ----C1---------------------------------------------------------------
fitFaceC1$C 

## ----face-plot, eval=FALSE, echo=FALSE--------------------------------
#  plot(fitFaceC1, xlab = "Pleasant-Unpleasant", ylab = "Tension-Sleep",
#       main = "Face Expression (Diagonal Restriction)")

## ----face-plot1, echo=FALSE, out.width='0.65\\textwidth'--------------
plot(fitFaceC1, xlab = "Pleasant-Unpleasant", ylab = "Tension-Sleep", 
     main = "Face Expression (Diagonal Restriction)")

## ----ordcheck---------------------------------------------------------
order(Z[,1], decreasing = TRUE)  ## point order PU (x-axis)
order(Z[,2], decreasing = TRUE)  ## point order TS (y-axis) 

## ----face2a-----------------------------------------------------------
fitFaceC2 <- smacofConstraint(Delta, type = "ordinal", 
  constraint = "diagonal", external = Z, constraint.type = "interval", 
  init = fitFace$conf)
round(fitFaceC2$stress, 3)

## ----face2b-----------------------------------------------------------
fitFaceC3 <- smacofConstraint(Delta, type = "ordinal", 
  constraint = "linear", external = Z, constraint.type = "ordinal", 
  init = fitFace$conf)
round(fitFaceC3$stress, 3)
fitFaceC3$C   ## C unrestricted

## ----face2-plot, eval=FALSE, echo=2:4---------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(fitFaceC2, xlab = "Pleasant-Unpleasant", ylab = "Tension-Sleep",
#       main = "Face Expression (Interval Transformation)")
#  plot(fitFaceC3, main = "Face Expression (Linear Restriction)")
#  par(op)

## ----face2-plot1, echo=FALSE, fig.width = 12, fig.height = 6----------
op <- par(mfrow = c(1,2))
plot(fitFaceC2, xlab = "Pleasant-Unpleasant", ylab = "Tension-Sleep", 
     main = "Face Expression (Interval Transformation)")
plot(fitFaceC3, main = "Face Expression (Linear Restriction)")
par(op)

## ----carconf, message=FALSE-------------------------------------------
library("prefmod")
carconf1 <- carconf[1:100, 1:6]
head(carconf1)

## ----ordinal-unfolding------------------------------------------------
unf_ord <- unfolding(carconf1, type = "ordinal")  ## fit ordinal unfolding
unf_ord

## ----carunf-plot, eval=FALSE, echo=2:4--------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(unf_ord, main = "Unfolding Configuration Car Preferences")
#  plot(unf_ord, plot.type = "Shepard",
#       main = "Shepard Diagram Car Preferences")
#  par(op)

## ----carunf-plot1, echo=FALSE, fig.width=12, fig.height=6-------------
op <- par(mfrow = c(1,2))
plot(unf_ord, main = "Unfolding Configuration Car Preferences")
plot(unf_ord, plot.type = "Shepard", 
     main = "Shepard Diagram Car Preferences")
par(op)

## ----row-unfolding, cache = TRUE--------------------------------------
startconf <- list(unf_ord$conf.row, unf_ord$conf.col) 
unf_cond <- unfolding(carconf1, type = "ordinal", conditionality = "row", 
                      eps = 6e-5, init = startconf)
unf_cond

## ----condunf-plot, eval=FALSE, echo=2:4-------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(unf_cond, main = "Conditional Unfolding Configuration Car Preferences")
#  plot(unf_cond, plot.type = "Shepard",
#       main = "Shepard Diagram Car Preferences", col.dhat = "gray")
#  par(op)

## ----condunf-plot1, echo=FALSE, fig.width=12, fig.height=6------------
op <- par(mfrow = c(1,2))
plot(unf_cond, main = "Conditional Unfolding Configuration Car Preferences")
plot(unf_cond, plot.type = "Shepard", 
     main = "Shepard Diagram Car Preferences", col.dhat = "gray")
par(op)

## ----unfcircle, cache=TRUE--------------------------------------------
unf_vals <- unfolding(indvalues)                      ## unrestricted
unf_valsc <- unfolding(indvalues, circle = "column")  ## restricted

## ----schwartzcir-plot, eval=FALSE, echo=FALSE-------------------------
#  op <- par(mfrow = c(1,2))
#  plot(unf_vals, main = "Values Unfolding Configuration (Unrestricted)", ylim = c(-1.5, 1.1))
#  abline(h = 0, v = 0, col = "gray", lty = 2)
#  plot(unf_valsc, main = "Values Unfolding Configuration (Circular Restriction)", ylim = c(-1.5, 1.1))
#  abline(h = 0, v = 0, col = "gray", lty = 2)
#  rad <- sqrt(sum(unf_valsc$conf.col[1,]^2))
#  draw.circle(0, 0, radius = rad, border = "cadetblue", lty = 2)
#  par(op)

## ----schwartzcir-plot1, echo=FALSE, message=FALSE, fig.width=12, fig.height=6----
op <- par(mfrow = c(1,2))
plot(unf_vals, main = "Values Unfolding Configuration (Unrestricted)", ylim = c(-1.5, 1.1))
abline(h = 0, v = 0, col = "gray", lty = 2)
plot(unf_valsc, main = "Values Unfolding Configuration (Circular Restriction)", ylim = c(-1.5, 1.1))
abline(h = 0, v = 0, col = "gray", lty = 2)
rad <- sqrt(sum(unf_valsc$conf.col[1,]^2))
draw.circle(0, 0, radius = rad, border = "cadetblue", lty = 2)
par(op)

## ----valcirc----------------------------------------------------------
tuv <- matrix(NA, nrow = ncol(PVQ40agg), ncol = 2) 
alpha <- -360/10
for (i in 1:10){ 
  alpha <- alpha+360/10
  tuv[i,1] <- cos(alpha*pi/180)
  tuv[i,2] <- sin(alpha*pi/180) 
}

## ----extunf-----------------------------------------------------------
delta <- (max(PVQ40agg) + 1) - PVQ40agg 
unf_pvq <- unfolding(delta, type = "ordinal")        ## unconstrained
unf_pvqext <- unfolding(delta, type = "ordinal", 
    fixed = "column", fixed.coord = tuv)             ## constrained

## ----unfext-plot1, echo=FALSE, message=FALSE, fig.width=12, fig.height=6----
op <- par(mfrow = c(1,2))
plot(unf_pvq, label.conf.columns = list(cex = 1), ylim = c(-1.1, 1.1), main = "PVQ Unrestricted")
plot(unf_pvqext, label.conf.columns = list(cex = 1), ylim = c(-1.1, 1.1), main = "PVQ Fixed Value Circle")
par(op)

## ----vmu--------------------------------------------------------------
fit_vmu <- vmu(PVQ40agg)        ## VMU fit
fit_vmu

## ----vmubi-plot, eval=FALSE-------------------------------------------
#  plot(fit_vmu, cex = c(1, 0.7))

## ----vmubi-plot1, echo=FALSE, message=FALSE, out.width='0.7\\textwidth'----
plot(fit_vmu, cex = c(1, 0.7)) 

## ----plato------------------------------------------------------------
PlatoD <- dist(t(Plato7))
fitPlato <- uniscale(PlatoD, verbose = FALSE)  ## unidimensional fit
round(sort(fitPlato$conf), 3)

## ---------------------------------------------------------------------
gopD <- gravity(GOPdtm)$gravdiss

## ----gravity-plot, eval=FALSE-----------------------------------------
#  fitGOP <- mds(gopD, type = "ordinal")    ## 2D ordinal MDS
#  plot(fitGOP, plot.type = "bubbleplot", bubscale = 20)

## ----gravity-plot1, echo=FALSE, out.width='0.7\\textwidth'------------
fitGOP <- mds(gopD, type = "ordinal")    ## 2D ordinal MDS
plot(fitGOP, plot.type = "bubbleplot", bubscale = 20)

## ----drift-plot, eval=FALSE-------------------------------------------
#  morseDrift <- driftVectors(morse2, type = "ordinal")
#  plot(morseDrift, main = "Drift Vectors Morse Codes")

## ----drift-plot1, echo=FALSE, out.width='0.7\\textwidth'--------------
morseDrift <- driftVectors(morse2, type = "ordinal") 
plot(morseDrift, main = "Drift Vectors Morse Codes")

## ----maryam-----------------------------------------------------------
artD <- sim2diss(VaziriXu$artificialR)
fitart <- mds(artD, type = "ordinal")    ## artificial condition
realD <- sim2diss(VaziriXu$realR)
fitnat <- mds(realD, type = "ordinal")   ## natural condition  

## ----maryam-plot, echo=FALSE, eval=FALSE------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(fitart, main = "Configuration Artificial")
#  plot(fitnat, main = "Configuration Natural")
#  par(op)

## ----maryam-plot1, echo=FALSE, fig.width=12, fig.height=6-------------
op <- par(mfrow = c(1,2))
plot(fitart, main = "Configuration Artificial")
plot(fitnat, main = "Configuration Natural")
par(op)

## ----procnat----------------------------------------------------------
fitproc <- Procrustes(fitart$conf, fitnat$conf)
fitproc

## ----proc-plot, eval=FALSE--------------------------------------------
#  plot(fitproc, legend = list(labels = c("artificial", "natural")))

## ----proc-plot1, echo=FALSE, out.width='0.6\\textwidth'---------------
plot(fitproc, legend = list(labels = c("artificial", "natural")))

## ----stress0----------------------------------------------------------
stress0(rectangles, init = rect_constr)  ## stress value 0 iterations

## ----proc2------------------------------------------------------------
fitrect <- mds(rectangles, type = "ordinal", init = rect_constr)
round(fitrect$stress, 3)
procRect <- Procrustes(rect_constr, fitrect$conf)
round(procRect$congcoef, 3)        ## congruence coefficient

## ----proc2-plot, eval=FALSE-------------------------------------------
#  plot(procRect, xlim = c(2, 7),
#       legend = list(labels = c("theoretical", "observed")))

## ----proc2-plot1, echo=FALSE, out.width='0.6\\textwidth'--------------
plot(procRect, xlim = c(2, 7), 
     legend = list(labels = c("theoretical", "observed")))

## ----simplebi---------------------------------------------------------
ext <- data.frame(PU = FaceScale[,1])
fitFace <- mds(FaceExp, type = "ordinal")    ## ordinal MDS fit
biFace <- biplotmds(fitFace, extvar = ext)   ## map external variable
coef(biFace)
round(biFace$R2vec, 3)

## ----simplebi-plot, eval=FALSE, echo=c(1, 3:9)------------------------
#  library("calibrate")
#  op <- par(mfrow = c(1,2))
#  plot(biFace, main = "Biplot Vector", vecscale = 0.8,
#       vec.conf = list(col = "brown"), xlim = c(-1.5, 1.5))
#  plot(fitFace, main = "Biplot Axis", xlim = c(-1.5, 1.5))
#  abline(h = 0, v = 0, lty = 2, col = "gray")
#  calFace <- calibrate(biFace$coef[, 1], ext$PU, tm = seq(-5, 6, by = 1),
#                       Fr = fitFace$conf, dp = TRUE, axiscol = "brown",
#                       axislab = "PU", labpos = 3, verb = FALSE)
#  par(op)

## ----simplebi-plot1, echo=FALSE, message=FALSE, fig.width = 12, fig.height = 6----
library("calibrate")
op <- par(mfrow = c(1,2))
plot(biFace, main = "Biplot Vector", vecscale = 0.8, 
     vec.conf = list(col = "brown"), xlim = c(-1.5, 1.5))
plot(fitFace, main = "Biplot Axis", xlim = c(-1.5, 1.5))
abline(h = 0, v = 0, lty = 2, col = "gray")
calFace <- calibrate(biFace$coef[, 1], ext$PU, tm = seq(-5, 6, by = 1), 
                     Fr = fitFace$conf, dp = TRUE, axiscol = "brown", 
                     axislab = "PU", labpos = 3, verb = FALSE)
par(op)

## ----neuralbi---------------------------------------------------------
library("MPsychoR")
data(NeuralActivity)       ## dissimilarity matrices
data(NeuralScales)         ## external variables
NeuralD <- Reduce("+", NeuralActivity)/length(NeuralActivity)
fitNeural <- mds(NeuralD, type = "mspline")
biNeural <- biplotmds(fitNeural, NeuralScales[,1:8])
round(biNeural$R2vec, 3)

## ----neuralbi-plot, eval=FALSE----------------------------------------
#  plot(biNeural, col = "gray", label.conf = list(col = "gray", cex = 0.7),
#       vec.conf = list(col = "brown"), main = "Neural Biplot",
#       xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))

## ----neuralbi-plot1, echo=FALSE, out.width='0.7\\textwidth'-----------
plot(biNeural, col = "gray", label.conf = list(col = "gray", cex = 0.7),
     vec.conf = list(col = "brown"), main = "Neural Biplot",
     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))


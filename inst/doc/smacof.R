## ----preliminaries, echo=FALSE, results='hide'------------------------
library("knitr")
opts_chunk$set(highlight=FALSE, prompt=TRUE, background='#FFFFFF')
options(replace.assign=TRUE, width=72, prompt="R> ")

## ----setup, include=FALSE, cache=FALSE--------------------------------
render_sweave()

## ----intelligence, message=FALSE--------------------------------------
library("smacof")
idiss <- sim2diss(intelligence[,paste0("T", 1:8)])  
fitrat <- mds(idiss)                                
fitint <- mds(idiss, type = "interval")             
fitord <- mds(idiss, type = "ordinal")              
fitspl <- mds(idiss, type = "mspline")              
round(c(fitrat$stress, fitint$stress, fitord$stress, fitspl$stress), 3)

## ----shepard-plot, eval=FALSE, echo=FALSE-----------------------------
#  op <- par(mfrow = c(2,2))
#  plot(fitrat, plot.type = "Shepard", main = "Ratio MDS")
#  plot(fitint, plot.type = "Shepard", main = "Interval MDS")
#  plot(fitord, plot.type = "Shepard", main = "Ordinal MDS")
#  plot(fitspl, plot.type = "Shepard", main = "Spline MDS")
#  par(op)

## ----shepard-plot1, echo=FALSE, fig.width=8, fig.height=8, out.width='0.8\\textwidth'----
op <- par(mfrow = c(2,2))
plot(fitrat, plot.type = "Shepard", main = "Ratio MDS")
plot(fitint, plot.type = "Shepard", main = "Interval MDS")
plot(fitord, plot.type = "Shepard", main = "Ordinal MDS")
plot(fitspl, plot.type = "Shepard", main = "Spline MDS")
par(op)

## ----lawlerfit--------------------------------------------------------
WishD <- sim2diss(wish, method = 7)  
fitWish <- mds(WishD)                

## ----rstarts, cache=TRUE----------------------------------------------
set.seed(123)
stressvec <- rep(NA, 100)
fitbest <- mds(WishD, init = "random")
stressvec[1] <- fitbest$stress
for(i in 2:100) {
  fitran <- mds(WishD, init = "random")
  stressvec[i] <- fitran$stress
  if (fitran$stress < fitbest$stress) fitbest <- fitran
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
rstress <- randomstress(n = 9, ndim = 2, nrep = 500, type = "ratio")

## ----mrstress---------------------------------------------------------
bound <- mean(rstress) - 2*sd(rstress)
round(bound, 3)

## ----permlaw, cache=TRUE----------------------------------------------
set.seed(123)
permLaw <- permtest(fitLaw, nrep = 500, verbose = FALSE)
permLaw

## ----wenmds-----------------------------------------------------------
data("Wenchuan", package = "MPsychoR")
Wdelta <- dist(t(Wenchuan))                                        
fitWen <- mds(Wdelta, type = "interval")        
round(fitWen$stress, 3)

## ----wenperm, cache=TRUE----------------------------------------------
set.seed(123)
permWen <- permtest(fitWen, data = Wenchuan, method.dat = "euclidean", 
  nrep = 1000, verbose = FALSE)
permWen

## ----wenperm-plot, eval=FALSE, echo = FALSE---------------------------
#  op <- par(mfrow = c(1,2))
#  plot(permWen, xlim = c(0.19, 0.30))
#  perm5 <- quantile(permWen$stressvec, probs = 0.05)
#  hist(permWen$stressvec, xlim = c(0.10, 0.40), xlab = "Stress Values",
#       main = "Wenchuan Permutation Histogram")
#  abline(v = perm5, lty = 2, col = "gray", lwd = 2)
#  abline(v = fitWen$stress, col = "gray")
#  par(op)

## ----wenperm-plot1, echo=FALSE, fig.width=12, fig.height=6------------
op <- par(mfrow = c(1,2))
plot(permWen, xlim = c(0.19, 0.30))
perm5 <- quantile(permWen$stressvec, probs = 0.05)    
hist(permWen$stressvec, xlim = c(0.10, 0.40), xlab = "Stress Values", 
     main = "Wenchuan Permutation Histogram")
abline(v = perm5, lty = 2, col = "gray", lwd = 2)             
abline(v = fitWen$stress, col = "gray")       
par(op)

## ----jackmds----------------------------------------------------------
data("Rogers", package = "MPsychoR")
RogersSub <- Rogers[,1:16]
RogersD <- dist(t(RogersSub))
fitRogers <- mds(RogersD, type = "ordinal")
jackRogers <- jackmds(fitRogers)
jackRogers

## ----jack-plot, eval=FALSE--------------------------------------------
#  plot(jackRogers, legend = TRUE, cex.legend = 0.8, inset = c(-0.3, 0))

## ----jack-plot1, echo=FALSE, out.width='0.8\\textwidth'---------------
plot(jackRogers, legend = TRUE, cex.legend = 0.8, inset = c(-0.3, 0))

## ----bootmds1, eval=FALSE---------------------------------------------
#  set.seed(123)
#  bootRogers <- bootmds(fitRogers, RogersSub, method.dat = "euclidean",
#    nrep = 500)
#  bootRogers

## ----boot-plot, eval=FALSE--------------------------------------------
#  plot(bootRogers)

## ----boot-plot1, echo=FALSE, out.width='0.7\\textwidth'---------------
plot(bootRogers)

## ----confmds----------------------------------------------------------
fitRogers2 <- mds(RogersD)               
confRogers <- confEllipse(fitRogers2)    

## ----cell-plot, eval=FALSE--------------------------------------------
#  plot(confRogers, eps = 0.01, ylim = c(-0.11, 0.11),
#    ell = list(lty = 1, col = "gray"))

## ----cell-plot1, echo=FALSE, out.width='0.7\\textwidth'---------------
plot(confRogers, eps = 0.01, ylim = c(-0.11, 0.11), 
  ell = list(lty = 1, col = "gray"))

## ----simplebi---------------------------------------------------------
fitFace <- mds(FaceExp, type = "ordinal")    
ext <- data.frame(PU = FaceScale[,1])
biFace <- biplotmds(fitFace, extvar = ext)   
coef(biFace)

## ----biaxprep---------------------------------------------------------
library("calibrate")
biFace2 <- biplotmds(fitFace, extvar = ext, scale = FALSE)           
coef(biFace2)  
PUc <- scale(ext, scale = FALSE)
tm <- seq(floor(min(ext)), ceiling(max(ext)), by = 1)           
tmc_PU <- tm - mean(tm)  
X <- fitFace$conf

## ----simplebi-plot, eval=FALSE----------------------------------------
#  plot(biFace, main = "Biplot Vector Representation", vecscale = 0.8,
#    xlim = c(-1.5, 1.5), vec.conf = list(col = "brown"), pch = 20, cex = 0.5)
#  plot(fitFace, main = "Biplot Axis Representation", xlim = c(-1.5, 1.5))
#  abline(h = 0, v = 0, lty = 2, col = "gray")
#  calPU <- calibrate(coef(biFace2), PUc, tm = tmc_PU, tmlab = tm, Fr = X,
#    dp = TRUE, axiscol = "brown", axislab = "PU", labpos = 3, verb = FALSE)

## ----simplebi-plot1, echo=FALSE, message=FALSE, out.width='0.6\\textwidth'----
plot(biFace, main = "Biplot Vector Representation", vecscale = 0.8, 
  xlim = c(-1.5, 1.5), vec.conf = list(col = "brown"), pch = 20, cex = 0.5)
plot(fitFace, main = "Biplot Axis Representation", xlim = c(-1.5, 1.5))
abline(h = 0, v = 0, lty = 2, col = "gray")
calPU <- calibrate(coef(biFace2), PUc, tm = tmc_PU, tmlab = tm, Fr = X,
  dp = TRUE, axiscol = "brown", axislab = "PU", labpos = 3, verb = FALSE)

## ----neuralbi---------------------------------------------------------
library("MPsychoR")
data("NeuralActivity")       
data("NeuralScales")         
NeuralD <- Reduce("+", NeuralActivity)/length(NeuralActivity)
fitNeural <- mds(NeuralD, type = "mspline")
biNeural <- biplotmds(fitNeural, NeuralScales[,1:8])
round(biNeural$R2vec, 3)

## ----neuralbi-plot, eval=FALSE, echo=FALSE----------------------------
#  plot(biNeural, col = "gray", label.conf = list(col = "gray", cex = 0.7),
#    vec.conf = list(col = "brown"), main = "Neural Biplot",
#    xlim = c(-1.55, 1.55), ylim = c(-1.55, 1.55))

## ----neuralbi-plot1, echo=FALSE, out.width='0.7\\textwidth'-----------
plot(biNeural, col = "gray", label.conf = list(col = "gray", cex = 0.7),
  vec.conf = list(col = "brown"), main = "Neural Biplot",
  xlim = c(-1.55, 1.55), ylim = c(-1.55, 1.55))

## ----facecons_int-----------------------------------------------------
fitFace <- mds(FaceExp, type = "ordinal")  
Z <- FaceScale[, c(1,3)]                 
fitFaceC1 <- smacofConstraint(FaceExp, type = "ordinal", 
  constraint = "diagonal", external = Z, constraint.type = "interval", 
  init = fitFace$conf)
round(fitFaceC1$C, 3)

## ----faceint-plot, eval=FALSE, echo=FALSE-----------------------------
#  Zc <- scale(Z, scale = FALSE)
#  X1 <- fitFaceC1$conf
#  tm <- seq(1, 9, by = 1)
#  tmc_PU <- tm - mean(Z[,1])
#  tmc_TS <- tm - mean(Z[,2])
#  op <- par("xpd" = TRUE, mar = c(5, 4, 7, 5) + 0.1)
#  plot(fitFaceC1, main = "Constrained MDS (Interval)", xlim = c(-1.4, 1.4), ylim = c(-1, 1))
#  calPU <- calibrate(c(1,0), Zc[,"PU"], tm = tmc_PU, tmlab = tm, Fr = X1, dp = FALSE,
#    axiscol = "brown", axislab = "PU", labpos = 3, verb = FALSE, shiftvec = c(0, par("usr")[4]),
#    shiftfactor = 1, where = 2, laboffset = c(-0.2, 0.15), cex.axislab = 0.9, reverse = TRUE)
#  calTS <- calibrate(c(0,1), Zc[,"TS"], tm = tmc_TS, tmlab = tm, Fr = X1, dp = FALSE,
#    axiscol = "brown", axislab = "TS", labpos = 3, verb = FALSE, shiftvec = c(par("usr")[2], 0),
#    shiftfactor = 1, where = 2, laboffset = c(0.2, -0.25), cex.axislab = 0.9)
#  par(op)

## ----faceint-plot1, echo=FALSE, out.width='0.65\\textwidth'-----------
Zc <- scale(Z, scale = FALSE)
X1 <- fitFaceC1$conf
tm <- seq(1, 9, by = 1)
tmc_PU <- tm - mean(Z[,1])
tmc_TS <- tm - mean(Z[,2])
op <- par("xpd" = TRUE, mar = c(5, 4, 7, 5) + 0.1)
plot(fitFaceC1, main = "Constrained MDS (Interval)", xlim = c(-1.4, 1.4), ylim = c(-1, 1))
calPU <- calibrate(c(1,0), Zc[,"PU"], tm = tmc_PU, tmlab = tm, Fr = X1, dp = FALSE, 
  axiscol = "brown", axislab = "PU", labpos = 3, verb = FALSE, shiftvec = c(0, par("usr")[4]), 
  shiftfactor = 1, where = 2, laboffset = c(-0.2, 0.15), cex.axislab = 0.9, reverse = TRUE)
calTS <- calibrate(c(0,1), Zc[,"TS"], tm = tmc_TS, tmlab = tm, Fr = X1, dp = FALSE, 
  axiscol = "brown", axislab = "TS", labpos = 3, verb = FALSE, shiftvec = c(par("usr")[2], 0), 
  shiftfactor = 1, where = 2, laboffset = c(0.2, -0.25), cex.axislab = 0.9)
par(op)

## ----fitface_ord------------------------------------------------------
fitFaceC2 <- smacofConstraint(FaceExp, type = "ordinal", 
  constraint = "diagonal", external = Z, constraint.type = "ordinal", 
  init = fitFace$conf)
round(fitFaceC2$C, 3)

## ----trafo-plot, eval=FALSE, echo=FALSE-------------------------------
#  op <- par(mfrow = c(1,2))
#  ind1 <- order(Z[,1])
#  plot(Z[ind1, 1], fitFaceC2$external[ind1, 1], type = "s", lwd = 2,
#    xlab = "Original Scale", ylab = "Transformed Scale",
#    main = "Pleasant-Unpleasant Transformation")
#  ind2 <- order(Z[,2])
#  plot(Z[ind2, 2], fitFaceC2$external[ind2, 2], type = "s", lwd = 2,
#    xlab = "Original Scale", ylab = "Transformed Scale",
#    main = "Tension-Sleep Transformation")
#  par(op)

## ----trafo-plot1, echo=FALSE, fig.width = 12, fig.height = 5----------
op <- par(mfrow = c(1,2))
ind1 <- order(Z[,1])
plot(Z[ind1, 1], fitFaceC2$external[ind1, 1], type = "s", lwd = 2,
  xlab = "Original Scale", ylab = "Transformed Scale", 
  main = "Pleasant-Unpleasant Transformation")
ind2 <- order(Z[,2])
plot(Z[ind2, 2], fitFaceC2$external[ind2, 2], type = "s", lwd = 2,
  xlab = "Original Scale", ylab = "Transformed Scale", 
  main = "Tension-Sleep Transformation")
par(op)

## ----faceord-plot, eval=FALSE, echo=FALSE-----------------------------
#  X2 <- fitFaceC2$conf
#  Zhat <- fitFaceC2$external
#  tm <- seq(-0.8, 0.8, by = 0.2)
#  op <- par("xpd" = TRUE, mar = c(5, 4, 7, 5) + 0.1)
#  plot(fitFaceC2, main = "Constrained MDS (Ordinal)", xlim = c(-1.4, 1.4), ylim = c(-1, 1))
#  calPU <- calibrate(c(1,0), Zhat[,"PU"], tm = tm, tmlab = tm, Fr = X2, dp = FALSE,
#    axiscol = "brown", axislab = "PU transformed", labpos = 3, verb = FALSE, shiftvec = c(0, par("usr")[4]),
#    shiftfactor = 1, where = 2, laboffset = c(-0.2, 0.15), cex.axislab = 0.9, reverse = TRUE)
#  calPU <- calibrate(c(0,1), Zhat[,"TS"], tm = tm, tmlab = tm, Fr = X2, dp = FALSE,
#    axiscol = "brown", axislab = "TS transformed", labpos = 3, verb = FALSE, shiftvec = c(par("usr")[2], 0),
#    shiftfactor = 1, where = 2, laboffset = c(0.25, -0.2), cex.axislab = 0.9)
#  par(op)

## ----faceord-plot1, echo=FALSE, out.width='0.65\\textwidth'-----------
X2 <- fitFaceC2$conf
Zhat <- fitFaceC2$external
tm <- seq(-0.8, 0.8, by = 0.2)
op <- par("xpd" = TRUE, mar = c(5, 4, 7, 5) + 0.1)
plot(fitFaceC2, main = "Constrained MDS (Ordinal)", xlim = c(-1.4, 1.4), ylim = c(-1, 1))
calPU <- calibrate(c(1,0), Zhat[,"PU"], tm = tm, tmlab = tm, Fr = X2, dp = FALSE,
  axiscol = "brown", axislab = "PU transformed", labpos = 3, verb = FALSE, shiftvec = c(0, par("usr")[4]), 
  shiftfactor = 1, where = 2, laboffset = c(-0.2, 0.15), cex.axislab = 0.9, reverse = TRUE)
calPU <- calibrate(c(0,1), Zhat[,"TS"], tm = tm, tmlab = tm, Fr = X2, dp = FALSE,
  axiscol = "brown", axislab = "TS transformed", labpos = 3, verb = FALSE, shiftvec = c(par("usr")[2], 0),
  shiftfactor = 1, where = 2, laboffset = c(0.25, -0.2), cex.axislab = 0.9)
par(op)

## ----face2b-----------------------------------------------------------
Z <- FaceScale
fitFaceC3 <- smacofConstraint(FaceExp, type = "ordinal",
  constraint = "unrestricted", external = Z, constraint.type = "interval",
  init = fitFace$conf)
round(fitFaceC3$C, 3)

## ----face2-plot, eval=FALSE, echo=FALSE-------------------------------
#  Zc <- scale(Z, scale = FALSE)
#  X <- fitFaceC3$conf
#  biFaceZ <- biplotmds(fitFaceC3, extvar = Z, scale = FALSE)
#  betas <- coef(biFaceZ)
#  tm <- seq(1, 9, by = 1)
#  tmc_PU <- tm - mean(Z[,1])
#  tmc_AR <- tm - mean(Z[,2])
#  tmc_TS <- tm - mean(Z[,3])
#  plot(fitFaceC3, main = "Constraint MDS Biplot", xlim = c(-1.3, 1.3), ylim = c(-1, 1))
#  abline(h = 0, v = 0, lty = 2, col = "gray")
#  calPU <- calibrate(betas[,"PU"], Zc[,"PU"], tm = tmc_PU, tmlab = tm, Fr = X,
#    dp = FALSE, axiscol = "brown", axislab = "PU", labpos = 3, verb = FALSE)
#  calAR <- calibrate(betas[,"AR"], Zc[,"AR"], tm = tmc_AR, tmlab = tm, Fr = X,
#    dp = FALSE, axiscol = "brown", axislab = "AR", labpos = 3, verb = FALSE)
#  calTS <- calibrate(betas[,"TS"], Zc[,"TS"], tm = tmc_TS, tmlab = tm, Fr = X,
#    dp = FALSE, axiscol = "brown", axislab = "TS", labpos = 3, verb = FALSE)

## ----face2-plot1, echo=FALSE, out.width='0.7\\textwidth'--------------
Zc <- scale(Z, scale = FALSE)
X <- fitFaceC3$conf
biFaceZ <- biplotmds(fitFaceC3, extvar = Z, scale = FALSE)
betas <- coef(biFaceZ)
tm <- seq(1, 9, by = 1)
tmc_PU <- tm - mean(Z[,1])
tmc_AR <- tm - mean(Z[,2])
tmc_TS <- tm - mean(Z[,3])
plot(fitFaceC3, main = "Constraint MDS Biplot", xlim = c(-1.3, 1.3), ylim = c(-1, 1))
abline(h = 0, v = 0, lty = 2, col = "gray")
calPU <- calibrate(betas[,"PU"], Zc[,"PU"], tm = tmc_PU, tmlab = tm, Fr = X,
  dp = FALSE, axiscol = "brown", axislab = "PU", labpos = 3, verb = FALSE)
calAR <- calibrate(betas[,"AR"], Zc[,"AR"], tm = tmc_AR, tmlab = tm, Fr = X,
  dp = FALSE, axiscol = "brown", axislab = "AR", labpos = 3, verb = FALSE)
calTS <- calibrate(betas[,"TS"], Zc[,"TS"], tm = tmc_TS, tmlab = tm, Fr = X,
  dp = FALSE, axiscol = "brown", axislab = "TS", labpos = 3, verb = FALSE)

## ----carconf, message=FALSE-------------------------------------------
library("prefmod")
carconf1 <- carconf[1:100, 1:6]
head(carconf1)

## ----ordinal-unfolding------------------------------------------------
unf_ord <- unfolding(carconf1, type = "ordinal") 
unf_ord

## ----carunf-plot, eval=FALSE, echo=2:5--------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(unf_ord, main = "Unfolding Configuration Car Preferences",
#    ylim = c(-1, 1))
#  plot(unf_ord, plot.type = "Shepard",
#    main = "Shepard Diagram Car Preferences")
#  par(op)

## ----carunf-plot1, echo=FALSE, fig.width=12, fig.height=6-------------
op <- par(mfrow = c(1,2))
plot(unf_ord, main = "Unfolding Configuration Car Preferences", 
  ylim = c(-1, 1))
plot(unf_ord, plot.type = "Shepard", 
  main = "Shepard Diagram Car Preferences")
par(op)

## ----row-unfolding, cache = TRUE--------------------------------------
startconf <- list(unf_ord$conf.row, unf_ord$conf.col) 
unf_cond <- unfolding(carconf1, type = "ordinal", conditionality = "row", 
  eps = 6e-5, init = startconf)
unf_cond

## ----condunf-plot, eval=FALSE, echo=2:5-------------------------------
#  op <- par(mfrow = c(1,2))
#  plot(unf_cond, main = "Conditional Unfolding Configuration Car Preferences")
#  plot(unf_cond, plot.type = "Shepard",
#    main = "Shepard Diagram Car Preferences", col.dhat = "gray",
#    xlim = c(0.9, 6.1))
#  par(op)

## ----condunf-plot1, echo=FALSE, fig.width=12, fig.height=6------------
op <- par(mfrow = c(1,2))
plot(unf_cond, main = "Conditional Unfolding Configuration Car Preferences")
plot(unf_cond, plot.type = "Shepard", 
  main = "Shepard Diagram Car Preferences", col.dhat = "gray", 
  xlim = c(0.9, 6.1))
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
unf_pvq <- unfolding(delta, type = "ordinal")        
unf_pvqext <- unfolding(delta, type = "ordinal", fixed = "column", 
  fixed.coord = tuv)             

## ----unfext-plot1, echo=FALSE, message=FALSE, fig.width=12, fig.height=6----
op <- par(mfrow = c(1,2))
plot(unf_pvq, label.conf.columns = list(cex = 1), ylim = c(-1.1, 1.1), main = "PVQ Unrestricted")
plot(unf_pvqext, label.conf.columns = list(cex = 1), ylim = c(-1.1, 1.1), main = "PVQ Fixed Value Circle")
par(op)

## ----unfcircle, cache=TRUE--------------------------------------------
unf_vals <- unfolding(indvalues)                      
unf_valsc <- unfolding(indvalues, circle = "column")  

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

## ----vmu--------------------------------------------------------------
fit_vmu <- vmu(PVQ40agg)        
fit_vmu

## ----vmubi-plot, eval=FALSE-------------------------------------------
#  plot(fit_vmu, col = c("black", "coral3"), cex = c(1, 0.7))

## ----vmubi-plot1, echo=FALSE, message=FALSE, out.width='0.7\\textwidth'----
plot(fit_vmu, col = c("black", "coral3"), cex = c(1, 0.7)) 

## ----plato------------------------------------------------------------
PlatoD <- dist(t(Plato7))
fitPlato <- uniscale(PlatoD, verbose = FALSE)  
round(sort(fitPlato$conf), 3)

## ---------------------------------------------------------------------
gopD <- gravity(GOPdtm)$gravdiss

## ----gravity-plot, eval=FALSE-----------------------------------------
#  fitGOP <- mds(gopD, type = "ordinal")
#  plot(fitGOP, plot.type = "bubbleplot", bubscale = 20)

## ----gravity-plot1, echo=FALSE, out.width='0.7\\textwidth'------------
fitGOP <- mds(gopD, type = "ordinal")    
plot(fitGOP, plot.type = "bubbleplot", bubscale = 20)

## ----drift-plot, eval=FALSE-------------------------------------------
#  morseDrift <- driftVectors(morse2, type = "ordinal")
#  plot(morseDrift, main = "Drift Vectors Morse Codes",
#    col.drift = "darkgray")

## ----drift-plot1, echo=FALSE, out.width='0.7\\textwidth'--------------
morseDrift <- driftVectors(morse2, type = "ordinal") 
plot(morseDrift, main = "Drift Vectors Morse Codes", 
  col.drift = "darkgray")

## ----maryam-----------------------------------------------------------
artD <- sim2diss(VaziriXu$artificialR)
fitart <- mds(artD, type = "ordinal")    
realD <- sim2diss(VaziriXu$realR)
fitnat <- mds(realD, type = "ordinal")   

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
stress0(rectangles, init = rect_constr)$stress  

## ----proc2-plot, eval=FALSE, echo=TRUE--------------------------------
#  fitrect <- mds(rectangles, type = "ordinal", init = rect_constr)
#  procRect <- Procrustes(rect_constr, fitrect$conf)
#  plot(procRect, legend = list(labels = c("theoretical", "observed")),
#    xlim = c(2, 7))

## ----proc2-plot1, echo=FALSE, out.width='0.6\\textwidth'--------------
fitrect <- mds(rectangles, type = "ordinal", init = rect_constr)
procRect <- Procrustes(rect_constr, fitrect$conf)
plot(procRect, legend = list(labels = c("theoretical", "observed")),
  xlim = c(2, 7))     


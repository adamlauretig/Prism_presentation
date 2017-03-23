(library(GISTools))
(library(rgdal))
(library(spdep))
(library(sp))
(library(spatstat))
file_path <- "/Users/adamlauretig/data/prism_stuff/prism_presentation/NYAIDS_data"
ny <- readOGR(dsn = file_path, layer = "NYAIDS", verbose=FALSE)
ny <- ny[-c(38),]

nygal <- poly2nb(ny) #Create the neighborhood object
nyQ1.gal <- nb2listw(nygal, zero.policy=T) #Create the weights object
lrate <- log1p(ny$Rate1000)

reg1 <- lm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH, data = ny)
summary(reg1)

W.eig <- eigenw(nyQ1.gal, quiet=NULL)
ny.lag.eig <- lagsarlm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH,
 data=ny, nyQ1.gal, method="eigen", quiet=TRUE, control=(pre_eig=W.eig), zero.policy = TRUE) # How to use pre-calculated eigenvalues

summary(ny.lag.eig)
w_mat <- nb2mat(nygal, zero.policy=T)
w_mat <- ifelse(w_mat != 0, 1, 0)
library(CARBayes) #version 4.7
set.seed(614)
car1 <- S.CARleroux(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH, "gaussian", 
            data=ny@data, W = w_mat,  burnin=20000, n.sample=100000, thin=10, 
            prior.nu2=NULL, prior.tau2=NULL, fix.rho=FALSE, rho=NULL, 
            verbose=FALSE)

car2 <- S.CARleroux(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH, "gaussian", 
                    data=ny@data, W = w_mat,  burnin=20000, n.sample=100000, thin=10, 
                    prior.nu2=c(10, 1), prior.tau2=c(10, 1), fix.rho=FALSE, rho=NULL, 
                    verbose=FALSE)

car2$summary.results

library(GWmodel) # v. 2.0-2
coordinates(ny)
dist_mat <- gw.dist(dp.locat = coordinates(ny), rp.locat =coordinates(ny), focus = 0, p = 2)
bw1 <- bw.gwr(lrate ~ Gini, approach="aic",adaptive=TRUE, data=ny, kernel = "gaussian", dMat=dist_mat)
gw.model1 <- gwr.basic(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH, 
  bw = bw1, data = ny, kernel = "gaussian", adaptive = TRUE, p = 2, 
  longlat = FALSE)
gw.model1
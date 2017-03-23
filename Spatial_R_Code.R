#####################################
# Workshop on Spatial Statistics in R
# Adam Lauretig
# March 24, 2017
#####################################

rm(list = ls())
options(stringsAsFactors = FALSE)

# This is the basic package for working with spatial data in R
install.packages("sp")

# This is a package from Brunsdon & Comber's book, good for quick visualizations 
install.packages("GISTools")  

# This is the package for loading shapefiles
install.packages("rgdal")

# This is the package for creating weights matrices
install.packages("spdep")

# This is the package for basic spatial stats
install.packages("spatstat")

# This is the package for Conditional Autoregressive models
install.packages("CARBayes")

# This is the package for Geographically Weighted Regression models
install.packages("GWmodel")


#####################################
library(sp)
library(spdep)
library(spatstat)
library(rgdal)
library(GISTools)

#############
# Unfortunately, we can't use the "~" here (readOGR doesn't like it)
# You'll need to change this for your computer, replace
file_path <- "/Users/adamlauretig/data/Prism_presentation/NYAIDS_data"
#############

# Loading the shapefile
ny <- readOGR(dsn = file_path, layer = "NYAIDS", verbose=FALSE) 

## Dropping the polygon with no neighbors (not best practice, but 
## useful in this case)
ny <- ny[-c(38),]


# A quick plot of our DV
choropleth(ny, log1p(ny$Rate1000))
############
### The Moran's I test
############

# Create the neighborhood object
nygal <- poly2nb(ny) 
# Create the Weights matrix
nyQ1.gal <- nb2listw(nygal, zero.policy=T) 
# Run the Test
moran.test(log1p(ny$Rate1000), nyQ1.gal, alternative="two.sided", zero.policy=T)

############
### The LISA plot and Map for Clusters
############

# Standardizing our values
locm <- localmoran(log1p(ny$Rate1000), nyQ1.gal, zero.policy=T, 
  alternative="two.sided")
Scrate <- (log1p(ny$Rate1000) - mean(log1p(ny$Rate1000)))/sd(log1p(ny$Rate1000))
lagScrate <- lag.listw(nyQ1.gal, Scrate, zero.policy=T)

# Plotting out the z-scores and lagged z-scores
plot(Scrate, lagScrate, col="purple", xlab="AIDS Rate z-score", 
     ylab="Lagged AIDS Rate z-score", main="Moran Scatterplot AIDS Rate",
     xlim=c(-5,5), ylim=c(-3,3), pch = 20)
abline(h=0, v=0)
abline(lm(lagScrate ~ Scrate), lty=2, lwd=2, col="blue")

# Now, making the LISA map

ny$hh <- (Scrate>= 0 & lagScrate>= 0) & (locm[,5]<= 0.05)
ny$ll <- (Scrate<= 0 & lagScrate<= 0) & (locm[,5]<= 0.05)
ny$hl <- (Scrate>= 0 & lagScrate<= 0) & (locm[,5]<= 0.05)
ny$lh <- (Scrate<= 0 & lagScrate>= 0) & (locm[,5]<= 0.05)
ny$ns <- locm[,5]> 0.05

# Create a single categorial variable summing up the
# five logical variables

ny$var <- 0
ny$var <- replace(ny$var, (ny$hh==TRUE), 1)
ny$var <- replace(ny$var, (ny$ll==TRUE), 2)
ny$var <- replace(ny$var, (ny$hl==TRUE), 3)
ny$var <- replace(ny$var, (ny$lh==TRUE), 4)
ny$var <- replace(ny$var, (ny$ns==TRUE), 5)

# Set the breaks for the thematic map classes

breaks <-seq(1,5,1)

# Set the corresponding labels for the thematic map classes

labels <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not Signif.")

# Determine which map class each observation falls into based on
# the value of ny$var

np <- findInterval(ny$var, breaks, all.inside=FALSE)

# Assign colors to each map class

colors <- c("red", "blue", "pink", "skyblue2", "gray90")

# Plot the SpatialPolygonsDataFrame using specified breaks
# and color scheme

plot(ny, col=colors[np], border = "white", lwd = .25)
mtext("LISA Map AIDS Rate; weights: Q1", cex=1.2, side=3, line=1)
legend("left", legend=labels, fill=colors, bty="n")

############
### Running a basic OLS regression
############

lrate <- log1p(ny$Rate1000) #logging our DV

reg1 <- lm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH, data = ny)
summary(reg1)

# Testing for autocorrelation in our residuals:
ny$resid <- reg1$fitted.values

moran.test(ny$resid, nyQ1.gal, alternative="two.sided", zero.policy=T)

############
### Where does the model fit poorly?
############

locm_resid <- localmoran(ny$resid, nyQ1.gal, zero.policy=T, 
  alternative="two.sided")
Scrate_resid <- (ny$resid - mean(ny$resid))/sd(ny$resid)
lagScrate_resid <- lag.listw(nyQ1.gal, Scrate, zero.policy=T)


ny$hh_resid <- (Scrate_resid>= 0 & lagScrate_resid>= 0) & 
  (locm_resid[,5]<= 0.05)
ny$ll_resid <- (Scrate_resid<= 0 & lagScrate_resid<= 0) & 
  (locm_resid[,5]<= 0.05)
ny$hl_resid <- (Scrate_resid>= 0 & lagScrate_resid<= 0) & 
  (locm_resid[,5]<= 0.05)
ny$lh_resid <- (Scrate_resid<= 0 & lagScrate_resid>= 0) & 
  (locm_resid[,5]<= 0.05)
ny$ns_resid <- locm[,5]> 0.05

# Create a single categorial variable summing up the
# five logical variables

ny$var_resid <- 0
ny$var_resid <- replace(ny$var_resid, (ny$hh_resid==TRUE), 1)
ny$var_resid <- replace(ny$var_resid, (ny$ll_resid==TRUE), 2)
ny$var_resid <- replace(ny$var_resid, (ny$hl_resid==TRUE), 3)
ny$var_resid <- replace(ny$var_resid, (ny$lh_resid==TRUE), 4)
ny$var_resid <- replace(ny$var_resid, (ny$n_resids==TRUE), 5)

# Set the breaks for the thematic map classes

breaks <-seq(1,5,1)

# Set the corresponding labels for the thematic map classes

labels <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not Signif.")

# Determine which map class each observation falls into based on
# the value of ny$var

np_resid <- findInterval(ny$var_resid, breaks, all.inside=FALSE)

# Assign colors to each map class

colors <- c("red", "blue", "pink", "skyblue2", "gray90")

# Plot the SpatialPolygonsDataFrame using specified breaks
# and color scheme

plot(ny, col=colors[np_resid], border="white", lwd = .25)
mtext("LISA Map AIDS Rate; weights: Q1", cex=1.2, side=3, line=1)
legend("left", legend=labels, fill=colors, bty="n")

##############
### Spatial Regressions
##############

# Precalculating the Eigenvalues, which speeds things up
W.eig <- eigenw(nyQ1.gal, quiet=NULL)

# The spatial lag model
ny.lag.eig <- lagsarlm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH,
  data=ny, nyQ1.gal, method="eigen", quiet=TRUE,
  control=(pre_eig=W.eig)) # How to use pre-calculated eigenvalues
summary(ny.lag.eig)


# The spatial error model
ny.err.eig <- errorsarlm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH,
   data=ny, nyQ1.gal, method="eigen", quiet=TRUE, control=(pre_eig=W.eig))
summary(ny.err.eig)

##############
### Did these resolve spatial autocorrelation?
##############

# Yes!
moran.test(ny.err.eig$resid, nyQ1.gal, alternative="two.sided", zero.policy=T)


##############
### Running a CAR model
##############
library(CARBayes)
# CARBayes requires a different weight matrix structure
w_mat <- nb2mat(nygal, zero.policy=T)

# And the matrix has to be symmetric
w_mat <- ifelse(w_mat != 0, 1, 0)

set.seed(614) # for reproducible sampling

# uninformative (default) priors
car1 <- S.CARleroux(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH, 
  "gaussian", data=ny@data, W = w_mat,  burnin=20000, n.sample=100000, thin=10, 
  prior.nu2=NULL, prior.tau2=NULL, fix.rho=FALSE, rho=NULL, 
  verbose=FALSE)
# can't just type summary()
car1$summary.results 


# weakly informative priors -- note the larger number of effective samples
car2 <- S.CARleroux(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH, 
  "gaussian", data=ny@data, W = w_mat,  burnin=20000, n.sample=100000, thin=10, 
  prior.nu2=c(10, 1), prior.tau2=c(10, 1), fix.rho=FALSE, rho=NULL, 
  verbose=FALSE)
# can't just type summary()
car2$summary.results 

##############
### Running a GWR model
##############

library(GWmodel) # v. 2.0-2
# Creating the distance matrix, which is different from the other two matrices
dist_mat <- gw.dist(dp.locat = coordinates(ny), rp.locat =coordinates(ny), 
  focus = 0, p = 2)

# Calculating the kernel bandwidth
bw1 <- bw.gwr(lrate ~ Gini, approach="aic",adaptive=TRUE, data=ny, 
  kernel = "gaussian", dMat=dist_mat)

# Running the GWR
gw.model1 <- gwr.basic(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH, 
  bw = bw1, data = ny, kernel = "gaussian", adaptive = TRUE, p = 2, 
  longlat = FALSE)

# regression table output
gw.model1

# But plotting this would make more sense
# Let's plot Gini, since that changed the most in the SAR/SEM vs. OLS models
# SDF is a spatial dataframe which links data to spatial location
# Gini is the variable
ny$gw_out <- gw.model1$SDF$Gini

# To get cutpoints
quantile(ny$gw_out)

# Color scheme to use
shades_to_use <- shading(breaks = c(4.57, 6.33, 7.49), 
  cols = brewer.pal(4, "Reds"))

#plotting the map
choropleth(sp = ny, v = ny$gw_out, border = "white", lwd = .25, 
  shading = shades_to_use)
#making the legend
legend_colors <- shades_to_use$cols
legend("left", legend = c("0-25", "25-50", "50-75", "75-100"), 
  fill = legend_colors,bg="transparent", box.col = "transparent", 
  border = "white", title = "GWR Beta values for Gini, \nby quantile")


###############
### Other Useful packages
###############
# Shapefile maps of the world, 1946-2012
install.packages("cshapes") 

# Mapping packages
# Displaying maps, depends on the next two packages
install.packages("maps") 
install.packages("mapproj")
install.packages("mapdata")

# Spatial Stats packages
install.packages("fields") # Curve, surface, and function fitting (Splines and Kriging)
install.packages("geoR") # Useful for kriging
install.packages("rgeos")# Useful for kriging

install.packages("CARBayesST") # Spatio-temporal CAR models
install.packages("McSpatial") # Frequentist approaches to Logit/GMM models
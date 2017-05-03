##' Calculation and visualization of gemplots (3-dimensional extension
##' of boxplot and bagplot).
##' 
##' Gemplots are the 3-dimensional extension of the 1-dimensional
##' boxplots and the 2-dimensional bagplot. The package calculates a
##' 3-dimensional median (\code{\link{depmed}}), an inner bag
##' (\code{\link{bag}}) that contains 50 percent of the data and an
##' outer bag (\code{\link{loop}}). Data points outside the outer bag
##' are flagged as outliers. The gemplot can be visualized in a
##' 3-dimensional interactice device (\code{\link{gem}}). The
##' calculations are based on the concept of halfspace location
##' depths. The halfspace location depths need to be calculated in a
##' user specified grid prior to the calculation of the gemplots
##' components (\code{\link{gridfun}}, \code{\link{hldepth}}).
##' Outlier detection for more than 3 dimensions is also implemented
##' in the above workflow. An aditional parameter k needs then to be
##' specified in the \code{\link{gridfun}} function.
##' 
##' @references
##' \itemize{
##' \item Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The
##'   bagplot: a bivariate boxplot. The American Statistician, 53(4),
##'   382-387.
##' \item Kruppa, J. and Jung, K. (2017) Automated multigroup
##'   outlier identification in molecular high-throughput data using
##'   bagplots and gemplots. BMC Bioinformatics, 18, 232 
##' }
##' @docType package
##' @name gemPlot
##' @examples
##' ### Two 3-dimensional example data sets D1 and D2
##' n <- 200
##' x1 <- rnorm(n, 0, 1)
##' y1 <- rnorm(n, 0, 1)
##' z1 <- rnorm(n, 0, 1)
##' D1 <- data.frame(cbind(x1, y1, z1))
##' x2 <- rnorm(n, 1, 1)
##' y2 <- rnorm(n, 1, 1)
##' z2 <- rnorm(n, 1, 1)
##' D2 <- data.frame(cbind(x2, y2, z2))
##' colnames(D1) <- c("x", "y", "z")
##' colnames(D2) <- c("x", "y", "z")
##' 
##' ### Placing outliers in D1 and D2
##' D1[17,] = c(4, 5, 6)
##' D2[99,] = -c(3, 4, 5)
##' 
##' ### Grid size and graphic parameters
##' grid.size <- 20
##' red <- rgb(200, 100, 100, alpha = 100, maxColorValue = 255)
##' blue <- rgb(100, 100, 200, alpha = 100, maxColorValue = 255)
##' yel <- rgb(255, 255, 102, alpha = 100, maxColorValue = 255)
##' white <- rgb(255, 255, 255, alpha = 100, maxColorValue = 255)
##' require(rgl)
##' material3d(color=c(red, blue, yel, white), alpha=c(0.5, 0.5, 0.5, 0.5), smooth=FALSE, specular="black")
##' 
##' ### Calucation and visualization of gemplot for D1
##' G <- gridfun(D1, grid.size=20)
##' G$H <- hldepth(D1, G, verbose=TRUE)
##' dm <- depmed(G)
##' B <- bag(D1, G)
##' L <- loop(D1, B, dm)
##' rgl.open()
##' points3d(D1[L$outliers==0,1], D1[L$outliers==0,2], D1[L$outliers==0,3], col=red)
##' text3d(D1[L$outliers==1,1], D1[L$outliers==1,2], D1[L$outliers==1,3], as.character(which(L$outliers==1)), col=yel)
##' spheres3d(dm[1], dm[2], dm[3], col=yel, radius=0.1)
##' gem(B$coords, B$hull, red)
##' gem(L$coords.loop, L$hull.loop, red)
##' axes3d(col="white")
##' 
##' ### Calucation and visualization of gemplot for D2
##' G <- gridfun(D2, grid.size=20)
##' G$H <- hldepth(D2, G, verbose=TRUE)
##' dm <- depmed(G)
##' B <- bag(D2, G)
##' L <- loop(D2, B, dm)
##' points3d(D2[L$outliers==0,1], D2[L$outliers==0,2], D2[L$outliers==0,3], col=red)
##' text3d(D2[L$outliers==1,1], D2[L$outliers==1,2], D2[L$outliers==1,3], as.character(which(L$outliers==1)), col=yel)
##' spheres3d(dm[1], dm[2], dm[3], col=yel, radius=0.1)
##' gem(B$coords, B$hull, blue)
##' gem(L$coords.loop, L$hull.loop, blue)
NULL

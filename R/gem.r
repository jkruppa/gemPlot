##' Specifies a k-dimensional array as grid for the calculation of the
##' halfspace location depths.
##'
##' D must have at least three columns. If D has three columns,
##' automatically a 3-dimensional grid is generated. If D has more
##' than three columns, k must be specified.
##' @title Specifies grid for the calculation of the halfspace location depths
##' @param D Data set with rows representing the individuals and
##'     columns representing the features. In the case of three
##'     dimensions, the colnames of D must be c("x", "y", "z").
##' @param grid.size Number of grid points in each dimension.
##' @param k Number of dimensions of the grid. Needs only be specified
##'     if D has more than columns.
##' @return
##' A list containing the following elements:
##' \describe{
##' \item{\emph{H}}{The k-dimensional array.}
##' }
##' In the case of a 3-dimensional array, additional elements are:
##' \describe{
##' \item{\emph{grid.x, grid.y, grid.z}}{The coordinates of the grid points at each dimension.}
##' }
##' In the case that the array has more than three dimensions, additional elements are:
##' \describe{
##' \item{\emph{grid.k}}{A matrix with the coordinates of the grid. Row represents dimensions and columns represent grid points.}
##' }
##' @author Jochen Kruppa, Klaus Jung
##' @export
gridfun <- function(D, grid.size, k) {
    grid.x <- seq(min(D$x)[1], max(D$x)[1], length.out = grid.size)
    grid.y <- seq(min(D$y)[1], max(D$y)[1], length.out = grid.size)
    grid.z <- seq(min(D$z)[1], max(D$z)[1], length.out = grid.size)
    H <- array(NA, c(grid.size, grid.size, grid.size))
    return(list(grid.x = grid.x, grid.y = grid.y, grid.z = grid.z, H = H))
}


##' Calculates the halfspace location depth for each point in a given grid.
##'
##' Calculation of the halfspace location depth at each grid point is
##' mandatory before calculating the depth median
##' (\code{\link{depmed}}), the bag (\code{\link{bag}}) and the loop
##' (\code{\link{loop}}). Ideally, the output is assigned to the array
##' H produced by \code{\link{gridfun}}.
##' @title Calculates the halfspace location depth
##' @param D Data set with rows representing the individuals and
##'     columns representing the features. In the case of three
##'     dimensions, the colnames of D must be c("x", "y", "z").
##' @param G List containing the grid information produced by
##'     \code{\link{gridfun}}.
##' @param verbose Logical. Indicates whether progress information is
##'     printed during calculation.
##' @return
##' \describe{
##' \item{\emph{H}}{An array of the same dimension as the array in argument G. The elements contain the halfspace location depth at the related grid location.}
##' }
##' 
##' @references Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The
##'     bagplot: a bivariate boxplot. The American Statistician,
##'     53(4), 382-387.
##' @author Jochen Kruppa, Klaus Jung
##' @export
##' @examples
##' ### A 3-dimensional example data set D1
##' n <- 200
##' x1 <- rnorm(n, 0, 1)
##' y1 <- rnorm(n, 0, 1)
##' z1 <- rnorm(n, 0, 1)
##' D1 <- data.frame(cbind(x1, y1, z1))
##' colnames(D1) <- c("x", "y", "z")
##' 
##' ### Specification of the grid and calculation of the halfspace location depth at each grid location.
##' G <- gridfun(D1, grid.size=20)
##' G$H <- hldepth(D1, G, verbose=TRUE)
hldepth <- function(D, G, verbose=TRUE){
    require(depth)
    n <- dim(D)[1]
    grid.size <- length(G$grid.x)
    perc <- 10
    H <- G$H
    for (i in 1:grid.size) {
        for (j in 1:grid.size) {
            for (k in 1:grid.size) {
                u <- c(G$grid.x[i], G$grid.y[j], G$grid.z[k])
                H[i,j,k] <- n * depth(u, D)
            }
        }
        if (100 * i/grid.size >= perc && verbose == TRUE) {
            message(paste("Calculation of halfspace location depths: ", round(100 * i/grid.size, 0), " % of grid points done", sep=""))
            perc <- perc + 10
        }
    }
    return(H)
}

##' Calculates the depth median.
##'
##' Calculates the depth median in a specified grid array with given
##' halfspace location depth at each grid location.
##' @title Calculates the depth median.
##' @param G List containing the grid information produced by
##'     \code{\link{gridfun}} and the halfspace location depths
##'     produced by \code{\link{hldepth}}.
##' @return An vector with a length equal to the number of dimension
##'     of the array in G, containing the coordinates of the depth
##'     median.
##' @author Jochen Kruppa, Klaus Jung
##' @references Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The
##'     bagplot: a bivariate boxplot. The American Statistician,
##'     53(4), 382-387.
##' @export
##' @examples
##' ### A 3-dimensional example data set D1
##' n <- 200
##' x1 <- rnorm(n, 0, 1)
##' y1 <- rnorm(n, 0, 1)
##' z1 <- rnorm(n, 0, 1)
##' D1 <- data.frame(cbind(x1, y1, z1))
##' colnames(D1) <- c("x", "y", "z")
##' 
##' ### Specification of the grid and calculation of the halfspace location depth at each grid location.
##' G <- gridfun(D1, grid.size=20)
##' G$H <- hldepth(D1, G, verbose=TRUE)
##' dm <- depmed(G) ### Calculation of the depth median
depmed <- function(G) {
    grid.size <- length(G$grid.x)
    maxi <- max(G$H)[1]
    H2 <- (G$H == maxi)
    i2 <- rep(0, grid.size)
    j2 <- rep(0, grid.size)
    k2 <- rep(0, grid.size)
    for (i in 1:grid.size) {
        for (j in 1:grid.size) {
            for (k in 1:grid.size) {
                if (H2[i,j,k] == TRUE) {
                    i2[i] <- i
                    j2[j] <- j
                    k2[k] <- k
                }
            }
        }
    }
    med.i <- median(G$grid.x[which(i2 != 0)])
    med.j <- median(G$grid.y[which(j2 != 0)])
    med.k <- median(G$grid.z[which(k2 != 0)])
    return(c(med.i, med.j, med.k))
}

##' Calculates the bag of a gemplot (i.e. the inner gemstone).
##'
##' Determines those grid points that belong to the bag, i.e. a convex
##' hull that contains 50 percent of the data. In the case of a
##' 3-dimensional data set, the bag can be visualized by an inner
##' gemstone that can be accompanied by an outer gemstone (\code{\link{loop}}).
##'
##' @title  Calculates the bag
##' @param D Data set with rows representing the individuals and
##'     columns representing the features. In the case of three
##'     dimensions, the colnames of D must be c("x", "y", "z").
##' @param G List containing the grid information produced by
##'     \code{\link{gridfun}} and the halfspace location depths calculated by
##'     \code{\link{hldepth}}.
##' 
##' @return A list containg the following elements:
##' \describe{
##' \item{\emph{coords}}{Coordinates of the grid points that belong to
##'     the bag. Each row represents a grid point and each column
##'     represents one dimension.}
##' \item{\emph{hull}}{A data matrix that
##'     contains the indices of the margin grid points of the bag that
##'     cover the convex hull by triangles. Each row represents one
##'     triangle. The indices correspond to the rows of coords.}
##' }
##' 
##' @references Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The
##'     bagplot: a bivariate boxplot. The American Statistician,
##'     53(4), 382-387.
##' @author Jochen Kruppa, Klaus Jung
##' @export
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
bag <- function(D, G) {
    require(geometry)
    n <- dim(D)[1]
    D.k <- rep(NA, n)
    for (i in 1:n) {
        I <- matrix(NA, 8, 3)
        I[1,] <- c(G$grid.x[max(which(G$grid.x<=D[i,1]))],
                   G$grid.y[max(which(G$grid.y<=D[i,2]))],
                   G$grid.z[max(which(G$grid.z<=D[i,3]))])
        I[2,] <- c(G$grid.x[max(which(G$grid.x<=D[i,1]))],
                   G$grid.y[max(which(G$grid.y<=D[i,2]))],
                   G$grid.z[min(which(G$grid.z>=D[i,3]))])
        I[3,] <- c(G$grid.x[max(which(G$grid.x<=D[i,1]))],
                   G$grid.y[min(which(G$grid.y>=D[i,2]))],
                   G$grid.z[max(which(G$grid.z<=D[i,3]))])
        I[4,] <- c(G$grid.x[max(which(G$grid.x<=D[i,1]))],
                   G$grid.y[min(which(G$grid.y>=D[i,2]))],
                   G$grid.z[min(which(G$grid.z>=D[i,3]))])
        I[5,] <- c(G$grid.x[min(which(G$grid.x>=D[i,1]))],
                   G$grid.y[max(which(G$grid.y<=D[i,2]))],
                   G$grid.z[max(which(G$grid.z<=D[i,3]))])
        I[6,] <- c(G$grid.x[min(which(G$grid.x>=D[i,1]))],
                   G$grid.y[max(which(G$grid.y<=D[i,2]))],
                   G$grid.z[min(which(G$grid.z>=D[i,3]))])
        I[7,] <- c(G$grid.x[min(which(G$grid.x>=D[i,1]))],
                   G$grid.y[min(which(G$grid.y>=D[i,2]))],
                   G$grid.z[max(which(G$grid.z<=D[i,3]))])
        I[8,] <- c(G$grid.x[min(which(G$grid.x>=D[i,1]))],
                   G$grid.y[min(which(G$grid.y>=D[i,2]))],
                   G$grid.z[min(which(G$grid.z>=D[i,3]))])
        I <- cbind(I, NA)
        for (t in 1:8) {
            index1 <- match(I[t,1], G$grid.x)
            index2 <- match(I[t,2], G$grid.y)
            index3 <- match(I[t,3], G$grid.z)
            I[,4] <- G$H[index1, index2, index3]
        }
        D.k[i] <- min(I[,4])
    }
    H2 <- (G$H>=max(which(cumsum(table(D.k))<=(n/2))))
    BAG <- matrix(NA, 0, 3)
    for (i in 1:grid.size) {
        for (j in 1:grid.size) {
            for (k in 1:grid.size) {
                if (H2[i,j,k]==TRUE){
                    BAG <- rbind(BAG, c(G$grid.x[i],
                                        G$grid.y[j],
                                        G$grid.z[k]))
                }
            }
        }
    }
    convH <- convhulln(BAG)
    return(list(coords=BAG, hull=convH))
}

##' Calculates the fence and the loop of a gemplot (i.e. the outer gemstone).
##'
##' The fence inflates the the bag relative to the depth median by the
##' factor inflation. Data points outside the bag and inside the fence
##' the loop or outer gemstone are flagged as outliers. Data points
##' outside the fence are marked as outliers. In the case of a
##' 3-dimensional data set, the loop can be visualized by an outer
##' gemstone around the inner gemstone or bag.
##' @title Calculates the fence and the loop
##' @param D Data set with rows representing the individuals and
##'     columns representing the features. In the case of three
##'     dimensions, the colnames of D must be c("x", "y", "z").
##' @param B List containing the information about the coordinates of
##'     the bag and the convex hull that forms the bag (determined by
##'     \code{\link{bag}}).
##' @param inflation A numeric value > 0 that specifies the inflation
##'     factor of the bag relative to the median (default = 3).
##' @param dm The coordinates of the depth median as produced by
##'     \code{\link{depmed}}.
##' @return A list containing the following elements:
##' \describe{
##' \item{\emph{coords.loop}}{Coordinates of the data points that are inside the convex hull around the loop.}
##' \item{\emph{hull.loop}}{A data matrix that contains the indices of the margin data points of the loop that cover the convex hull by triangles. Each row represnts one triangle. The indices correspond to the rows of coords.loop.}
##' \item{\emph{coords.fence}}{Coordinates of the grid points that are inside the fence but outside the bag.}
##' \item{\emph{hull.fence}}{A data matrix that contains the indices of the margin grid points of the fence that cover the convex hull around the fence by triangles. Each row represnts one triangle. The indices correspond to the rows of coords.fence.}
##' \item{\emph{outliers}}{A vector of length equal to the sample size. Data points that are inside the fence are labelled by 0 and values outside the fence (i.e. outliers) are labelled by 1.}
##' }
##' 
##' @references Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The
##'     bagplot: a bivariate boxplot. The American Statistician,
##'     53(4), 382-387.
##' @author Jochen Kruppa, Klaus Jung
##' @export
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
loop <- function(D, B, dm, inflation = 3) {
    n <- dim(D)[1]
    index.F <- sort(intersect(as.vector(B$hull), as.vector(B$hull)))
    FENCE <- B$coords[index.F,]
    MED.MAT <- t(matrix(c(dm[1], dm[2], dm[3]), 3, dim(FENCE)[1]))
    FENCE <- MED.MAT + 3 * (FENCE - MED.MAT)
    colnames(FENCE) <- colnames(D)
                convH <- convhulln(FENCE)
    outliers <- rep(0, n)
    for (i in 1:n) {
        Z <- rbind(FENCE, D[i,])
        convH.Z <- convhulln(Z)
        if (!is.na(match(dim(FENCE)[1]+1, convH.Z))){
            outliers[i] <- 1
        }
    }
    LOOP <- D[which(outliers==0),]
    convH2 <- convhulln(LOOP)
    return(list(coords.loop=LOOP, hull.loop=convH2, coords.fence=FENCE, hull.fence=convH, outliers=outliers))
}

##' Plots a gemstone to an interactive graphics device.
##'
##' Only applicable to 3-dimensional data sets. Transparent colors are
##' recommended for outer gemstone of the gemplot. Further graphical
##' parameters can be set using material3d() of the rgl-package.
##' 
##' @title Plots a gemstone to an interactive graphics device
##' @param coords Matrix with coordinates of the grid or of data
##'     points that belong to the gemstone, calculated by either
##'     \code{\link{bag}} or \code{\link{loop}}. Each row represents a
##'     grid point and each column represents one dimension.
##' @param hull Matrix with indices of triangles that cover a convex
##'     hull arround the gemstone. Each row represents one triangle
##'     and the indices refer to the rows of coords.
##' @param clr Specifies the color of the gemstone.
##' @return NULL
##' @author Jochen Kruppa, Klaus Jung
##' @export
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
##' spheres3d(dm[1], dm[2], dm[3], col=white, radius=0.1)
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
##' points3d(D2[L$outliers==0,1], D2[L$outliers==0,2], D2[L$outliers==0,3], col=blue)
##' text3d(D2[L$outliers==1,1], D2[L$outliers==1,2], D2[L$outliers==1,3], as.character(which(L$outliers==1)), col=yel)
##' spheres3d(dm[1], dm[2], dm[3], col=white, radius=0.1)
##' gem(B$coords, B$hull, blue)
##' gem(L$coords.loop, L$hull.loop, blue)
gem <- function(coords, hull, clr) {
    require(geometry)
    for (i in 1:dim(hull)[1]) {
        x <- coords[hull[i,],1]
        y <- coords[hull[i,],2]
        z <- coords[hull[i,],3]
        triangles3d(x, y, z, col=clr)
    }
}


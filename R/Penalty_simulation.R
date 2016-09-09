#' @title Spatial simulation
#' 
#' @description Simulate a geometrical figure on a 2D spatial field
#' 
#' @param n.obs the number of repetitions.
#' @param xmax the length of the spatial field.
#' @param ymax the width of the spatial field
#' @param radius the maximum distance to the center for which pixels will set to 1. The others will be set to 0.
#' @param center the position relative to which the distance will be computed.
#' @param coords.centered should the coordinates be centered around 0,0? Else they will start at 0,0 and take positve values. May affect the computation of the distance.
#' @param distance the type of distance to be used.
#' 
#' @details The available distances are those of the \code{dist} function from the \emph{stats} package.
#' 
#' @return a list containing:
#' \itemize{
#' \item{"coords"} the coordinates of each point on the field
#' \item{"center"} the center used to compute the distance
#' \item{"distCenter"} the distance of each point to the center
#' \item{"X"} the distance thresholded by the radius
#' }
#' 
#' @examples 
#' image(simForm(100, 10, 10, 2)$X)
#' image(simForm(100, 10, 10, 2, distance = "maximum")$X)
#' 
#' @export
simForm <- function(n.obs, xmax, ymax, radius, center = NULL, coords.centered = TRUE,
                    distance = "euclidean"){
  
  if(is.null(center)){
    if(coords.centered == TRUE){
      center <- c(0,0)
    }else{
      center <- c(xmax/2,ymax/2)
    }
  }
  
  coords <- scale(expand.grid(1:xmax, 1:ymax), center = coords.centered, scale = FALSE)
  n.coord <- nrow(coords) 
  
  distCenter <- apply(coords, 1, function(x){dist( rbind(x,center), method = distance)})
  beta <- distCenter<radius
  
  return(list(coords = coords,
              center = center,
              distCenter = matrix(distCenter, nrow = xmax, ncol = ymax),
              X = matrix(beta, nrow = xmax, ncol = ymax)
  ))
}
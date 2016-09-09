#' @export
`plot.plvmfit` <- function(x, ...) {
  plot(x$regularizationPath, ...)
}  

#' @export
`plot.regPath` <- function(x, lambda = getLambda(x), type = NULL, 
                           coefficient = "penalized", row = NULL,
                           add.line = TRUE, line.size = 2,
                           add.point = TRUE, point.shape = 4, point.size = 2,
                           add.best = TRUE, color.selected = TRUE) {
  
  if(is.null(type)){
    if(!is.null(x$performance)){
      type <- "criterion"
    }else{
      type <- "path"
    }
  }
  
  if(type == "path"){
    path <- getPath(x, coefficient = coefficient, lambda = lambda, row = row)
    df.Path <- data.table::melt(path, 
                                measure=names(path)[-1], 
                                value.name = "value", variable.name = "coefficient")
    
    ggPath <- ggplot(df.Path, aes_string(y = "value", 
                                         x = lambda, 
                                         group = "coefficient", 
                                         col = "coefficient")
    )
    if(add.line){ggPath <- ggPath + geom_line(size = line.size)}
    if(add.point){ggPath <- ggPath + geom_point(size = point.size, shape = point.shape)}
    
    if(!is.null(x$performance)){
      if(color.selected){
        df.Path$selected <- df.Path$coefficient %in% names(x$optimum$coef)
        names.selected <- unique(df.Path$coefficient[df.Path$selected])
        names.Nselected <- unique(df.Path$coefficient[!df.Path$selected])
        n.selected <- length(names.selected)
        n.Nselected <- length(names.Nselected)
        color.selected <- rgb(green = seq(0.3,0.7,length.out = n.selected), 0, 0)
        color.Nselected <- rgb(red = seq(0.3,0.7,length.out = n.Nselected), 0, 0)
        
        color.order <- as.character(unique(df.Path$coefficient))
        color.order[color.order %in% names.selected] <- color.selected
        color.order[color.order %in% names.Nselected] <- color.Nselected
        ggPath <- ggPath + scale_color_manual(values = color.order)
      }
      
      if(add.best){
        ggPath <- ggPath + geom_vline(size = line.size/2, 
                                      xintercept = x$optimum[[lambda]], 
                                      linetype = 2, color = "blue")
      }
      
    }
    
    return(ggPath)
    
  }else if(type == "criterion"){
    performance <- getPerformance(x)[order(getPerformance(x)[[lambda]]),]
    
    df <- data.frame(lambda = performance[[lambda]],
                     criterion = performance$value,
                     optimum = performance$optimum)
    names(df)[1] <- lambda
    names(df)[2] <- x$optimum$criterion
    
    ggPerf <- ggplot(df, aes_string(x = lambda, y = names(df)[2]))
    if(add.line){ggPerf <- ggPerf + geom_line(size = line.size)}
    if(add.point){ggPerf <- ggPerf + geom_point(size = point.size, aes_string(color = "optimum"))}
    
    return(ggPerf)
  }else{
    stop("type must be \"path\" or \"criterion\" \n")
  }
}

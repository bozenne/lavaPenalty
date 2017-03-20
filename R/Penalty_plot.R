#' @export
`plot.plvmfit` <- function(x, coef, ...) {

    if(is.null(x$regularizationPath)){
        stop("not implemented yet!")
    }else{

        if(!missing(coef)){
            if(identical(coef, "penalized")){
                table.penalty <- penalty(x, nuclear = FALSE, type = NULL)
                coef <- table.penalty[penalty == "lasso",link]            
            }else if(identical(coef, "Npenalized")){
                table.penalty <- penalty(x, nuclear = FALSE, type = NULL)
                coef <- setdiff(coef(x),table.penalty[penalty == "lasso",link])
            }
        }
        
    }
    plot(x$regularizationPath, coef = coef, ...)
}  

#' @export
`plot.regPath` <- function(x, coef, lambda = "lambda1", 
                           type = NULL, row = NULL, 
                           xlim = NULL, ylim = NULL,
                           add.line = TRUE, line.size = 2,
                           add.point = TRUE, point.shape = 4, point.size = 2,
                           add.best = TRUE, color.selected = TRUE, plot = TRUE) {
    if(is.null(type)){
        if(!is.null(x$performance)){
            type <- "criterion"
        }else{
            type <- "path"
        }
    }else if(type %in% c("criterion","path") == FALSE){
        stop("wrong specification of argument \'type\' \n",
             "valid type: \"criterion\" \"path\" \n")
    }

    if(type == "path"){
        if(length(lambda)!=1){
            stop("Argument \'lambda\' must have length one \n")
        }
        dt.path <- getPath(x, coef = coef, lambda = lambda, row = row,
                           keep.index = FALSE, keep.indexChange = FALSE, keep.optimum = FALSE)
        dtL.path <- data.table::melt(dt.path, 
                                     measure=setdiff(names(dt.path),lambda), 
                                     value.name = "value", variable.name = "coefficient")

        ggPath <- ggplot() + coord_cartesian(xlim = xlim, ylim = ylim)
        if(add.line){
            ggPath <- ggPath + geom_line(aes_string(y = "value", x = lambda, group = "coefficient", col = "coefficient"),
                                         size = line.size, data = dtL.path)}
        if(add.point){
            ggPath <- ggPath + geom_point(aes_string(y = "value", x = lambda, group = "coefficient", col = "coefficient"),
                                          size = point.size, shape = point.shape, data = dtL.path)
        }

        if(!is.null(x$criterion)){
            if(color.selected){
                name.coef <- setdiff(names(x$path),
                                     c("index","lambda1.abs","lambda1","lambda2.abs","lambda2","indexChange",x$criterion,"cv","optimum")
                                     )
                
                coefOptimum <- x$path[optimum==TRUE,name.coef,with = FALSE]
                names.selected <- names(coefOptimum)[abs(coefOptimum)>1e-12]
                names.Nselected <- names(coefOptimum)[abs(coefOptimum)<1e-12]
                n.selected <- length(names.selected)
                n.Nselected <- length(names.Nselected)
                color.selected <- rgb(green = seq(0.3,0.7,length.out = n.selected), 0, 0)
                color.Nselected <- rgb(red = seq(0.3,0.7,length.out = n.Nselected), 0, 0)
        
                color.order <- as.character(unique(dtL.path$coefficient))
                color.order[color.order %in% names.selected] <- color.selected
                color.order[color.order %in% names.Nselected] <- color.Nselected
                ggPath <- ggPath + scale_color_manual(values = color.order)
            }
      
            if(add.best){                
              ggPath <- ggPath + geom_vline(size = line.size/2, 
                                            xintercept = x$path[optimum==TRUE][[lambda]], 
                                            linetype = 2, color = "blue")
          }
      
        }

        ls.res <- list(plot = ggPath,
                       data = dtL.path)        
    
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
    
    ls.res <- list(plot = ggPerf,
                   data = df)
           
  }

    if(plot){
        print(ls.res$plot)
    }
    return(invisible(ls.res))
}

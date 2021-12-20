#' This is the function to run pathway analysis
#' @export

run_path_analysis<- function(lf.c, gsets ){
    fc.kegg.p <- gage(lf.c, gsets , ref = NULL, samp = NULL)
    
 
    #return(path.ids.2)
    return(fc.kegg.p)
}

#' This is the function to run pathway analysis
#' @import gage
#'
#' @export
run_path_analysis<- function(lf.c, gsets ){
    fc.kegg.p <- gage::gage(lf.c, gsets , ref = NULL, samp = NULL)
    return(fc.kegg.p)
}

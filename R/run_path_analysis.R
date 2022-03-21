#' This is the function to run pathway analysis
#' @import gage
#'
#' @export
run_path_analysis<- function(lf.c, gsets, ref = ref, samp = samp, compareval ){
    fc.kegg.p <- gage::gage(exprs = lf.c, gsets=gsets , ref = NULL, samp = NULL,compare= compareval)
    return(fc.kegg.p)
}

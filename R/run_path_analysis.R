#' This is the function to run pathway analysis
#' @export

run_path_analysis<- function(lf.c, gsets ){
    fc.kegg.p <- gage(lf.c, gsets , ref = NULL, samp = NULL)
    sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
        + !is.na(fc.kegg.p$greater[, "q.val"])
    path.ids.2 <- rownames(fc.kegg.p$greater)[sel]
    if(length(fc.kegg.p) > 2){
    sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
        + !is.na(fc.kegg.p$less[,"q.val"])
    path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
    path.ids.2 <- substr(c(path.ids.2, path.ids.l), 1, 8)
    }
    return(path.ids.2)
}

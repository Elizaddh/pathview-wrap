#' this is a function to run pathview.2 function
#' @import pathview


#' this is a function to run pathview.2 function
#' @export


pathview.2 <- function( run ,diff.tool, gene.data= cnts, ref, samp, gsets, pathway.id,  both.dirs = list(gene = T, cpd = T), plot.gene.data, outname)
{
  library(pathview)
  if(is.null(pathway.id)==FALSE){
    pathview( gene.data = gene.data,  pathway.id = pathway.id , out.suffix=outname)
    
  }
  else {
    if (run =="complete")
    {
      logfoldchange <- rundifftool(diff.tool, gene.data, ref, samp, outname)
      print("diff tool run successful")
      fc.kegg.p <- run_path_analysis(logfoldchange, gsets)#, gene.data = gene.data, ref, samp, plot.gene.data = T )
      print("gage run successful")
      print("now pathview")
      
      path.ids.2<- rownames(fc.kegg.p$greater)[fc.kegg.p$greater[, "q.val"] < 0.1 &
                                                 + !is.na(fc.kegg.p$greater[, "q.val"])]
      
      if(length(fc.kegg.p) > 2){
        print(length(fc.kegg.p)) 
        path.ids.l <- rownames(fc.kegg.p$less)[fc.kegg.p$less[, "q.val"] < 0.1 &
                                                 + !is.na(fc.kegg.p$less[,"q.val"])]
        path.ids.2 <- substr(c(path.ids.2[1:3], path.ids.l[1:3]), 1, 8)
      }
      #visualize pathway  
      pv.out.list <- sapply(na.omit(path.ids.2[1:6]), function(pid) pathview( gene.data = logfoldchange, 
                                                                              pathway.id = pid, out.suffix=diff.tool))
    }
    else{
      fc.kegg.p <- run_path_analysis(gene.data, gsets)#, gene.data = gene.data, ref, samp, plot.gene.data = T  )
      print("now pathview")
      path.ids.2<- rownames(fc.kegg.p$greater)[fc.kegg.p$greater[, "q.val"] < 0.1 &
                                                 + !is.na(fc.kegg.p$greater[, "q.val"])]
      
      if(length(fc.kegg.p) > 2){
        path.ids.l <- rownames(fc.kegg.p$less)[fc.kegg.p$less[, "q.val"] < 0.1 &
                                                 + !is.na(fc.kegg.p$less[,"q.val"])]
        path.ids.2 <- substr(c(path.ids.2[1:3], path.ids.l[1:3]), 1, 8)
        
        #visualize pathway  
        pv.out.list <- sapply(na.omit(path.ids.2[1:6]), function(pid) pathview( gene.data =  gene.data, 
                                                                                pathway.id = pid, out.suffix=diff.tool))}
      
      
    }
    
  }
  
  
  
  
  #plot data
  
  if (plot.gene.data==T ){#& is.null(pathway.id))  ){
    gs=unique(unlist(gsets[rownames(fc.kegg.p$greater)[1:3]]))
    essData=essGene(gs,gene.data , ref =ref, samp =samp)
    for (gs in rownames(fc.kegg.p$greater)[1:3]) {
      outname = gsub(" |:|/", "_", substr(gs, 10, 100))
      geneData(genes = gsets[[gs]], exprs = essData, ref = ref,
               samp = samp, outname = outname, txt = T, heatmap = T,
               Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
    }
    if(length(fc.kegg.p) > 2){
      gs=unique(unlist(gsets[rownames(fc.kegg.p$lesser)[1:3]]))
      essData=essGene(gs,gene.data , ref =ref, samp =samp)
      for (gs in rownames(fc.kegg.p$lesser)[1:3]) {
        outname = gsub(" |:|/", "_", substr(gs, 10, 100))
        geneData(genes = gsets[[gs]], exprs = essData, ref = ref,
                 samp = samp, outname = outname, txt = T, heatmap = T,
                 Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
      }
      
    }
  }
}

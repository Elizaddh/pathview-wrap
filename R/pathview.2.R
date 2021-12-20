#' @export

pathview.2 <- function( run ,diff.tool, gene.data= cnts, ref, samp, gsets, pathway.id,  both.dirs = list(gene = T, cpd = T))
{
    #library(pathview)
    if (run =="complete")
    {
        logfoldchange <- rundifftool(diff.tool, gene.data, ref, samp)
        print("diff tool run successful")
        path.ids <- run_path_analysis(logfoldchange, gsets )
        print("gage run successful")
        print("now pathview")
        if ( length(path.ids) <= 6 )
        {pv.out.list <- sapply(path.ids, function(pid) pathview::pathview( gene.data = logfoldchange, 
                                    pathway.id = pid, out.suffix=diff.tool))}
        if(length(path.ids) > 6)
        {pv.out.list <- sapply(path.ids[1:6], function(pid) pathview::pathview( gene.data = logfoldchange, 
                                                                 pathway.id = pid, out.suffix=diff.tool))}
     
    }
    
    else {
        if(is.null(pathway.id)){
            print("starting with gage analysis")
            path.ids <- run_path_analysis(gene.data, gsets )
            print("now pathview")
            if ( length(path.ids) <6 )
            {pv.out.list <- sapply(path.ids, function(pid) pathview::pathview( gene.data = gene.data, 
                                                                     pathway.id = pid, out.suffix=diff.tool))}
            if(length(path.ids) > 6)
            {
                print("pids are more than 6")
                pv.out.list <- sapply(path.ids[1:6], function(pid) pathview::pathview( gene.data =gene.data, 
                                                                          pathway.id = pid,species = "hsa", out.suffix=diff.tool))}
            
        }
        else{
            pathview::pathview( gene.data = gene.data,  pathway.id = pathway.id , out.suffix=diff.tool)}
            
            
        }
        
        
    }
    

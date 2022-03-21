#' this is a function to run different diff expression analysis tools
#' @import DESeq2
#' @import edgeR
#' @import limma
#' @import EnhancedVolcano
#'
#' @export
rundifftool <- function(diff.tool, gene.data, ref, samp, outname){
 
    grp.idx = NULL
    grp.idx[ref] <- "reference"
    grp.idx[samp] <- "sample"
    if (diff.tool %in% c("deseq2", "limma", "edgeR"))
    {
   if (diff.tool == "deseq2")
    {
        library(DESeq2)
        print("Deseq2 is running")
        coldat=DataFrame(grp=factor(grp.idx))
        dds <- DESeq2::DESeqDataSetFromMatrix(gene.data, colData=coldat, design = ~ grp)
        dds <- DESeq2::DESeq(dds)
        deseq2.res <- results(dds)
        #direction of fc, depends on levels(coldat$grp), the first level
        #taken as reference (or control) and the second one as experiment.
        deseq2.fc=deseq2.res$log2FoldChange
        names(deseq2.fc)=rownames(deseq2.res)
        exp.fc=deseq2.fc
        #pdf("Volcano_deseq2.pdf", width = 14,height= 14)
        tiff(paste0(outname,"Volcano_deseq2.tiff"), units="in", width=15, height=15, res=300)
        plot(EnhancedVolcano::EnhancedVolcano(deseq2.res, x = 'log2FoldChange', y = 'pvalue', lab = rownames(deseq2.res)))
        #plot(x = 1:10, y = 1:10)
        dev.off()
        print("this is the png file")
   }
    if (diff.tool == "edgeR"){
    
    dgel <- edgeR::DGEList(counts=gene.data, group=factor(grp.idx))
    dgel <- edgeR::calcNormFactors(dgel)
    dgel <- edgeR::estimateCommonDisp(dgel)
    dgel <- edgeR::estimateTagwiseDisp(dgel)
    et <- edgeR::exactTest(dgel)
    edger.fc=et$table$logFC
    names(edger.fc)=rownames(et$table)
    exp.fc=edger.fc
    tiff(paste0(outname,"Volcano_edgeR.tiff"), units="in", width=15, height=15, res=300)
    plot(EnhancedVolcano::EnhancedVolcano(et$table, x = 'logFC', y = "PValue", lab = rownames(et$table)))
     dev.off()
    
    }
    if(diff.tool == "limma")
        {
        
       
       dgel2 <- edgeR::DGEList(counts=gene.data, group=factor(grp.idx))
       dgel2 <- edgeR::calcNormFactors(dgel2)
       library(limma)
       design <- limma::model.matrix(~grp.idx)
       log2.cpm <- limma::voom(dgel2,design)
       fit <- limma::lmFit(log2.cpm,design)
       fit <- limma::eBayes(fit)
       limma.res=limma::topTable(fit,coef=2,n=Inf,sort="p")
       limma.fc=limma.res$logFC
 
       names(limma.fc)=limma.res$ID
       exp.fc=limma.fc
    }
        
    }
    else
    {
        print("The diff tool is not avaliable")
    }
    return(exp.fc)
   
}

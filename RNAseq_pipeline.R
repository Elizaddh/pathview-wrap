#!/usr/bin/env Rscript
my_pipeline <- function(fq.dir  = "/Users/edhungel/Documents/Research/UNCC/mouse/fresh/mouse_raw",  ref.dir = NA, phenofile = "/Users/edhungel/Documents/Research/UNCC/mouse/fresh/mouse_raw/pheno.txt", outdir = "/Users/edhungel/Documents/Research/UNCC/mouse/fresh/mouse_raw/results", endness="SE",  species = "Mus musculus",corenum = 2, diff.tool = "DESEQ2", compare = "unpaired"){ 
	library(stringr)

	#################################################################
	##
	#Check if files/folders  exists and create if not 
	##
	#################################################################

	if (!file.exists(outdir)){ 
		# default output file
		dir.create(outdir)
	}
	result.dir <- outdir
	#print("The results will be organized in ",result.dir)

	# make sure the second column is class and first column is sample name 
	# make sure file is tab seperated
	if (!file.exists(phenofile)){ ###TO DO make sure reference is first aplhanumerically# 
		print("Please provide phenofile with Class information")
	}
	coldata <- read.table(phenofile, sep = "\t", header = T)
	if(colnames(coldata)[2]!="Class"){
		print("Please make sure class information is in cloumn 2 with colname 'Class' . ")
	}
	coldata$Class <- as.factor(coldata$Class)
	ref <- which(coldata[, 2] ==  levels(coldata[, 2])[1])
	samp <- which(coldata[, 2] ==  levels(coldata[, 2])[2])

	##TO DO write something to automatically determine paired infromation, rev/fr etc

	#check and create dir for organizing results
	checkcretdir <- function(parentname, dirname){
		if(!file.exists(file.path(parentname, dirname))) {
			dir.create(file.path(parentname, dirname))
		} 
		assign(dirname,value = file.path(parentname, dirname), envir = .GlobalEnv)
	}

	folder_to_create<- list("fastqc_results", "trimmomatic_results","gage_results", "differential_analysis","aligned_bam","pathway_analysis" )
	trim_dir <-list ( "trimmomatic_log", "unpaired")
	diff_dir <-list ("DESeq2","edgeR")
	pathway_types <- list("KEGG", "GO")
	kegg_types <- list("signalling", "metabolism", "disease", "sig_n_met")
	go_types <- list("biological_process", "molecular_function", "cellular_component")
	lapply(folder_to_create, checkcretdir, parentname= result.dir  )
	lapply(trim_dir, checkcretdir, parentname= file.path(result.dir ,"trimmomatic_results")  )
	lapply(diff_dir, checkcretdir, parentname= file.path(result.dir , "differential_analysis")  )
	lapply(pathway_types, checkcretdir, parentname= file.path(result.dir , "gage_results")  )
	lapply(kegg_types, checkcretdir, parentname= file.path(result.dir ,"gage_results","KEGG" )  )
	lapply(go_types, checkcretdir, parentname= file.path(result.dir ,"gage_results","GO" )  )

	#just to make sure rest of codes are same
	qc.dir <- fastqc_results
	diff.dir <- differential_analysis
	trim.dir <- trimmomatic_results
	gage.dir <- gage_results
	trim.log <- trimmomatic_log
	pathway.dir <- pathway_analysis
	edger.dir <- edgeR
	deseq2.dir <- DESeq2
	kegg.dir <- KEGG
	go.dir <- GO


	#References
	#if only species name is given and both geneAnnotation and genome is NULL
	if( is.na(ref.dir)){
		ref_info <- read.table("data/species_genome_annotation_pkg", sep = "\t", header = T, na.strings=c(""," ","NA")) #this file is supplied with script

		species_no <- which(ref_info$species==species)
		annotate_pkg <- ref_info$annotation[species_no]
		genome_pkg <- ref_info$genome[species_no]

	###make sure both annot and genome package is installed for the species
	# (set of genome and annotation pkg come from developers list)

		pkg.on = require(annotate_pkg, character.only = TRUE, lib.loc = .libPaths()[1])
		if (!pkg.on) {
			if (!requireNamespace("BiocManager", quietly=TRUE))
				install.packages("BiocManager")
			BiocManager::install(annotate_pkg, suppressUpdates =TRUE, lib.loc = .libPaths()[1] )
			pkg.on = require(annotate_pkg, character.only = TRUE, lib.loc = .libPaths()[1])
			if (!pkg.on)
				stop(paste("Fail to install/load gene annotation package ",annotate_pkg, "!", sep = ""))
		}
		geneAnnotation <-  file.path(.libPaths()[1],annotate_pkg, "extdata"  , paste0(annotate_pkg, ".sqlite" ) ) 
		genomeFile <- genome_pkg
	} else {
		genomeFile <- list.files(ref.dir, ".fa$", full.names= T)
		geneAnnotation <- list.files(ref.dir, ".gtf$", full.names = T) #could be changed to include one of gtf, gff etc, check with quasR package
		}


	### To run qAlign we need samplefile
	## TO DO Make sure this works for all types of file name and single and paired end data, bunch of bam files and bunch of fastq files, partially complete
	##############################################################################

	if( endness== "SE"){
		pinfo_string <- ".fastq"
	}else{
		pinfo_string <- "_1.fastq"
	}
	library(stringr)
	FileName <- grep(pinfo_string,list.files(fq.dir, full.names=T) ,value =T)
	FileName <- str_replace_all(file.path(trim.dir,  basename(FileName)),pinfo_string, paste0("_paired", pinfo_string))
	sampleFile <- file.path(result.dir, "sampleFile.txt")
	SampleName <-  str_remove_all(basename(FileName), paste0("_paired",pinfo_string))
	if(endness == "SE") {
	write.table(file =sampleFile,sep = "\t", as.data.frame( cbind(FileName, SampleName)) ,quote =F ,  col.names=T, row.names=F)
	} else{
		FileName1 <- FileName
		FileName2 <- str_replace_all(FileName1, "_1.fastq", "_2.fastq")
		sampleFile <- file.path(result.dir, "sampleFile.txt")
		write.table(file =sampleFile,sep = "\t", as.data.frame( cbind(FileName1, FileName2, SampleName)) ,quote =F ,  col.names=T, row.names=F)
	}

	###############################################################################
	##Make sure all these libraries are actually needed
	library(ggplot2)
	library(Rsamtools)
	library(GenomicFeatures)
	library(QuasR)
	library(parallel)
	library(fastqcr)
	library(stringr)
	##to detect how many cores a machine has
	##detectCores() can be used but better let user suppy corenumbers
	#
	#
	#############################################################################
	#2. RUN THE ANALYSIS # FASTQC
	############################################################### 
	#fastqc(fq.dir, qc.dir, threads =8,fastqc.path =system( "which fastqc", intern = T))
	qc <- qc_aggregate(qc.dir)
	#work here #TO DO  make sure figures look good
	tiff(file.path(qc.dir,"total_seq.tiff") ,units="in", width=15, height=15, res=300)
	ggplot(qc, aes(x=sample, y=tot.seq))+ geom_bar(stat="identity")  + coord_flip()
	dev.off()
	tiff(file.path(qc.dir,"qc_heatmap.tiff"), units="in", width=15, height=15, res=300)
	ggplot(qc, aes(x=sample, y =status)) + geom_tile()
	dev.off()
	print("Multiqc report")
	## Building Multi QC Reports # Error: pandoc version 1.12.3 or higher is required and was not found (see the help page ?rmarkdown::pandoc_available).
	#qc_report(qc.dir, result.file =file.path(qc.dir,"multi-qc-report"))
	#system("multiqc qc.dir") # html are better

	#############################################################################
	#2. RUN THE ANALYSIS # QUALITY TRIMMING 
	############################################################################# 
	#TO DO Download trimmomatic and use that file path to do the processing #make sure right adapters are used

	run_trimmomaticR <-function(samplename){
		print("trim function is called")
		if (endness=="PE"){
			cmd <- paste0("java -jar tools/Trimmomatic-0.39/trimmomatic-0.39.jar  PE " , file.path(fq.dir , "samplename_to_sed_1.fastq"), " ", file.path(fq.dir , "samplename_to_sed_2.fastq"), " ",file.path(trim.dir , "samplename_to_sed_paired_1.fastq"), " ",file.path(trim.dir , "unpaired/samplename_to_sed_unpaired_1.fastq"), " ",file.path(trim.dir , "samplename_to_sed_paired_2.fastq"), " ",file.path(trim.dir , "unpaired/samplename_to_sed_unpaired_2.fastq"),   " ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:4:20 TRAILING:5 MINLEN:25 >> ", file.path(trim.log , "trimmomatic.samplename_to_sed.log"))
		} else {
			cmd <- paste0("java -jar tools/Trimmomatic-0.39/trimmomatic-0.39.jar  SE " ,  file.path(fq.dir , "samplename_to_sed.fastq"), " ",file.path(trim.dir , "samplename_to_sed_paired.fastq"),  " ILLUMINACLIP:TruSeq3-SE:2:30:10 SLIDINGWINDOW:10:20 >> ", file.path(trim.log , "trimmomatic.samplename_to_sed.log") ) #may be change outfile name, can be miss leading
		}
		cmd <- stringr::str_replace_all(cmd, "samplename_to_sed", samplename)
		system(cmd)
	}

	#call function for quality trimming
	cl <- makeCluster(corenum)
	clusterExport(cl,c("fq.dir","endness", "trim.dir", "trim.log"), envir = environment())#.GlobalEnv)
	ans <- parSapply(cl , read.csv( sampleFile , header =T, sep ="\t")$SampleName  , run_trimmomaticR)
	print("the run trimm is complete")
	stopCluster(cl)
	#
	print("the cluster are done" )

	#############################################################################
	#2. RUN THE ANALYSIS # Alignmnet and counting 
	############################################################################# 

	cl2 <- makeCluster(corenum)
	print("Alignment is running")
	if (endness=="PE"){
		aligned_proj <- QuasR::qAlign(sampleFile,paired ="fr", clObj=cl2, alignmentsDir =aligned_bam ,  genome=genomeFile,geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner ="Rhisat2" ) 
	print("align done")
	} else {

		aligned_proj <- QuasR::qAlign(sampleFile,paired ="no", clObj=cl2, alignmentsDir =aligned_bam ,  genome=genomeFile,geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner ="Rhisat2" )
	}
	saveRDS(aligned_proj, file.path(result.dir , "alltrimmedalignedobj.RDS"))

	##move this chunck to up where u have processed annotation and genome file##
	txdb <- try(loadDb(geneAnnotation), silent = T)
	if (class(txdb)==  "TxDb"){
		if(!grepl("chr", seqlevels(txdb)[1])){ #check if this is necessary
		 newSeqNames <- paste('Chr', seqlevels(txdb), sep = '')
		 names(newSeqNames) <- seqlevels(txdb)
		 txdb <- renameSeqlevels( txdb, newSeqNames )
		 seqlevels(txdb)
	}
	}else{
	library(GenomicFeatures)

	chrLen <- scanFaIndex(genomeFile)
	chrominfo <- data.frame(chrom = as.character(seqnames(chrLen)),
				length = width(chrLen),
				is_circular = rep(FALSE, length(chrLen)))
	txdb <- makeTxDbFromGFF(file = geneAnnotation, format = "gtf",
				chrominfo = chrominfo,
				dataSource = "Ensembl",
				organism = species)
	}
	geneLevels <- qCount(aligned_proj, txdb, reportLevel ="gene", clObj=cl2)
	saveRDS(geneLevels, file.path(result.dir, "combinedcount.trimmed.RDS"))


	library(dplyr)
	library(tibble)
	library(biomaRt)

	###########################################################################################################################################
	## RUN THE ANALYSIS C. DIFFERENTIAL ANALYSIS
	###########################################################################################################################################
	library(DESeq2)
	#setwd("/scratch/edhungel/Rsubread/script")
	geneData_my <- as.data.frame(readRDS(file.path(result.dir, "combinedcount.trimmed.RDS"))) #TO DO maybe I can just use geneLevels variable

	library(gage)
	library(pathview)
	cnts <- geneData_my[,-1]
	kegg.gs.species <- kegg.gsets(species)
	orgcode<- kegg.species.code(species)
	data(bods)
	if(!all(rownames(cnts)%in% unlist(unname(kegg.gs.species$kg.sets)))){ #check if the use of "all" is appropriate
	  rownames(cnts)<- str_remove(rownames(cnts),"\\.[0-9]+$" )
	  cnts<- mol.sum(cnts, id.map = "ENSEMBL", gene.annotpkg =bods[ which(bods[,3]==orgcode)]) #converting to entrez # what if gene id is not ensembl and what if arabidopsis thaliana id.map might be ath or else thing
	}

	grp.idx <-NULL
	grp.idx[ref] <- "reference"
	grp.idx[samp] <- "sample"

	if(diff.tool == "DESEQ2"){
		coldat=DataFrame(grp=factor(grp.idx))
		dds <- DESeqDataSetFromMatrix(cnts, colData=coldat, design =~ grp)
		dds <- DESeq(dds)
		deseq2.res <- results(dds)
		#direction of fc, depends on levels(coldat$grp), the first level
		#taken as reference (or control) and the second one as experiment.
		deseq2.fc=deseq2.res$log2FoldChange
		names(deseq2.fc)=rownames(deseq2.res)
		exp.fc=deseq2.fc
		table(is.na(deseq2.res$padj))
		write.table(deseq2.res , file.path(deseq2.dir, "DESEQ2_logfoldchange.txt"), col.names =TRUE, row.names =TRUE, quote =FALSE)
		tiff(file.path(deseq2.dir, "Volcano_deseq2.tiff"), units="in", width=15, height=15, res=300)
		plot(EnhancedVolcano::EnhancedVolcano(deseq2.res, x ='log2FoldChange', y ='pvalue', lab =rownames(deseq2.res)))
		dev.off()
	}
	if(diff.tool=="edgeR"){
	###########################################################################################################################################
		library(edgeR)
		dgel <- edgeR::DGEList(counts=cnts, group=factor(grp.idx))
		dgel <- edgeR::calcNormFactors(dgel)
		dgel <- edgeR::estimateCommonDisp(dgel)
		dgel <- edgeR::estimateTagwiseDisp(dgel)
		et <- edgeR::exactTest(dgel)
		edger.fc=et$table$logFC
		names(edger.fc)=rownames(et$table)
		exp.fc=edger.fc
		write.table(et , file.path(edger.dir, "EDGER_logfoldchange.txt"), col.names =TRUE, row.names =TRUE, quote =FALSE)
		tiff(file.path(edger.dir, "edgeR_Volcano_edgeR.tiff"), units="in", width=15, height=15, res=300)
		plot(EnhancedVolcano::EnhancedVolcano(et$table, x ='logFC', y="PValue", lab=rownames(et$table)))
		dev.off()
	}
	###########################################################################################################################################
	#RUN THE ANALYSIS # PATHWAY ANALYSIS AND VISUALIZATON
	###########################################################################################################################################
	library(gage)

	kegg.gs<- kegg.gsets(species)
	###########################################################################################################################################
	#TO DO
	#make sure the gene id of gsets match with cnts genes id #partial complete , ensemble will be converted to entrez
	###########################################################################################################################################
	#TO THINK and futher test
	#maybe run gage with cnts also  data and do analysis 
	####################################################################################################################
	run_gset_analysis <- function(gsets, work.dir, same.dir , compare){
		setwd(work.dir)
		
		fc.kegg.p <- gage(exp.fc, gsets = gsets, ref = NULL, samp = NULL, same.dir = same.dir, compare = compare)
		sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
		path.ids <- rownames(fc.kegg.p$greater)[sel]
		write.table(fc.kegg.p$greater, file = file.path(work.dir , "fc.kegg.p.greater.txt"), sep = "\t")
		if (same.dir == T){
			sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
			path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
			write.table(fc.kegg.p$less, file = file.path(work.dir , "fc.kegg.p.less.txt"), sep = "\t")
		        path.ids <- c(path.ids[1:3], path.ids.l[1:3])
	
		}

		path.ids <- substr(path.ids, 1, 8)
		#path.ids <- gsub("[^0-9.-]", "", sapply(stringr::str_split(path.ids, " ", 2),"[[",1))
		#visualize top 3 pathways
		#run pathview only for KEGG pathways
		if (same.dir==F){
			pv.out.list <- sapply(na.omit(path.ids[1:6]), function(pid) pathview( gene.data = exp.fc, pathway.id = pid, species = species, out.suffix=diff.tool))
		}
		#kegg.sig<-sigGeneSet(fc.kegg.p, outname=paste0(species,"kegg.sig",basename(work.dir)), pdf.size=c(10,10), heatmap = T ) #wont give heatmap
		#write.table(kegg.sig$greater, file = file.path(gage.dir , "kegg.sig.txt"), sep = "\t")
		 #kegg.esg.up <- esset.grp(fc.kegg.p$greater, exp.fc, gsets = kegg.gs$kg.sets[path_type_ids], ref = NULL, samp = NULL, test4up = T, output = T, outname = "kegg.up", make.plot = T) # currently heatmap is not supported 
		 gs=unique(unlist(gsets[rownames(fc.kegg.p$greater)[1:3]]))
		 essData=essGene(gs, cnts, ref =NULL, samp =NULL)
		 head(essData, 4)
		 ref1=ref
		 samp1=samp
		 for (gs in rownames(fc.kegg.p$greater)[1:3]) {
		 outname = gsub(" |:|/", "_", substr(gs, 10, 100))
		 geneData(genes = kegg.gs[[gs]], exprs = essData, ref =ref, samp = samp, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T) 
		}
		 for (gs in rownames(fc.kegg.p$greater)[1:3]) {
		 outname = gsub(" |:|/", "_", substr(gs, 10, 100))
		 outname = paste(outname, "all", sep=".")
		 geneData(genes = kegg.gs$kegg.gs[[gs]], exprs = cnts, ref = ref, samp = samp, outname =outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
		 }

	}

	#TO DO try using lapply for all function call of kegg pathways

	run_gset_analysis(kegg.gs$kg.sets[kegg.gs$sigmet.idx],sig_n_met,same.dir = F, compare=compare ) 
  	run_gset_analysis(kegg.gs$kg.sets[kegg.gs$dise.idx], disease, same.dir = F, compare=compare)
  	run_gset_analysis(kegg.gs$kg.sets[kegg.gs$sig.idx], signalling, same.dir = F, compare=compare)
  	run_gset_analysis(kegg.gs$kg.sets[kegg.gs$met.idx], metabolism, same.dir = F,compare=compare)
	
	############################################################################################  

	common_name_species <- korg[korg[,4]==species][5]
	go.gs <- go.gsets(common_name_species)
	go.bp<- go.gs$go.sets[go.gs$go.subs$BP]
	go.mf<-go.gs$go.sets[go.gs$go.subs$MF]
	go.cc<- go.gs$go.sets[go.gs$go.subs$CC]


	run_gset_analysis(go.bp,biological_process,same.dir = T, compare=compare ) 
  	run_gset_analysis(go.mf,molecular_function, same.dir = T, compare=compare)
  	run_gset_analysis(go.cc, cellular_component, same.dir = T, compare=compare)
}

##make your own function call here
my_pipeline (fq.dir  = "/Users/edhungel/Documents/Research/UNCC/mouse/fresh/mouse_raw/mouse_test",  ref.dir = "~/Documents/Research/UNCC/Data/Reference/mouse", phenofile = "/Users/edhungel/Documents/Research/UNCC/mouse/fresh/mouse_raw/mouse_test/pheno.txt",  outdir = "/Users/edhungel/Documents/Research/UNCC/mouse/fresh/mouse_raw/mouse_test/results2",  endness="SE",  species = "Mus musculus",corenum=8, diff.tool = "DESEQ2", compare = "unpaired")

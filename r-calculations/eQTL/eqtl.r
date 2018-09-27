emeraLD2R <- function(path, bin = "/Users/kevinjiang/Desktop/dev/nci-webtools-dceg-vQTL/r-calculations/eQTL/emeraLD"){
	require(data.table)
	if(!file.exists(bin)) stop(paste0("bin = '", bin, "' file does not exist"))
	function(region, matrix.out = TRUE, info.out = TRUE){
		require(data.table)
		opts <- paste(c("--matrix", "--stdout --extra")[c(matrix.out, TRUE)], collapse = " ")
		pc <- "| tr ':' '\t'"
		chr <- strsplit(region, ":")[[1]][1]
		vcfile <- gsub("\\$chr", chr, path)
		if(!file.exists(vcfile)) stop(paste0(vcfile, " does not exist"))
		out <- suppressMessages(fread(
			input  = paste(bin, "-i", vcfile, "--region", region, opts, pc), 
			header = FALSE, showProgress = FALSE
		))
		info <- NULL
		if(info.out){
			info <- out[,1:5]
			colnames(info) <- c("chr", "pos", "id", "ref", "alt")
		}
		if(matrix.out) out <- as.matrix(out[,-(1:5)]); colnames(out) <- NULL
		list("Sigma" = out, "info" = info)
	}
}

eqtl <- function(workDir, genoFile, exprFile, assocFile) {
  setwd(workDir)

  library(tidyverse)
  library(hrbrthemes)
  library(scales)
  library(ggrepel)
  library(forcats)
  library(jsonlite)

  gdatafile <- paste0('uploads/', genoFile)
  edatafile <- paste0('uploads/', exprFile)
  qdatafile <- paste0('uploads/', assocFile)

  gdata <- read_delim(gdatafile,delim = "\t",col_names = T)
  edata <- read_delim(edatafile,delim = "\t",col_names = T)
  qdata <- read_delim(qdatafile,delim = "\t",col_names = T,col_types = cols(variant_id='c'))
  qdata <- qdata %>% arrange(pval_nominal,desc(abs(slope)),abs(tss_distance)) %>% group_by(gene_id,variant_id,rsnum,ref,alt) %>% slice(1) %>% ungroup()

  chromosome <- qdata$chr[1]

  rcdatafile <- paste0('Recombination_Rate_CEU/CEU-',chromosome,'-final.txt.gz')
  rcdata <- read_delim(rcdatafile,delim = "\t",col_names = T)
  colnames(rcdata) <- c('pos','rate','map','filtered')
  rcdata$pos <- as.integer(rcdata$pos)

  # calculate gene expression boxplots
  tmp <- qdata %>% 
    group_by(gene_id,gene_symbol) %>% 
    arrange(pval_nominal) %>% 
    slice(1) %>% 
    ungroup() %>% 
    arrange(pval_nominal) %>% 
    slice(1:15)

  edata_boxplot <- edata %>% gather(Sample,exp,-(chr:gene_id)) %>% 
    right_join(tmp %>% select(gene_id,gene_symbol) %>% 
    unique())

  edata_boxplot <- edata_boxplot %>% 
    left_join(edata_boxplot %>% group_by(gene_id) %>% summarise(mean=mean(exp))) %>% 
    left_join(tmp %>% select(gene_id,pval_nominal)) %>% 
    mutate(gene_symbol=fct_reorder(gene_symbol,(pval_nominal)))

  gene_expression_data <- list(setNames(as.data.frame(edata_boxplot),c("chr","start","end","gene_id","Sample","exp","gene_symbol","mean","pval_nominal")))

  # calculate locus zoom plot

  tmp <- qdata %>% arrange(pval_nominal,desc(abs(slope)),abs(tss_distance)) %>% slice(1)
  default_gene <- tmp$gene_id
  default_vairnat <- tmp$variant_id
  default_rsnum <- tmp$rsnum
  defaul_info <-tmp %>% select(gene_id:alt)

  qdata_region <- qdata %>% filter(gene_id==default_gene)
  rcdata_region <- rcdata %>% filter(pos<=max(qdata_region$pos),pos>=min(qdata_region$pos))
  qdata_top_annotation <- qdata_region %>% filter(variant_id==default_vairnat)

  # source('emeraLD2R.r')
  in_path <- '/Users/kevinjiang/Desktop/dev/nci-webtools-dceg-vQTL/r-calculations/eQTL/chr1_149039120_152938045.vcf.gz'
  regionLD <- paste0(chromosome,":",min(qdata_region$pos),"-",max(qdata_region$pos))
  getLD <- emeraLD2R(path = in_path) 
  ld_data <- getLD(region=regionLD)

  index <- which(ld_data$info$id==default_rsnum|str_detect(ld_data$info$id,paste0(";",default_rsnum))|str_detect(ld_data$info$id,paste0(default_rsnum,";")))
  ld_info <- as.data.frame(ld_data$Sigma[,index])
  colnames(ld_info) <- "R2"
  rownames(ld_info) <- ld_data$info$id

  qdata_region$R2 <- (ld_info[qdata_region$rsnum,"R2"])^2

  locus_zoom_data <- list(setNames(as.data.frame(qdata_region),c("gene_id","gene_symbol","variant_id","rsnum","chr","pos","ref","alt","tss_distance","pval_nominal","slope","slope_se","R2")))

  # return outputs in list
  dataSource <- append(gene_expression_data, locus_zoom_data)
  return(dataSource)
}

# function call parameters will be inputted by R Wrapper
eqtl(workingDirectory, genotypeFile, expressionFile, associationFile)
library(tidyverse)
library(ggplot2)
library(cowplot)
# library(hrbrthemes)
library(scales)
library(ggsci)
library(ggrepel)

# sed Delete / Remove ^M Carriage Return
formatM_input<- function(tfile){
  filetype <- summary(file(tfile))$class
  if(filetype == "gzfile"){
    cmd1=paste0("bgzip -d -c ",tfile,' > ',tfile,'.txt')  
    cmd2=paste0("sed -i '' 's/\r//' ",tfile,'.txt')
    cmd3=paste0('cat ',tfile,'.txt |bgzip > ',tfile,' &&  rm -rf ',tfile,'.txt')
    system(cmd1)
    system(cmd2)
    system(cmd3)
  }else {
    cmd=paste0("sed -i '' 's/\r//' ",tfile)  
    system(cmd)
    
  }
}



coloc_QC <- function(gwasfile=NULL,gwasfile_pub=FALSE, qtlfile=NULL, qtlfile_pub=FALSE, ldfile=NULL,ldfile_pub=FALSE, leadsnp=NULL,leadpos=NULL, distance=100000,zscore_gene=NULL,output_plot_prefix=NULL, output_prefix="./ezQTL_input", logfile="ezQTL.log"){
  
  ainfo <- paste0("ezQTL Analysis Starting QC at ",Sys.time())
  cat(ainfo,file=logfile,sep="\n",append = FALSE)
  
  if(is.null(gwasfile) & is.null(qtlfile) & is.null(ldfile)){
    errinfo <- "ezQTL need at least one of the following files: QTL association file, GWAS association file, and LD matrix file;"
    cat(errinfo,file=logfile,sep="\n",append = T)
    stop("ezQTL QC failed: no file input")
  }
  
  # check gwas file
  if(!is.null(gwasfile)){
    # format file, remove the CR/^M characters from windows 
    if(!gwasfile_pub){ formatM_input(gwasfile) }
    # read orignal file 
    gwas <- read_delim(gwasfile,delim = '\t',col_names = T)
    # check file format, and select only requeired columns
    gwas_colnames <- c('chr','pos','ref','alt','rsnum','pvalue','zscore','effect','se')
    gwas_colname_diff <- gwas_colnames[!gwas_colnames %in% colnames(gwas)]
    if(length(gwas_colname_diff)!=0){
      errinfo <- paste0("ERROR: The following GWAS column names for ezQTL are not existed: ",paste0(gwas_colname_diff,collapse = ', '))
      cat(errinfo,file=logfile,sep="\n",append = T)
      stop("ezQTL QC failed: GWAS data format")
    }else{
      gwas <-gwas %>% select(one_of(gwas_colnames))
      gwas_nchr <- gwas %>% count(chr) %>% dim() %>% .[[1]]
      if(gwas_nchr!=1){
        errinfo <- "ERROR: Number of chromosome in GWAS file large than 1"
        cat(errinfo,file=logfile,sep="\n",append = T)
        stop("ezQTL QC failed: GWAS data format")
      }
      
      cat("\nGWAS summary",file=logfile,sep="\n",append = T)
      gwastmp <- gwas %>% select(rsnum,chr,pos,ref,alt) %>% unique() 
      cat(paste0('# number of variants included: ',dim(gwastmp)[1]),file=logfile,sep="\n",append = T)
      gwastmp <- gwastmp %>% filter(str_detect(rsnum,'^rs'))
      cat(paste0('# number of variants with rsnum: ',dim(gwastmp)[1]),file=logfile,sep="\n",append = T)
      gwastmp <- gwastmp  %>% filter((ref %in% c('A','T') & alt %in% c('A','T'))|(ref %in% c('C','G') & alt %in% c('C','G')))
      cat(paste0('# number of variants with rsnum and A/T, C/G alleles: ',dim(gwastmp)[1]),file=logfile,sep="\n",append = T)
    }
    
  }
  
  # check qtl file
  if(!is.null(qtlfile)){
    # format file, remove the CR/^M characters from windows 
    if(!qtlfile_pub){ formatM_input(qtlfile) }
    # read orignal file 
    qtl <- read_delim(qtlfile,delim = '\t',col_names = T,col_types = cols('variant_id'='c'))
    # check file format, and select only requeired columns
    qtl_colnames <- c('gene_id','gene_symbol','variant_id','rsnum','chr','pos','ref','alt','tss_distance','pval_nominal','slope','slope_se')
    qtl_colname_diff <- qtl_colnames[!qtl_colnames %in% colnames(qtl)]
    if(length(qtl_colname_diff)!=0){
      errinfo <- paste0("ERROR: The following QTL column names for ezQTL are not existed: ",paste0(qtl_colname_diff,collapse = ', '))
      cat(errinfo,file=logfile,sep="\n",append = T)
      stop("ezQTL QC failed: QTL data format")
    }else{
      qtl <-qtl %>% select(one_of(qtl_colnames))
      
      qtl_nchr <- qtl %>% count(chr) %>% dim() %>% .[[1]]
      if(qtl_nchr!=1){
        errinfo <- "ERROR: Number of chromosome in QTL file large than 1"
        cat(errinfo,file=logfile,sep="\n",append = T)
        stop("ezQTL QC failed: QTL data format")
      }
      
      cat("\nQTL summary",file=logfile,sep="\n",append = T)
      ntraits <- qtl %>% count(gene_id) %>% dim() %>% .[[1]]
      cat(paste0('# number of QTL traits: ',ntraits),file=logfile,sep="\n",append = T)
      qtltmp <- qtl %>% select(variant_id,rsnum,chr,pos,ref,alt) %>% unique() 
      cat(paste0('# number of variants included: ',dim(qtltmp)[1]),file=logfile,sep="\n",append = T)
      qtltmp <- qtltmp %>% filter(str_detect(rsnum,'^rs'))
      cat(paste0('# number of variants with rsnum: ',dim(qtltmp)[1]),file=logfile,sep="\n",append = T)
      qtltmp <- qtltmp  %>% filter((ref %in% c('A','T') & alt %in% c('A','T'))|(ref %in% c('C','G') & alt %in% c('C','G')))
      cat(paste0('# number of variants with rsnum and A/T, C/G alleles: ',dim(qtltmp)[1]),file=logfile,sep="\n",append = T)
      
      
      
      if(all(is.na(qtl$gene_symbol)) | all(qtl$gene_symbol=="NA")){
        qtl <- qtl %>% mutate(gene_symbol=gene_id)
      }else{
        if(any(is.na(qtl$gene_symbol))){
          qtl <- qtl %>% mutate(gene_symbol=paste0(gene_id,'|',gene_symbol))
        }
      }
      
      dup_qtl_id <- qtl %>% count(gene_id,gene_symbol) %>% count(gene_symbol) %>% filter(n>1)
      if(dim(dup_qtl_id)[1]>1){
        cat(paste0('Duplicated gene symbol in the qtl files. Use the gene id and gene symbol together as the new id'),file=logfile,sep="\n",append = T)
        qtl <- qtl %>% mutate(gene_symbol=paste0(gene_id,"|",gene_symbol))
      }
      
    }
    
  }
  
  # check ld file
  if(!is.null(ldfile)){
    # format file, remove the CR/^M characters from windows 
    if(!ldfile_pub){ formatM_input(ldfile) }
    # read orignal file 
    ld.matrix <- read_delim(ldfile,delim = '\t',col_names = F,col_types = cols('X1'='c')) %>% rename(chr=X1,pos=X2,rsnum=X3,ref=X4,alt=X5)
    # check file format, and select only requeired columns
    ld.matrix.size=dim(ld.matrix)
    if(ld.matrix.size[1]<7 | ld.matrix.size[2]<7 | (ld.matrix.size[2]-ld.matrix.size[1])!=5){
      errinfo <- paste0("ERROR: The dimensions of the LD matrix file are not correct. Please check the help page for the detail.")
      cat(errinfo,file=logfile,sep="\n",append = T)
      stop("ezQTL QC failed: LD data format")
    }else{
      # pre-porcessing data
      ld.info <- ld.matrix %>% select(chr,pos,rsnum,ref,alt) %>% mutate(Seq=seq_along(chr))
      ld.matrix <- ld.matrix %>% select(-c(chr,pos,rsnum,ref,alt)) %>% as.matrix
      
      cat("\nLD summary",file=logfile,sep="\n",append = T)
      ldinfotmp <- ld.info %>% select(rsnum,chr,pos,ref,alt) %>% unique() 
      cat(paste0('# number of variants included: ',dim(ldinfotmp)[1]),file=logfile,sep="\n",append = T)
      ldinfotmp <- ldinfotmp %>% filter(str_detect(rsnum,'^rs'))
      cat(paste0('# number of variants with rsnum: ',dim(ldinfotmp)[1]),file=logfile,sep="\n",append = T)
      ldinfotmp <- ldinfotmp  %>% filter((ref %in% c('A','T') & alt %in% c('A','T'))|(ref %in% c('C','G') & alt %in% c('C','G')))
      cat(paste0('# number of variants with rsnum and A/T, C/G alleles: ',dim(ldinfotmp)[1]),file=logfile,sep="\n",append = T)
    }
  }
  
  
  if(!is.null(leadpos)){
    if(!is.null(gwasfile)){
      leadsnp <- gwas %>% filter(pos==leadpos)%>% pull(rsnum) %>% unique() 
    }else{
      if(!is.null(qtlfile)){
        leadsnp <- qtl %>% filter(pos==leadpos) %>% pull(rsnum) %>% unique()
      }
    }
  }
  
  if(is.null(leadsnp)){
    # default, using the most significant GWAS variant as the leadsnp
    if(!is.null(gwasfile)){
      leadsnp <- gwas %>% arrange(pvalue) %>% slice(1) %>% pull(rsnum)  
    }else{
      if(!is.null(qtlfile)){
        leadsnp <- qtl %>% arrange(pval_nominal) %>% slice(1) %>% pull(rsnum)  
      }
    }
  }
  
  if(!is.null(leadsnp))  {
    cat(paste0("\nReference SNP using for filtering SNPs based on the snp: ",leadsnp,' and distance: ',distance),file=logfile,sep="\n",append = T)
    if(!is.null(gwasfile)){
      leadchr <- gwas %>% filter(rsnum==leadsnp) %>% pull(chr)
      if(is.null(leadpos)) {leadpos <- gwas %>% filter(rsnum==leadsnp) %>% pull(pos)}
      
    }else{
      if(!is.null(qtlfile)) {
        leadchr <- qtl %>% filter(rsnum==leadsnp) %>% slice(1) %>% pull(chr) 
        if(is.null(leadpos)) {leadpos <- qtl %>% filter(rsnum==leadsnp) %>% slice(1) %>% pull(pos)}
      }
    }
    
    leadpos1 <- leadpos - distance
    leadpos2 <- leadpos + distance
    
    if(!is.null(gwasfile)) {gwas <- gwas %>% filter(chr==leadchr,pos>=leadpos1,pos<=leadpos2)}
    if(!is.null(qtlfile)) {qtl <- qtl %>% filter(chr==leadchr,pos>=leadpos1,pos<=leadpos2)}
    
    #ld.info <- ld.info %>% filter(chr==leadchr,pos>=leadpos1,pos<=leadpos2)
    if(!is.null(ldfile)){
      leadindex <- ld.info %>% filter(rsnum==leadsnp) %>% pull(Seq)
      #!is.na(leadindex) | 
      if(length(leadindex)!=0){
        ldtmp <- tibble(LD=ld.matrix[,leadindex])
        ld_leadsnp <- TRUE
      }else{
        ld_leadsnp <- FALSE
        ldtmp <- tibble(LD=rep(NA,dim(ld.matrix)[1]))
      }
      
      ld.data <- bind_cols(
        ld.info,
        ldtmp
      ) %>% filter(chr==leadchr,pos>=leadpos1,pos<=leadpos2)
    }
    
  }else {
    cat("No Reference SNP using for filtering SNPs based on the distance",file=logfile,sep="\n",append = T)
  }
  
  ## QTL plots
  if(!is.null(qtlfile)){
    # qtl minmal p-visualize
    qtl_min <- qtl %>% 
      group_by(gene_id,gene_symbol) %>%
      arrange(pval_nominal) %>%
      slice(1) %>%
      ungroup() %>% 
      arrange((pval_nominal)) %>% 
      slice(1:50) %>% 
      arrange(desc(pval_nominal))
    
    pall0 <- qtl %>% filter(gene_id %in% qtl_min$gene_id) %>% 
      mutate(gene_id=factor(gene_id,levels = qtl_min$gene_id)) %>% 
      ggplot(aes(-log10(pval_nominal),gene_id,fill=-log10(pval_nominal)))+
      geom_point(pch=21,stroke=0.2,col="black",size=3)+
      geom_point(data = qtl_min,pch=21,stroke=0.5,col="red",size=3)+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XYx",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill="minmal p-value",y='',x='QTL P-value (-log10)')+
      theme(legend.position = 'none',panel.background = element_blank())+
      scale_fill_viridis_c(direction = -1)+
      scale_x_continuous(breaks = pretty_breaks(),expand = expansion(mult = c(0.01, .12)))+
      ggrepel::geom_text_repel(data=qtl_min,aes(label=rsnum),hjust=0,vjust=1,nudge_y = -0.1,segment.size=0.2,segment.color="black",segment.curvature = -0.1,segment.ncp = 3,segment.angle = 30)+
      #coord_cartesian(clip = 'off')+
      panel_border(color = 'black')
    
    xleng <- 2+dim(qtl_min)[1]*0.5
    xleng <- if_else(xleng>12,12,xleng)
    ggsave(filename = paste0(output_plot_prefix,"_QC_QTLminP.svg"),plot = pall0,width = 12,height = xleng)
  }else{
    pall0 <- ggplot()
    ggsave(filename = paste0(output_plot_prefix,"_QC_QTLminP.svg"),plot = pall0,width = 12,height = 4)
  }
  
  
  # Summary of the dataset --------------------------------------------------
  
  # The following function work for all three major inputs ------------------
  if(!is.null(gwasfile) & !is.null(qtlfile) & !is.null(ldfile)){
    
    gwas <- gwas %>% mutate(chr=as.character(chr))
    qtl <- qtl %>% mutate(chr=as.character(chr))
    ld.info <- ld.info %>% mutate(chr=as.character(chr))
    
    snpall <- inner_join(
      gwas %>% select(rsnum,chr,pos_gwas=pos,ref_gwas=ref,alt_gwas=alt) %>% unique(),
      qtl %>% select(rsnum,chr,pos_qtl=pos,ref_qtl=ref,alt_qtl=alt) %>% unique()
    ) %>% 
      inner_join(
        ld.info %>% select(rsnum,chr,pos_ld=pos,ref_ld=ref,alt_ld=alt) %>% unique()
      )
    
    cat("\nOverlapped SNPs",file=logfile,sep="\n",append = T)
    cat(paste0('# number of overlapped snps by rsnum: ',dim(snpall)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpall %>% filter(pos_gwas!=pos_ld|pos_gwas!=pos_qtl|pos_qtl!=pos_ld | ref_gwas!=ref_ld|ref_gwas!=ref_qtl|ref_qtl!=ref_ld|alt_gwas!=alt_ld|alt_gwas!=alt_qtl|alt_qtl!=alt_ld)
    cat(paste0('# number of overlapped snps by rsnum, but have different pos, ref, or alt:  ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    if(dim(snpalltmp)[1]>0){
      pertmp <- dim(snpalltmp)[1]/dim(snpall)[1]
      pertmp2 <- percent_format(accuracy = 0.01)(pertmp)
      cat(paste0('\nWarning: found ',pertmp2,' SNPs with same rsnum, but different information. Please check the input datasets. Make sure input the data have the same genome build and allels information. Check the SNP Matching information for detail (snp_not_match.txt). However, ezQTL still can use the overlapped rsnum as the ID for colocalization analyses.'),file=logfile,sep="\n",append = T)
      snpalltmp %>% write_delim('snp_not_match.txt',delim = '\t',col_names = T)
      if(pertmp > 0.9){
        cat(paste0('Warning: align GWAS and QTL dataset based on LD information'),file=logfile,sep="\n",append = T)
        snpall <- snpall %>% select(rsnum,chr,pos=pos_ld,ref=ref_ld,alt=alt_ld)
        gwas <- gwas %>% filter(rsnum %in% snpall$rsnum) %>% select(-pos,-ref,-alt) %>% left_join(snpall) %>% select(one_of(gwas_colnames))
        qtl <- qtl %>% filter(rsnum %in% snpall$rsnum) %>% select(-pos,-ref,-alt) %>% left_join(snpall) %>% select(one_of(qtl_colnames))
        
        snpall <- inner_join(
          gwas %>% select(rsnum,chr,pos_gwas=pos,ref_gwas=ref,alt_gwas=alt) %>% unique(),
          qtl %>% select(rsnum,chr,pos_qtl=pos,ref_qtl=ref,alt_qtl=alt) %>% unique()
        ) %>% 
          inner_join(
            ld.info %>% select(rsnum,chr,pos_ld=pos,ref_ld=ref,alt_ld=alt) %>% unique()
          )
      }
      
    }
    
    cat("\nFinal Set of SNPs for ezQTL analyses",file=logfile,sep="\n",append = T)
    snpall <- snpall %>% filter(!(pos_gwas!=pos_ld|pos_gwas!=pos_qtl|pos_qtl!=pos_ld | ref_gwas!=ref_ld|ref_gwas!=ref_qtl|ref_qtl!=ref_ld|alt_gwas!=alt_ld|alt_gwas!=alt_qtl|alt_qtl!=alt_ld))
    cat(paste0('# number of variants included: ',dim(snpall)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpall %>% filter(str_detect(rsnum,'^rs'))
    cat(paste0('# number of variants with rsnum: ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpalltmp  %>% filter((ref_qtl %in% c('A','T') & alt_qtl %in% c('A','T'))|(ref_qtl %in% c('C','G') & alt_qtl %in% c('C','G')))
    cat(paste0('# number of variants with rsnum and A/T, C/G alleles: ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    
    # output the final dataset ------------------------------------------------
    #orsnum <- gwas %>% filter(rsnum %in% qtl$rsnum,rsnum %in% ld.data$rsnum) %>% pull(rsnum)
    orsnum <- snpall %>% pull(rsnum)
    
    # suggestive snp as the key
    qtl_ref <- qtl %>% filter(rsnum %in% orsnum) %>% arrange(pval_nominal) %>% slice(1)
    gwas_ref <- gwas %>% filter(rsnum %in% orsnum) %>% arrange(pvalue) %>% slice(1)
    
    leadsnp=gwas_ref$rsnum
    
    gwas_ref <- paste0(gwas_ref$rsnum, ' (',gwas_ref$chr,':',gwas_ref$pos,':',gwas_ref$ref,':',gwas_ref$alt,') GWAS P=', gwas_ref$pvalue)
    qtl_ref <- paste0(qtl_ref$rsnum, ' (',qtl_ref$chr,':',qtl_ref$pos,':',qtl_ref$ref,':',qtl_ref$alt,') QTL P=', qtl_ref$pval_nominal,' for ',qtl_ref$gene_id,":",qtl_ref$gene_symbol)
    
    cat('\nConsider the following snp as the reference snp in ezQTL: ',file=logfile,sep="\n",append = T)
    cat(paste0('Best GWAS overlapped snp: ',gwas_ref),file=logfile,sep="\n",append = T)
    cat(paste0('Best QTL overlapped snp: ',qtl_ref),file=logfile,sep="\n",append = T)
    
    ## additional filter for the QTL traits
    # recalculate the LD
    if(!ld_leadsnp){
      leadindex <- ld.info %>% filter(rsnum==leadsnp) %>% pull(Seq)
      #!is.na(leadindex) | 
      if(length(leadindex)!=0){
        ldtmp <- tibble(LD=ld.matrix[,leadindex])
      }else{
        ldtmp <- tibble(LD=rep(NA,dim(ld.matrix)[1]))
      }
      
      ld.data <- bind_cols(
        ld.info,
        ldtmp
      ) %>% filter(chr==leadchr,pos>=leadpos1,pos<=leadpos2)
      
    }
    
    
    
    ## overlap plot 1
    qtl2 <- qtl %>% group_by(rsnum,pos) %>% arrange(pval_nominal) %>% summarise(pvalue=min(pval_nominal),n=n()) %>% ungroup()
    ndifgene <- qtl %>% count(gene_id) %>% count(n) %>% dim() %>% .[[1]]
    
    pcol <- c("FALSE"='#cccccc',"TRUE"='#1a9850')
    p1 <- qtl2 %>% mutate(tmp=rsnum %in% orsnum) %>% arrange(tmp) %>% 
      ggplot(aes(pos/1000000,-log10(pvalue),fill=if_else(rsnum %in% orsnum,"TRUE","FALSE")))+
      geom_point(aes(size=n),pch=21,stroke=0.3,col="black")+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XY",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill=paste0("Overlapped Variants: ",length(orsnum)),x='',y='Minimal QTL P-value (-log10)',size="QTL-n")+
      theme(legend.position = 'top',panel.background = element_blank(),legend.background = element_blank())+
      scale_fill_manual(values = pcol)+
      scale_x_continuous(breaks = pretty_breaks())+
      panel_border(color = 'black')+
      ggrepel::geom_text_repel(data = qtl2 %>% filter(rsnum==leadsnp),aes(label=rsnum))+
      geom_point(data = qtl2 %>% filter(rsnum==leadsnp),pch=2,col="red",size=2,fill=NA,stroke=1)
    
    if(ndifgene > 1 ){
      p1 <- p1 +  scale_size_binned(breaks = pretty_breaks())
    }
    
    p2 <- gwas %>% 
      ggplot(aes(pos/1000000,-log10(pvalue),fill=if_else(rsnum %in% orsnum,"TRUE","FALSE")))+
      geom_point(pch=21,stroke=0.2,col="black",size=2.5)+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XY",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill="Overlapped Variants",x='',y='GWAS P-value (-log10)')+
      theme(legend.position = 'none',panel.background = element_blank())+
      scale_fill_manual(values = pcol)+
      scale_x_continuous(breaks = pretty_breaks())+
      panel_border(color = 'black')+
      ggrepel::geom_text_repel(data = gwas %>% filter(rsnum==leadsnp),aes(label=rsnum) )+
      geom_point(data = gwas %>% filter(rsnum==leadsnp),pch=2,col="red",size=2,stroke=1)
    
    
    p3 <- ld.data %>% mutate(tmp=rsnum %in% orsnum) %>% arrange(tmp) %>% 
      ggplot(aes(pos/1000000,(LD)^2,fill=rsnum %in% orsnum))+
      geom_point(pch=21,stroke=0.2,col="black",size=2.5)+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XY",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill="Overlapped Variants",x=paste0('Chromosome ',leadchr,'\n'),y="LD Matrix (R^2)")+
      theme(legend.position = 'none',panel.background = element_blank())+
      scale_fill_manual(values = pcol)+
      scale_x_continuous(breaks = pretty_breaks())+
      panel_border(color = 'black')+
      ggrepel::geom_text_repel(data = ld.data %>% filter(rsnum==leadsnp),aes(label=rsnum))+
      geom_point(data = ld.data %>% filter(rsnum==leadsnp),pch=2,col="red",size=2,stroke=1)
    
    pall1 <- plot_grid(p1+theme(plot.margin = margin(b = 2)),p2+theme(plot.margin = margin(b = 2)),p3+theme(plot.margin = margin(b = 2)),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(1.3,1,1))
    
    # effect size plot # 
    
    if(is.null(zscore_gene)){
      gene_ids <- qtl %>% arrange(pval_nominal) %>% select(gene_symbol) %>% slice(1) %>% pull(gene_symbol)  
    }else{
      gene_ids <- zscore_gene
    }
    
    qdata <- left_join(
      qtl %>% filter(rsnum %in% orsnum,gene_symbol %in% gene_ids|gene_id %in% gene_ids) %>% mutate(Z=slope/slope_se) %>% mutate(Gene=paste0(gene_id,'/',gene_symbol)) %>% select(rsnum, pos, ref,alt, Z,Gene),
      gwas %>% filter(rsnum %in% orsnum) %>% select(rsnum,pos,zscore)
    ) %>% 
      left_join(
        ld.data %>% select(rsnum,pos,LD)
      ) %>% 
      mutate(label=if_else(rsnum==leadsnp,rsnum,'')) %>% 
      mutate(ambiguous_snp=ifelse((ref %in% c('A','T') & alt %in% c('A','T'))|(ref %in% c('C','G') & alt %in% c('C','G')), "Y","N"))
    
    
    #unique(qdata$Gene)
    
    pall2 <- qdata %>% ggplot(aes(zscore,Z,fill=LD))+
      facet_wrap(~Gene,ncol = 1)+
      geom_abline(slope = 1,col='#cccccc',linetype=2)+
      geom_hline(yintercept = 0,col='#cccccc')+
      geom_vline(xintercept = 0,col='#cccccc')+
      geom_point(pch=21,stroke=0.2,col="black",size=2.5)+
      geom_point(data=qdata %>% filter(ambiguous_snp=="Y"),pch=4,stroke=0.5,col="black",size=1)+
      # theme_ipsum_rc(axis_title_just = 'm',grid = F,axis = F,ticks = T,axis_title_size= 14)+
      labs(fill="LD to GWAS Leadsnp\n",x='GWAS-Zscore',y='QTL-Zscore')+
      theme(legend.position = 'top',legend.key.width = unit(2,'cm'),panel.spacing = unit(0.5, "lines"),panel.background = element_blank(),panel.grid = element_blank(),strip.background = element_blank())+
      scale_fill_gsea(limits=c(-1,1))+
      panel_border(color = 'black')+
      geom_point(data = qdata %>% filter(rsnum==leadsnp),pch=2,col="black",size=2,stroke=0.5)+
      ggrepel::geom_text_repel(aes(label=label))
    
    if(is.null(output_plot_prefix)){
      return(list(pall1,pall2))
    }else{
      ggsave(filename = paste0(output_plot_prefix,"_QC_overlapping.svg"),plot = pall1,width = 12,height = 12)
      ggsave(filename = paste0(output_plot_prefix,"_QC_zscore.svg"),plot = pall2,width = 10,height = 8)
    }
    
    
    ## output final dataset
    qtl <- qtl %>% filter(rsnum %in% orsnum)
    total_traits <- qtl %>% count(gene_id,gene_symbol) %>% dim() %>% .[[1]]
    total_traits <- ceiling(total_traits*0.2)
    total_traits <- if_else(total_traits>50,50,total_traits)
    qtl_rm <- qtl %>% mutate(diff=abs(pos-leadpos),direction=if_else(pos>leadpos,"+","-")) %>%
      filter(diff!=0) %>% 
      group_by(gene_symbol,gene_id,direction) %>%
      arrange(diff) %>% 
      slice(1:total_traits) %>% 
      ungroup() %>% 
      count(gene_id,gene_symbol) %>% 
      filter(n<2*total_traits) 
    
    if(dim(qtl_rm)[1]>0){
      qtl <- qtl %>% filter(!(gene_id %in% qtl_rm$gene_id))
      cat(paste0('# number of QTL traits are removed due to low number of SNPs: ',dim(qtl_rm)[1]),file=logfile,sep="\n",append = T)
    }
    
    qtl %>% write_delim(file=paste0(output_prefix,"_qtl.txt"),delim = '\t',col_names = T)
    
    gwas <- gwas %>% filter(rsnum %in% orsnum)
    gwas %>% write_delim(file=paste0(output_prefix,"_gwas.txt"),delim = '\t',col_names = T)
    
    ld.info <- ld.info %>% filter(rsnum %in% orsnum)
    ld.matrix <- ld.matrix[ld.info$Seq, ld.info$Seq]
    bind_cols(ld.info %>% select(-Seq),as_tibble(ld.matrix)) %>% write_delim(file=paste0(output_prefix,"_ld.gz"),delim = '\t',col_names = FALSE)
  }
  
  # The following function work for qtl and gwas only ------------------
  if(!is.null(gwasfile) & !is.null(qtlfile) & is.null(ldfile)){
    
    gwas <- gwas %>% mutate(chr=as.character(chr))
    qtl <- qtl %>% mutate(chr=as.character(chr))
    
    snpall <- inner_join(
      gwas %>% select(rsnum,chr,pos_gwas=pos,ref_gwas=ref,alt_gwas=alt) %>% unique(),
      qtl %>% select(rsnum,chr,pos_qtl=pos,ref_qtl=ref,alt_qtl=alt) %>% unique()
    ) 
    
    cat("\nOverlapped SNPs",file=logfile,sep="\n",append = T)
    cat(paste0('# number of overlapped snps by rsnum: ',dim(snpall)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpall %>% filter(pos_gwas!=pos_qtl | ref_gwas!=ref_qtl| alt_gwas!=alt_qtl)
    cat(paste0('# number of overlapped snps by rsnum, but have different pos, ref, or alt:  ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    if(dim(snpalltmp)[1]>0){
      pertmp <- dim(snpalltmp)[1]/dim(snpall)[1]
      pertmp2 <- percent_format(accuracy = 0.01)(pertmp)
      cat(paste0('\nWarning: found ',pertmp2,' SNPs with same rsnum, but different information. Please check the input datasets. Make sure input the data have the same genome build and allels information. Check the SNP Matching information for detail (snp_not_match.txt). However, ezQTL still can use the overlapped rsnum as the ID.'),file=logfile,sep="\n",append = T)
      snpalltmp %>% write_delim('snp_not_match.txt',delim = '\t',col_names = T)
      if(pertmp > 0.9){
        cat(paste0('Warning: align QTL based on GWAS information'),file=logfile,sep="\n",append = T)
        snpall <- snpall %>% select(rsnum,chr,pos=pos_gwas,ref=ref_gwas,alt=alt_gwas)
        qtl <- qtl %>% filter(rsnum %in% snpall$rsnum) %>% select(-pos,-ref,-alt) %>% left_join(snpall) %>% select(one_of(qtl_colnames))
        
        snpall <- inner_join(
          gwas %>% select(rsnum,chr,pos_gwas=pos,ref_gwas=ref,alt_gwas=alt) %>% unique(),
          qtl %>% select(rsnum,chr,pos_qtl=pos,ref_qtl=ref,alt_qtl=alt) %>% unique()
        )
      }
      
    }
    
    cat("\nFinal Set of SNPs for ezQTL analyses",file=logfile,sep="\n",append = T)
    snpall <- snpall %>% filter(!(pos_gwas!=pos_qtl | ref_gwas!=ref_qtl| alt_gwas!=alt_qtl))
    cat(paste0('# number of variants included: ',dim(snpall)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpall %>% filter(str_detect(rsnum,'^rs'))
    cat(paste0('# number of variants with rsnum: ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpalltmp  %>% filter((ref_gwas %in% c('A','T') & alt_gwas %in% c('A','T'))|(ref_gwas %in% c('C','G') & alt_gwas %in% c('C','G')))
    cat(paste0('# number of variants with rsnum and A/T, C/G alleles: ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    
    # output the final dataset ------------------------------------------------
    #orsnum <- gwas %>% filter(rsnum %in% qtl$rsnum,rsnum %in% ld.data$rsnum) %>% pull(rsnum)
    orsnum <- snpall %>% pull(rsnum)
    
    # suggestive snp as the key
    qtl_ref <- qtl %>% filter(rsnum %in% orsnum) %>% arrange(pval_nominal) %>% slice(1)
    gwas_ref <- gwas %>% filter(rsnum %in% orsnum) %>% arrange(pvalue) %>% slice(1)
    
    leadsnp=gwas_ref$rsnum
    
    gwas_ref <- paste0(gwas_ref$rsnum, ' (',gwas_ref$chr,':',gwas_ref$pos,':',gwas_ref$ref,':',gwas_ref$alt,') GWAS P=', gwas_ref$pvalue)
    qtl_ref <- paste0(qtl_ref$rsnum, ' (',qtl_ref$chr,':',qtl_ref$pos,':',qtl_ref$ref,':',qtl_ref$alt,') QTL P=', qtl_ref$pval_nominal,' for ',qtl_ref$gene_id,":",qtl_ref$gene_symbol)
    
    cat('\nConsider the following snp as the reference snp in ezQTL: ',file=logfile,sep="\n",append = T)
    cat(paste0('Best GWAS overlapped snp: ',gwas_ref),file=logfile,sep="\n",append = T)
    cat(paste0('Best QTL overlapped snp: ',qtl_ref),file=logfile,sep="\n",append = T)
    
    ## overlap plot 1
    qtl2 <- qtl %>% group_by(rsnum,pos) %>% arrange(pval_nominal) %>% summarise(pvalue=min(pval_nominal),n=n()) %>% ungroup()
    ndifgene <- qtl %>% count(gene_id) %>% count(n) %>% dim() %>% .[[1]]
    
    pcol <- c("FALSE"='#cccccc',"TRUE"='#1a9850')
    p1 <- qtl2 %>% mutate(tmp=rsnum %in% orsnum) %>% arrange(tmp) %>% 
      ggplot(aes(pos/1000000,-log10(pvalue),fill=if_else(rsnum %in% orsnum,"TRUE","FALSE")))+
      geom_point(aes(size=n),pch=21,stroke=0.3,col="black")+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XY",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill=paste0("Overlapped Variants: ",length(orsnum)),x='',y='Minimal QTL P-value (-log10)',size="QTL-n")+
      theme(legend.position = 'top',panel.background = element_blank(),legend.background = element_blank())+
      scale_fill_manual(values = pcol)+
      scale_x_continuous(breaks = pretty_breaks())+
      panel_border(color = 'black')+
      ggrepel::geom_text_repel(data = qtl2 %>% filter(rsnum==leadsnp),aes(label=rsnum))+
      geom_point(data = qtl2 %>% filter(rsnum==leadsnp),pch=2,col="red",size=2,fill=NA,stroke=1)
    
    if(ndifgene > 1 ){
      p1 <- p1 +  scale_size_binned(breaks = pretty_breaks())
    }
    
    p2 <- gwas %>% 
      ggplot(aes(pos/1000000,-log10(pvalue),fill=if_else(rsnum %in% orsnum,"TRUE","FALSE")))+
      geom_point(pch=21,stroke=0.2,col="black",size=2.5)+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XY",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill="Overlapped Variants",x='',y='GWAS P-value (-log10)')+
      theme(legend.position = 'none',panel.background = element_blank())+
      scale_fill_manual(values = pcol)+
      scale_x_continuous(breaks = pretty_breaks())+
      panel_border(color = 'black')+
      ggrepel::geom_text_repel(data = gwas %>% filter(rsnum==leadsnp),aes(label=rsnum) )+
      geom_point(data = gwas %>% filter(rsnum==leadsnp),pch=2,col="red",size=2,stroke=1)
    
    pall1 <- plot_grid(p1+theme(plot.margin = margin(b = 2)),p2+theme(plot.margin = margin(b = 2)),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(1.3,1))
    
    # effect size plot # 
    
    if(is.null(zscore_gene)){
      gene_ids <- qtl %>% arrange(pval_nominal) %>% select(gene_symbol) %>% slice(1) %>% pull(gene_symbol)  
    }else{
      gene_ids <- zscore_gene
    }
    
    qdata <- left_join(
      qtl %>% filter(rsnum %in% orsnum,gene_symbol %in% gene_ids|gene_id %in% gene_ids) %>% mutate(Z=slope/slope_se) %>% mutate(Gene=paste0(gene_id,'/',gene_symbol)) %>% select(rsnum, pos, ref,alt, Z,Gene),
      gwas %>% filter(rsnum %in% orsnum) %>% select(rsnum,pos,zscore)
    ) %>% 
      mutate(label=if_else(rsnum==leadsnp,rsnum,'')) %>% 
      mutate(ambiguous_snp=ifelse((ref %in% c('A','T') & alt %in% c('A','T'))|(ref %in% c('C','G') & alt %in% c('C','G')), "Y","N"))
    
    
    #unique(qdata$Gene)
    
    pall2 <- qdata %>% ggplot(aes(zscore,Z))+
      facet_wrap(~Gene,ncol = 1)+
      geom_abline(slope = 1,col='#cccccc',linetype=2)+
      geom_hline(yintercept = 0,col='#cccccc')+
      geom_vline(xintercept = 0,col='#cccccc')+
      geom_point(pch=21,stroke=0.2,col="black",size=2.5,fill="#cccccc")+
      geom_point(data=qdata %>% filter(ambiguous_snp=="Y"),pch=4,stroke=0.5,col="black",size=1)+
      # theme_ipsum_rc(axis_title_just = 'm',grid = F,axis = F,ticks = T,axis_title_size= 14)+
      labs(x='GWAS-Zscore',y='QTL-Zscore')+
      theme(legend.position = 'top',legend.key.width = unit(2,'cm'),panel.spacing = unit(0.5, "lines"),panel.background = element_blank(),panel.grid = element_blank(),strip.background = element_blank())+
      scale_fill_gsea(limits=c(-1,1))+
      panel_border(color = 'black')+
      geom_point(data = qdata %>% filter(rsnum==leadsnp),pch=2,col="black",size=2,stroke=0.5)+
      ggrepel::geom_text_repel(aes(label=label))
    
    if(is.null(output_plot_prefix)){
      return(list(pall1,pall2))
    }else{
      ggsave(filename = paste0(output_plot_prefix,"_QC_overlapping.svg"),plot = pall1,width = 12,height = 12)
      ggsave(filename = paste0(output_plot_prefix,"_QC_zscore.svg"),plot = pall2,width = 10,height = 8)
    }
    
    
    ## output final dataset
    qtl <- qtl %>% filter(rsnum %in% orsnum)
    total_traits <- qtl %>% count(gene_id,gene_symbol) %>% dim() %>% .[[1]]
    total_traits <- ceiling(total_traits*0.2)
    total_traits <- if_else(total_traits>50,50,total_traits)
    qtl_rm <- qtl %>% mutate(diff=abs(pos-leadpos),direction=if_else(pos>leadpos,"+","-")) %>%
      filter(diff!=0) %>% 
      group_by(gene_symbol,gene_id,direction) %>%
      arrange(diff) %>% 
      slice(1:total_traits) %>% 
      ungroup() %>% 
      count(gene_id,gene_symbol) %>% 
      filter(n<2*total_traits) 
    
    if(dim(qtl_rm)[1]>0){
      qtl <- qtl %>% filter(!(gene_id %in% qtl_rm$gene_id))
      cat(paste0('# number of QTL traits are removed due to low number of SNPs: ',dim(qtl_rm)[1]),file=logfile,sep="\n",append = T)
    }
    
    qtl %>% write_delim(file=paste0(output_prefix,"_qtl.txt"),delim = '\t',col_names = T)
    
    gwas <- gwas %>% filter(rsnum %in% orsnum)
    gwas %>% write_delim(file=paste0(output_prefix,"_gwas.txt"),delim = '\t',col_names = T)
  }
  
  # the following function work for qtl and ld file only ------------------
  if(is.null(gwasfile) & !is.null(qtlfile) & !is.null(ldfile)){
    cat('\nWarning: GWAS file is missing, QTL and LD file are detected. Only Locus LD and Locus Alignment may work.',file=logfile,sep="\n",append = T)
    #gwas <- gwas %>% mutate(chr=as.character(chr))
    qtl <- qtl %>% mutate(chr=as.character(chr))
    ld.info <- ld.info %>% mutate(chr=as.character(chr))
    
    snpall <- inner_join(
      qtl %>% select(rsnum,chr,pos_qtl=pos,ref_qtl=ref,alt_qtl=alt) %>% unique(),
      ld.info %>% select(rsnum,chr,pos_ld=pos,ref_ld=ref,alt_ld=alt) %>% unique()
    )
    cat("\nOverlapped SNPs",file=logfile,sep="\n",append = T)
    cat(paste0('# number of overlapped snps by rsnum: ',dim(snpall)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpall %>% filter(pos_qtl!=pos_ld | ref_qtl!=ref_ld|alt_qtl!=alt_ld)
    cat(paste0('# number of overlapped snps by rsnum, but have different pos, ref, or alt:  ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    if(dim(snpalltmp)[1]>0){
      pertmp <- dim(snpalltmp)[1]/dim(snpall)[1]
      pertmp2 <- percent_format(accuracy = 0.01)(pertmp)
      cat(paste0('\nWarning: found ',pertmp2,' SNPs with same rsnum, but different information. Please check the input datasets. Make sure input the data have the same genome build and allels information. Check the SNP Matching information for detail (snp_not_match.txt). However, ezQTL still can use the overlapped rsnum as the ID.'),file=logfile,sep="\n",append = T)
      snpalltmp %>% write_delim('snp_not_match.txt',delim = '\t',col_names = T)
      if(pertmp > 0.9){
        cat(paste0('Warning: align GWAS and QTL dataset based on LD information'),file=logfile,sep="\n",append = T)
        snpall <- snpall %>% select(rsnum,chr,pos=pos_ld,ref=ref_ld,alt=alt_ld)
        qtl <- qtl %>% filter(rsnum %in% snpall$rsnum) %>% select(-pos,-ref,-alt) %>% left_join(snpall) %>% select(one_of(qtl_colnames))
        
        snpall <- inner_join(
          qtl %>% select(rsnum,chr,pos_qtl=pos,ref_qtl=ref,alt_qtl=alt) %>% unique(),
          ld.info %>% select(rsnum,chr,pos_ld=pos,ref_ld=ref,alt_ld=alt) %>% unique()
        )
        
      }
      
    }
    
    cat("\nFinal Set of SNPs for ezQTL analyses",file=logfile,sep="\n",append = T)
    snpall <- snpall %>% filter(!(pos_qtl!=pos_ld | ref_qtl!=ref_ld|alt_qtl!=alt_ld))
    cat(paste0('# number of variants included: ',dim(snpall)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpall %>% filter(str_detect(rsnum,'^rs'))
    cat(paste0('# number of variants with rsnum: ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpalltmp  %>% filter((ref_qtl %in% c('A','T') & alt_qtl %in% c('A','T'))|(ref_qtl %in% c('C','G') & alt_qtl %in% c('C','G')))
    cat(paste0('# number of variants with rsnum and A/T, C/G alleles: ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    
    # output the final dataset ------------------------------------------------
    #orsnum <- gwas %>% filter(rsnum %in% qtl$rsnum,rsnum %in% ld.data$rsnum) %>% pull(rsnum)
    orsnum <- snpall %>% pull(rsnum)
    
    # suggestive snp as the key
    qtl_ref <- qtl %>% filter(rsnum %in% orsnum) %>% arrange(pval_nominal) %>% slice(1)
    leadsnp=qtl_ref$rsnum
    qtl_ref <- paste0(qtl_ref$rsnum, ' (',qtl_ref$chr,':',qtl_ref$pos,':',qtl_ref$ref,':',qtl_ref$alt,') QTL P=', qtl_ref$pval_nominal,' for ',qtl_ref$gene_id,":",qtl_ref$gene_symbol)
    
    cat('\nConsider the following snp as the reference snp in ezQTL: ',file=logfile,sep="\n",append = T)
    cat(paste0('Best QTL overlapped snp: ',qtl_ref),file=logfile,sep="\n",append = T)
    
    
    # recalculate the LD
    if(!ld_leadsnp){
      leadindex <- ld.info %>% filter(rsnum==leadsnp) %>% pull(Seq)
      #!is.na(leadindex) | 
      if(length(leadindex)!=0){
        ldtmp <- tibble(LD=ld.matrix[,leadindex])
      }else{
        ldtmp <- tibble(LD=rep(NA,dim(ld.matrix)[1]))
      }
      
      ld.data <- bind_cols(
        ld.info,
        ldtmp
      ) %>% filter(chr==leadchr,pos>=leadpos1,pos<=leadpos2)
      
    }
    
    
    ## overlap plot 1
    qtl2 <- qtl %>% group_by(rsnum,pos) %>% arrange(pval_nominal) %>% summarise(pvalue=min(pval_nominal),n=n()) %>% ungroup()
    ndifgene <- qtl %>% count(gene_id) %>% count(n) %>% dim() %>% .[[1]]
    
    pcol <- c("FALSE"='#cccccc',"TRUE"='#1a9850')
    p1 <- qtl2 %>% mutate(tmp=rsnum %in% orsnum) %>% arrange(tmp) %>% 
      ggplot(aes(pos/1000000,-log10(pvalue),fill=if_else(rsnum %in% orsnum,"TRUE","FALSE")))+
      geom_point(aes(size=n),pch=21,stroke=0.3,col="black")+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XY",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill=paste0("Overlapped Variants: ",length(orsnum)),x='',y='Minimal QTL P-value (-log10)',size="QTL-n")+
      theme(legend.position = 'top',panel.background = element_blank(),legend.background = element_blank())+
      scale_fill_manual(values = pcol)+
      scale_x_continuous(breaks = pretty_breaks())+
      panel_border(color = 'black')+
      ggrepel::geom_text_repel(data = qtl2 %>% filter(rsnum==leadsnp),aes(label=rsnum))+
      geom_point(data = qtl2 %>% filter(rsnum==leadsnp),pch=2,col="red",size=2,fill=NA,stroke=1)
    
    if(ndifgene > 1 ){
      p1 <- p1 +  scale_size_binned(breaks = pretty_breaks())
    }
    
    p3 <- ld.data %>% mutate(tmp=rsnum %in% orsnum) %>% arrange(tmp) %>% 
      ggplot(aes(pos/1000000,(LD)^2,fill=rsnum %in% orsnum))+
      geom_point(pch=21,stroke=0.2,col="black",size=2.5)+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XY",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill="Overlapped Variants",x=paste0('Chromosome ',leadchr,'\n'),y="LD Matrix (R^2)")+
      theme(legend.position = 'none',panel.background = element_blank())+
      scale_fill_manual(values = pcol)+
      scale_x_continuous(breaks = pretty_breaks())+
      panel_border(color = 'black')+
      ggrepel::geom_text_repel(data = ld.data %>% filter(rsnum==leadsnp),aes(label=rsnum))+
      geom_point(data = ld.data %>% filter(rsnum==leadsnp),pch=2,col="red",size=2,stroke=1)
    
    pall1 <- plot_grid(p1+theme(plot.margin = margin(b = 2)),p3+theme(plot.margin = margin(b = 2)),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(1.3,1))
    
    if(is.null(output_plot_prefix)){
      return(list(pall1,pall2))
    }else{
      ggsave(filename = paste0(output_plot_prefix,"_QC_overlapping.svg"),plot = pall1,width = 12,height = 10)
      #ggsave(filename = paste0(output_plot_prefix,"_QC_zscore.svg"),plot = pall2,width = 10,height = 8)
    }
    
    
    ## output final dataset
    qtl <- qtl %>% filter(rsnum %in% orsnum)
    qtl_rm <- qtl %>% mutate(diff=abs(pos-leadpos),direction=if_else(pos>leadpos,"+","-")) %>%
      filter(diff!=0) %>% 
      group_by(gene_symbol,gene_id,direction) %>%
      arrange(diff) %>% 
      slice(1:50) %>% 
      ungroup() %>% 
      count(gene_id,gene_symbol) %>% 
      filter(n<100) 
    
    if(dim(qtl_rm)[1]>0){
      qtl <- qtl %>% filter(!(gene_id %in% qtl_rm$gene_id))
      cat(paste0('# number of QTL traits are removed due to low number of SNPs: ',dim(qtl_rm)[1]),file=logfile,sep="\n",append = T)
    }
    
    qtl %>% write_delim(file=paste0(output_prefix,"_qtl.txt"),delim = '\t',col_names = T)
    
    ld.info <- ld.info %>% filter(rsnum %in% orsnum)
    ld.matrix <- ld.matrix[ld.info$Seq, ld.info$Seq]
    bind_cols(ld.info %>% select(-Seq),as_tibble(ld.matrix)) %>% write_delim(file=paste0(output_prefix,"_ld.gz"),delim = '\t',col_names = FALSE)
    
  }
  
  # the following function work for gwas and ld only ------------------
  if(!is.null(gwasfile) & is.null(qtlfile) & !is.null(ldfile)){
    cat('\nWarning: QTL file is missing, GWAS and LD file are detected. Only Locus LD and Locus Alignment may work.',file=logfile,sep="\n",append = T)
    
    gwas <- gwas %>% mutate(chr=as.character(chr))
    ld.info <- ld.info %>% mutate(chr=as.character(chr))
    
    snpall <- inner_join(
      gwas %>% select(rsnum,chr,pos_gwas=pos,ref_gwas=ref,alt_gwas=alt) %>% unique(),
      ld.info %>% select(rsnum,chr,pos_ld=pos,ref_ld=ref,alt_ld=alt) %>% unique()
    )
    
    cat("\nOverlapped SNPs",file=logfile,sep="\n",append = T)
    cat(paste0('# number of overlapped snps by rsnum: ',dim(snpall)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpall %>% filter(pos_gwas!=pos_ld| ref_gwas!=ref_ld|alt_gwas!=alt_ld)
    cat(paste0('# number of overlapped snps by rsnum, but have different pos, ref, or alt:  ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    if(dim(snpalltmp)[1]>0){
      pertmp <- dim(snpalltmp)[1]/dim(snpall)[1]
      pertmp2 <- percent_format(accuracy = 0.01)(pertmp)
      cat(paste0('\nWarning: found ',pertmp2,' SNPs with same rsnum, but different information. Please check the input datasets. Make sure input the data have the same genome build and allels information. Check the SNP Matching information for detail (snp_not_match.txt). However, ezQTL still can use the overlapped rsnum as the ID for colocalization analyses.'),file=logfile,sep="\n",append = T)
      snpalltmp %>% write_delim('snp_not_match.txt',delim = '\t',col_names = T)
      if(pertmp > 0.9){
        cat(paste0('Warning: align GWAS and QTL dataset based on LD information'),file=logfile,sep="\n",append = T)
        snpall <- snpall %>% select(rsnum,chr,pos=pos_ld,ref=ref_ld,alt=alt_ld)
        gwas <- gwas %>% filter(rsnum %in% snpall$rsnum) %>% select(-pos,-ref,-alt) %>% left_join(snpall) %>% select(one_of(gwas_colnames))
        snpall <- inner_join(
          gwas %>% select(rsnum,chr,pos_gwas=pos,ref_gwas=ref,alt_gwas=alt) %>% unique(),
          ld.info %>% select(rsnum,chr,pos_ld=pos,ref_ld=ref,alt_ld=alt) %>% unique()
        ) 
      }
    }
    
    cat("\nFinal Set of SNPs for ezQTL analyses",file=logfile,sep="\n",append = T)
    snpall <- snpall %>% filter(!(pos_gwas!=pos_ld| ref_gwas!=ref_ld|alt_gwas!=alt_ld))
    cat(paste0('# number of variants included: ',dim(snpall)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpall %>% filter(str_detect(rsnum,'^rs'))
    cat(paste0('# number of variants with rsnum: ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    snpalltmp <- snpalltmp  %>% filter((ref_gwas %in% c('A','T') & alt_gwas %in% c('A','T'))|(ref_gwas %in% c('C','G') & alt_gwas %in% c('C','G')))
    cat(paste0('# number of variants with rsnum and A/T, C/G alleles: ',dim(snpalltmp)[1]),file=logfile,sep="\n",append = T)
    
    # output the final dataset ------------------------------------------------
    #orsnum <- gwas %>% filter(rsnum %in% qtl$rsnum,rsnum %in% ld.data$rsnum) %>% pull(rsnum)
    orsnum <- snpall %>% pull(rsnum)
    
    # suggestive snp as the key
    gwas_ref <- gwas %>% filter(rsnum %in% orsnum) %>% arrange(pvalue) %>% slice(1)
    leadsnp=gwas_ref$rsnum
    gwas_ref <- paste0(gwas_ref$rsnum, ' (',gwas_ref$chr,':',gwas_ref$pos,':',gwas_ref$ref,':',gwas_ref$alt,') GWAS P=', gwas_ref$pvalue)
    cat('\nConsider the following snp as the reference snp in ezQTL: ',file=logfile,sep="\n",append = T)
    cat(paste0('Best GWAS overlapped snp: ',gwas_ref),file=logfile,sep="\n",append = T)
    
    ## additional filter for the QTL traits
    # recalculate the LD
    if(!ld_leadsnp){
      leadindex <- ld.info %>% filter(rsnum==leadsnp) %>% pull(Seq)
      #!is.na(leadindex) | 
      if(length(leadindex)!=0){
        ldtmp <- tibble(LD=ld.matrix[,leadindex])
      }else{
        ldtmp <- tibble(LD=rep(NA,dim(ld.matrix)[1]))
      }
      
      ld.data <- bind_cols(
        ld.info,
        ldtmp
      ) %>% filter(chr==leadchr,pos>=leadpos1,pos<=leadpos2)
      
    }
    
    ## overlap plot 1
    pcol <- c("FALSE"='#cccccc',"TRUE"='#1a9850')
    
    p2 <- gwas %>% 
      ggplot(aes(pos/1000000,-log10(pvalue),fill=if_else(rsnum %in% orsnum,"TRUE","FALSE")))+
      geom_point(pch=21,stroke=0.2,col="black",size=2.5)+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XY",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill="Overlapped Variants",x='',y='GWAS P-value (-log10)')+
      theme(legend.position = 'none',panel.background = element_blank())+
      scale_fill_manual(values = pcol)+
      scale_x_continuous(breaks = pretty_breaks())+
      panel_border(color = 'black')+
      ggrepel::geom_text_repel(data = gwas %>% filter(rsnum==leadsnp),aes(label=rsnum) )+
      geom_point(data = gwas %>% filter(rsnum==leadsnp),pch=2,col="red",size=2,stroke=1)
    
    
    p3 <- ld.data %>% mutate(tmp=rsnum %in% orsnum) %>% arrange(tmp) %>% 
      ggplot(aes(pos/1000000,(LD)^2,fill=rsnum %in% orsnum))+
      geom_point(pch=21,stroke=0.2,col="black",size=2.5)+
      theme_ipsum_rc(axis_title_just = 'm',grid = "XY",axis = F,ticks = T,axis_title_size= 14)+
      labs(fill="Overlapped Variants",x=paste0('Chromosome ',leadchr,'\n'),y="LD Matrix (R^2)")+
      theme(legend.position = 'none',panel.background = element_blank())+
      scale_fill_manual(values = pcol)+
      scale_x_continuous(breaks = pretty_breaks())+
      panel_border(color = 'black')+
      ggrepel::geom_text_repel(data = ld.data %>% filter(rsnum==leadsnp),aes(label=rsnum))+
      geom_point(data = ld.data %>% filter(rsnum==leadsnp),pch=2,col="red",size=2,stroke=1)
    
    pall1 <- plot_grid(p2+theme(plot.margin = margin(b = 2)),p3+theme(plot.margin = margin(b = 2)),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(1.3,1))
    
    
    if(is.null(output_plot_prefix)){
      return(list(pall1,pall2))
    }else{
      ggsave(filename = paste0(output_plot_prefix,"_QC_overlapping.svg"),plot = pall1,width = 12,height = 12)
    }
    
    ## output final dataset
    gwas <- gwas %>% filter(rsnum %in% orsnum)
    gwas %>% write_delim(file=paste0(output_prefix,"_gwas.txt"),delim = '\t',col_names = T)
    
    ld.info <- ld.info %>% filter(rsnum %in% orsnum)
    ld.matrix <- ld.matrix[ld.info$Seq, ld.info$Seq]
    bind_cols(ld.info %>% select(-Seq),as_tibble(ld.matrix)) %>% write_delim(file=paste0(output_prefix,"_ld.gz"),delim = '\t',col_names = FALSE)
  }
  
  # the following function work for ld only ------------------
  if(is.null(gwasfile) & is.null(qtlfile) & !is.null(ldfile)){
    cat('\nWarning: QTL and GWAS files are missing, only LD file are detected. Only Locus LD may work.',file=logfile,sep="\n",append = T)
    bind_cols(ld.info %>% select(-Seq),as_tibble(ld.matrix)) %>% write_delim(file=paste0(output_prefix,"_ld.gz"),delim = '\t',col_names = FALSE)
  }
  
  # the following function work for gwas only ------------------
  if(!is.null(gwasfile) & is.null(qtlfile) & is.null(ldfile)){
    cat('\nWarning: QTL and LD files are missing, only GWAS file are detected. Only Locus Alignment may work.',file=logfile,sep="\n",append = T)
    gwas %>% write_delim(file=paste0(output_prefix,"_gwas.txt"),delim = '\t',col_names = T)
  }
  
  # the following function work for qtl only ------------------
  if(is.null(gwasfile) & !is.null(qtlfile) & is.null(ldfile)){
    cat('\nWarning: GWAS and LD files are missing, only QTL file are detected. Only Locus Alignment may work.',file=logfile,sep="\n",append = T)
    qtl %>% write_delim(file=paste0(output_prefix,"_qtl.txt"),delim = '\t',col_names = T)
  }
  
}





hycoloc_barplot <- function(hydata,output_plot=NULL,plot_width=NULL,plot_height=NULL){
  
  hydata <- hydata %>% mutate(posterior_prob=as.numeric(posterior_prob),candidate_snp=as.character(candidate_snp))
  
  ngene <- nrow(hydata)
  
  if(ngene>50){
    p <- hydata %>% 
      mutate(posterior_prob=if_else(is.na(posterior_prob),0,posterior_prob)) %>% 
      mutate(label=if_else(is.na(candidate_snp)|posterior_prob<0.5,NA_character_,paste0(gene_symbol,"/",candidate_snp))) %>% 
      #mutate(label=paste0(gene_symbol,"/",candidate_snp)) %>% 
      mutate(gene_symbol=fct_reorder(gene_symbol,posterior_prob)) %>% 
      ggplot(aes(gene_symbol,posterior_prob,fill=posterior_prob))+
      geom_hline(yintercept = 0.5,col='#2CA02CFF',size=0.5)+
      geom_col(width = 0.5,col='black',size=0.25)+
      scale_fill_material("deep-purple")+
      geom_text_repel(aes(label=label),hjust=1,vjust=0,nudge_y = 0.1,segment.color="black",segment.curvature = -0.1,segment.ncp = 3,segment.angle = 20)+
      scale_y_continuous(limits = c(0,1),breaks = pretty_breaks(),expand = c(0,0))+
      labs(x='QTL Trait',y='Posterior Probability\n')+
      theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Yy',axis = "XY",axis_col = 'black',base_family='Roboto Condensed')+
      theme(axis.text.x = element_blank(),legend.position = 'none',panel.grid = element_line(linetype = 5))
    
  }else{
    
    p <- hydata %>% 
      mutate(posterior_prob=if_else(is.na(posterior_prob),0,posterior_prob)) %>% 
      mutate(label=if_else(is.na(candidate_snp)|posterior_prob<0.5,NA_character_,candidate_snp)) %>% 
      #mutate(label=paste0(gene_symbol,"/",candidate_snp)) %>% 
      mutate(gene_symbol=fct_reorder(gene_symbol,posterior_prob)) %>% 
      ggplot(aes(gene_symbol,posterior_prob,fill=posterior_prob))+
      geom_hline(yintercept = 0.5,col='#2CA02CFF',size=0.5)+
      geom_col(width = 0.5,col='black',size=0.25)+
      scale_fill_material("deep-purple")+
      geom_text_repel(aes(label=label),hjust=1,vjust=0,nudge_y = 0.1,segment.color="black",segment.curvature = -0.1,segment.ncp = 3,segment.angle = 20)+
      scale_y_continuous(limits = c(0,1),breaks = pretty_breaks(),expand = c(0,0))+
      labs(x='QTL Trait',y='Posterior Probability\n')+
      theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Yy',axis = "XY",axis_col = 'black',base_family='Roboto Condensed')+
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = 'none',panel.grid = element_line(linetype = 5))
    
  }
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 4+ngene*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
}




hycoloc_boxplot <- function(hydata_score,output_plot=NULL,plot_width=NULL,plot_height=NULL){
  if(dim(hydata_score)[1]==0){
    labeltext="No colocalization traits found by HyPrColoc. No SNP score will be calculated"
    p <- ggplot()+geom_text(aes(x=1,y=1,label=labeltext),family="Roboto Condensed",size=6)+theme_nothing()
    ngene <- 2
  }else{
    
    ngene <- hydata_score %>% select(gene_symbol) %>% unique() %>% nrow()
    
    if(ngene<1){
      errinfo <- "ERROR: No coloclaization trait found between GWAS and QTL"
      return(errinfo)
    }else {
      p <- hydata_score %>% 
        mutate(label=if_else(snpscore>0.05,rsnum,NA_character_)) %>% 
        ggplot(aes(snpscore,gene_symbol))+
        geom_vline(xintercept = 0.05,col='#2CA02CFF',size=0.5)+
        geom_boxplot(outlier.size = 2,outlier.shape = 21,outlier.fill = "#D62728FF",size=0.3)+
        geom_text_repel(aes(label=label),hjust=1,vjust=0,nudge_y = 0.1,segment.color="black",segment.curvature = -0.1,segment.ncp = 3,segment.angle = 20)+
        scale_x_continuous(breaks = pretty_breaks())+
        labs(x='SNP Score',y='QTL Trait')+
        theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'XY',axis_col = 'black',base_family='Roboto Condensed')+
        theme(legend.position = 'none',panel.grid.major.x = element_line(linetype = 5))+
        panel_border()
    }
  }
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 10
    yleng <- 2+ngene
    yleng <- if_else(yleng>10,10,yleng)
    
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
}




ecaviar_visualize <- function(ecdata,output_plot_prefix=NULL,plot_width=NULL,plot_height=NULL){
  
  ecdata2 <- ecdata %>% 
    mutate(CLPP=as.numeric(CLPP),CLPP2=as.numeric(CLPP2)) %>% 
    select(gene_symbol,rsnum,Leadsnp,leadsnp_included,CLPP,CLPP2) %>% 
    pivot_longer(cols = c(CLPP,CLPP2),names_to = 'Group',values_to = 'CLPP') %>% 
    filter(!is.na(CLPP)) %>% 
    mutate(Group=if_else(Group=="CLPP","window=100kb","window=50 SNPs")) 
  
  
  genelevel <- ecdata2 %>% group_by(gene_symbol) %>% arrange(desc(CLPP)) %>% slice(1) %>% ungroup() %>% arrange(CLPP) %>% pull(gene_symbol)
  nudge_y <- 0.03*max(ecdata2$CLPP)
  ngene <- length(genelevel)
  
  p1 <- ecdata2 %>% 
    group_by(gene_symbol,Group) %>% 
    arrange(desc(CLPP)) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(label=if_else(is.na(rsnum)|CLPP<0.01,NA_character_,rsnum)) %>% 
    mutate(gene_symbol=factor(gene_symbol,levels = genelevel)) %>% 
    ggplot(aes(gene_symbol,CLPP,fill=CLPP))+
    facet_wrap(~Group,ncol = 1)+
    geom_hline(yintercept = 0.01,col='#2CA02CFF',size=0.5)+
    geom_col(width = 0.5,col='black',size=0.25)+
    scale_fill_material("deep-purple")+
    geom_text_repel(aes(label=label),hjust=1,vjust=0,nudge_y = nudge_y,segment.color="black",segment.curvature = -0.1,segment.ncp = 3,segment.angle = 20)+
    scale_y_continuous(breaks = pretty_breaks(),expand = c(0,0))+
    labs(x='QTL Trait',y='Colocalization Posterior Probability (CLPP)\n')+
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Yy',axis = "XY",axis_col = 'black',base_family='Roboto Condensed')+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = 'none',panel.grid = element_line(linetype = 5))
  
  p2 <- ecdata2 %>% 
    mutate(gene_symbol=factor(gene_symbol,levels = genelevel)) %>% 
    mutate(label=if_else(CLPP>0.01,rsnum,NA_character_)) %>% 
    ggplot(aes(CLPP,gene_symbol))+
    geom_vline(xintercept = 0.01,col='#2CA02CFF',size=0.5)+
    facet_wrap(~Group,ncol = 1)+
    geom_boxplot(outlier.size = 2,outlier.shape = 21,outlier.fill = "#D62728FF",size=0.3)+
    geom_text_repel(aes(label=label),hjust=1,vjust=0,nudge_y = 0.1,segment.color="black",segment.curvature = -0.1,segment.ncp = 3,segment.angle = 20,force = 5)+
    scale_x_continuous(breaks = pretty_breaks())+
    labs(x='Colocalization Posterior Probability (CLPP)',y='QTL Trait')+
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'XY',axis_col = 'black',base_family='Roboto Condensed')+
    theme(legend.position = 'none',panel.grid.major.x = element_line(linetype = 5))+
    panel_border()
  
  if(length(genelevel)>50){
    p1 <- p1 + theme(axis.text.x = element_blank())
    p2 <- p2 + theme(axis.text.y = element_blank())
  }
  
  if(is.null(output_plot_prefix)){
    return(list(p1,p2))
  }else{
    output_plot <- paste0(output_plot_prefix,"_barplot.svg")
    xleng <- 4+ngene*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- 12
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p1,width = plot_width,height = plot_height)
    
    output_plot <- paste0(output_plot_prefix,"_boxplot.svg")
    xleng <- 10
    yleng <- 2*(2+ngene)
    yleng <- if_else(yleng>16,16,yleng)
    
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p2,width = plot_width,height = plot_height)
  }
  
  
}




coloc_visualize <- function(hydata,ecdata,output_plot=NULL,plot_width=NULL,plot_height=NULL){
  
  hydata <- hydata %>% mutate(posterior_prob=as.numeric(posterior_prob),candidate_snp=as.character(candidate_snp))
  hydata <- hydata %>%  mutate(posterior_prob=if_else(is.na(posterior_prob),0,posterior_prob)) 
  
  ecdata2 <- ecdata %>% 
    mutate(CLPP=as.numeric(CLPP),CLPP2=as.numeric(CLPP2)) %>% 
    mutate(CLPP=if_else(is.na(CLPP),0,CLPP),CLPP2=if_else(is.na(CLPP2),0,CLPP2)) %>% 
    select(gene_id,gene_symbol,rsnum,Leadsnp,leadsnp_included,CLPP,CLPP2) %>% 
    pivot_longer(cols = c(CLPP,CLPP2),names_to = 'Group',values_to = 'CLPP') %>% 
    filter(!is.na(CLPP)) %>% 
    mutate(Group=if_else(Group=="CLPP","window=100kb","window=50 SNPs")) 
  
  cdata <- hydata %>% 
    select(gene_id,gene_symbol,posterior_prob) %>% 
    left_join(
      ecdata2 %>% group_by(gene_id,gene_symbol) %>% summarise(CLPP=max(CLPP)) %>% arrange(desc(CLPP))
    ) %>% 
    arrange(posterior_prob,CLPP)
  
  if(all(is.na(cdata$gene_symbol))){
    hydata <- hydata %>% mutate(gene_symbol=gene_id)
    ecdata <- ecdata %>% mutate(gene_symbol=gene_id)
    ecdata2 <- ecdata2 %>% mutate(gene_symbol=gene_id)
    cdata <- cdata %>% mutate(gene_symbol=gene_id)
  }else{
    if(any(is.na(cdata$gene_symbol))){
      hydata <- hydata %>% mutate(gene_symbol=paste0(gene_id,'|',gene_symbol))
      ecdata <- ecdata %>% mutate(gene_symbol=paste0(gene_id,'|',gene_symbol))
      ecdata2 <- ecdata2 %>% mutate(gene_symbol=paste0(gene_id,'|',gene_symbol))
      cdata <- cdata %>% mutate(gene_symbol=paste0(gene_id,'|',gene_symbol))
    }
    
  }
  
  genelevel <- cdata %>% pull(gene_symbol)
  
  p1 <- hydata %>% 
    mutate(posterior_prob=if_else(is.na(posterior_prob),0,posterior_prob)) %>% 
    mutate(label=if_else(is.na(candidate_snp)|posterior_prob<0.5,NA_character_,candidate_snp)) %>% 
    #mutate(label=paste0(gene_symbol,"/",candidate_snp)) %>% 
    mutate(gene_symbol=fct_reorder(gene_symbol,posterior_prob)) %>% 
    ggplot(aes(gene_symbol,posterior_prob,fill=posterior_prob))+
    geom_hline(yintercept = 0.5,col='#2CA02CFF',size=0.5)+
    geom_col(width = 0.5,col='black',size=0.25)+
    scale_fill_material("deep-purple")+
    geom_text_repel(aes(label=label),hjust=1,vjust=0,nudge_y = 0.1,segment.color="black",segment.curvature = -0.1,segment.ncp = 3,segment.angle = 20)+
    scale_y_continuous(limits = c(0,1),breaks = pretty_breaks(),expand = c(0,0))+
    labs(x='QTL Trait',y='HyPrColoc: Posterior Probability\n')+
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Yy',axis = "XY",axis_col = 'black',base_family='Roboto Condensed')+
    theme(axis.text.x = element_blank(),axis.title.x = element_blank(),legend.position = 'none',panel.grid = element_line(linetype = 5))
  
  
  nudge_y <- 0.03*max(ecdata2$CLPP)
  ngene <- length(genelevel)
  
  p2 <- ecdata2 %>% 
    group_by(gene_symbol) %>% 
    arrange(desc(CLPP)) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(label=if_else(is.na(rsnum)|CLPP<0.01,NA_character_,rsnum)) %>% 
    mutate(gene_symbol=factor(gene_symbol,levels = genelevel)) %>% 
    ggplot(aes(gene_symbol,CLPP,fill=CLPP))+
    geom_hline(yintercept = 0.01,col='#2CA02CFF',size=0.5)+
    geom_col(width = 0.5,col='black',size=0.25)+
    scale_fill_material("deep-purple")+
    geom_text_repel(aes(label=label),hjust=1,vjust=0,nudge_y = nudge_y,segment.color="black",segment.curvature = -0.1,segment.ncp = 3,segment.angle = 20)+
    scale_y_continuous(breaks = pretty_breaks(),expand = c(0,0))+
    labs(x='QTL Trait',y='eCAVIAR: Colocalization Posterior Probability\n')+
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Yy',axis = "XY",axis_col = 'black',base_family='Roboto Condensed')+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = 'none',panel.grid = element_line(linetype = 5))
  
  p <- plot_grid(p1+theme(plot.margin = margin(b = 2)),p2+theme(plot.margin = margin(t = 2)),align = 'v',axis = 'lr',ncol = 1)
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 4+ngene*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- 12
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
}








locus_quantification_cor <- function(qdata,qtldata,output_plot=NULL,plot_width=NULL,plot_height=NULL){
  
  if(any(is.na(unique(qtldata$gene_symbol)))){
    qtldata <- qtldata %>% mutate(gene_symbol=paste0(gene_id,'\n',gene_symbol))
  }
  
  
  id2symbol <- qtldata %>% select(gene_id,gene_symbol) %>% unique()
  qdata <- qdata %>% filter(gene_id %in% id2symbol$gene_id)
  
  mqdata <- qdata %>% select(-chr,-start,-end,-gene_id) %>% as.matrix()
  rownames(mqdata) <- qdata %>% select(gene_id) %>% left_join(id2symbol) %>% pull(gene_symbol)
  mqdata <- t(mqdata)
  # Pearson correlation
  pear <- corrr::correlate(mqdata, method = "pearson", quiet= TRUE) %>%
    corrr::shave(upper = TRUE) %>%
    corrr::stretch() %>%
    stats::na.omit() %>%
    dplyr::rename(pearson_r = "r")
  # Spearman correlation
  spear <- corrr::correlate(mqdata, method = "spearman", quiet= TRUE) %>%
    corrr::shave(upper = TRUE) %>%
    corrr::stretch() %>%
    stats::na.omit() %>%
    dplyr::rename(spearman_r = "r")
  
  df <- dplyr::full_join(pear, spear, by = c("x", "y")) %>%
    ggasym::asymmetrise(x, y)
  
  
  p <- ggplot(df) +
    ggasym::geom_asymmat(aes(x = x, y = y,
                     fill_tl = pearson_r, fill_br = spearman_r)) +
    ggasym::scale_fill_tl_distiller(type = "div", palette = "RdYlBu",
                            na.value = "grey90",
                            guide = guide_colourbar(direction = "vertical",barheight = unit(6,'cm'),
                                                    order = 1,
                                                    title.position = "top")) +
    ggasym::scale_fill_br_distiller(type = "div", palette = "RdYlBu",
                            na.value = "grey90",
                            guide = guide_colourbar(direction = "vertical",barheight = unit(6,'cm'),
                                                    order = 2,
                                                    title.position = "top")) +
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,axis_text_size = 10,grid = FALSE,axis = FALSE,ticks = FALSE)+
    theme(panel.background = element_rect(fill = "grey25"),axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
          panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_equal() +
    panel_border(size = 0.5,colour = 'black')+
    labs(x = "Spearman Correlations", y = "Pearson Correlations",
         title = "Comparing Pearson and Spearman Correlations",
         fill_tl = "Pearson r", fill_br = "Spearman r")
  
  ngene <- dim(mqdata)[2]
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 2+ngene*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- xleng
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
}







locus_quantification_dis<- function(qdata,qtldata,genesets=NULL,output_plot=NULL,plot_width=NULL,plot_height=NULL){
  
  if(any(is.na(unique(qtldata$gene_symbol)))){
    qtldata <- qtldata %>% mutate(gene_symbol=paste0(gene_id,'\n',gene_symbol))
  }
  
  
  id2symbol <- qtldata %>% select(gene_id,gene_symbol) %>% unique()
  
  if(is.null(genesets)){
    genesets <- qtldata %>% group_by(gene_symbol) %>% arrange(pval_nominal) %>% slice(1) %>% ungroup() %>% arrange(pval_nominal) %>% slice(1:15) %>% pull(gene_symbol)
  }
  genesets2 <- id2symbol %>% filter(gene_symbol %in% genesets) %>% pull(gene_id)
  qdata2 <- qdata %>% 
    filter(gene_id %in% genesets2) %>% 
    left_join(id2symbol) %>% 
    pivot_longer(cols = -c(chr,start,end,gene_id,gene_symbol),names_to = 'sample',values_to = 'expression') %>% 
    mutate(expression=log2(expression+1e-3))
  
  genelevel <- qdata2 %>% group_by(gene_symbol) %>% summarise(mexp=median(expression,na.rm = TRUE)) %>% arrange(mexp) %>% pull(gene_symbol)
  
  p <- qdata2 %>% 
    mutate(gene_symbol=factor(gene_symbol,levels = genelevel)) %>% 
    ggplot(aes(x=expression,y=gene_symbol,fill=stat(x)))+
    #stat_density_ridges(quantile_lines = TRUE, quantiles = 2,fill='#cccccc')+
    ggridges::geom_density_ridges_gradient() + 
    scale_fill_viridis_c(option = 'C')+
    scale_x_continuous(expand = c(0, 0),breaks = pretty_breaks()) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
    labs(x="Quantification (log2)",y='QTL Trait')+
    theme(legend.position = 'none')
  
  ngene <- length(genesets)
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 10
    yleng <- 4+ngene*0.4
    yleng <- if_else(yleng>15,15,yleng)
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
}




# Locus LD ----------------------------------------------------------------

IntRegionalPlot <- function(chr=NULL, left=NULL, right=NULL, association_file=NULL, trait=NULL, LDfile, genome_build="GRCh37",gtf_tabix_file, Distance = 50000,Distance_max = 2000000, output_file = NULL,
                            slide_length = -1, threadN = 1, ldstatistics = "rsquare", leadsnp = NULL,leadsnp_pos = NULL, threshold = NULL, 
                            link2gene = NULL, triangleLD = TRUE, link2LD = NULL, leadsnpLD = TRUE, label_gene_name = FALSE, 
                            colour02 = "gray", colour04 = "cyan", colour06 = "green", colour08 = "yellow", 
                            colour10 = "red", leadsnp_shape = 23, leadsnp_colour = "black", leadsnp_fill = "purple", 
                            leadsnp_size = 1.5, marker2highlight = NULL, marker2label = NULL, marker2label_angle = 60, 
                            marker2label_size = 1) {
  
  
  if(is.null(association_file) && is.null(leadsnp) && is.null(leadsnp_pos)){
    stop("You have to input the SNP in the locus informaiton for the LD visualization.")
  }
  
  ld.matrix <- read_delim(LDfile,delim = '\t',col_names = F,col_types = cols('X1'='c')) %>% rename(chr=X1,pos=X2,rsnum=X3,ref=X4,alt=X5)
  ld.info <- ld.matrix %>% dplyr::select(chr,pos,rsnum,ref,alt) %>% mutate(Seq=seq_along(chr))
  ld.matrix <- ld.matrix %>% dplyr::select(-c(chr,pos,rsnum,ref,alt)) %>% as.matrix
  rownames(ld.matrix) <- ld.info$Seq
  colnames(ld.matrix) <- ld.info$Seq
  
  if(is.null(association_file) && !(leadsnp %in% ld.info$rsnum) && !(leadsnp_pos %in% ld.info$pos)){
    stop("The input SNP informaiton does not existed in LD file")
  }
  
  if(is.null(association_file)){
    
    threshold=9
    if(!is.null(leadsnp)) {association <- ld.info %>% dplyr::mutate(p=if_else((rsnum==leadsnp),1e-10,0.05)) %>% dplyr::select(Marker=rsnum,Locus=chr,Site=pos,p)}
    if(!is.null(leadsnp_pos)) {association <- ld.info %>% dplyr::mutate(p=if_else((pos==leadsnp_pos),1e-10,0.05)) %>% dplyr::select(Marker=rsnum,Locus=chr,Site=pos,p)}
    
  }else{
    association <- read_delim(association_file,delim = '\t',col_names = T,n_max = 5)
    if("variant_id" %in% colnames(association)){
      association <- read_delim(association_file,delim = '\t',col_names = T,col_types = cols('variant_id'='c'))  
    }else{
      association <- read_delim(association_file,delim = '\t',col_names = T)
    }
    
    if('pvalue' %in% colnames(association)){  association <- association %>% select(Marker=rsnum,Locus=chr,Site=pos,p=pvalue) }
    if('pval_nominal' %in% colnames(association)){  
      if(is.null(trait)){
        trait <- association %>% arrange(pval_nominal) %>% slice(1) %>% pull(gene_symbol)
      }
      association <- association %>% filter(gene_symbol == trait) %>% select(Marker=rsnum,Locus=chr,Site=pos,p=pval_nominal) 
    }
  }
  
  
  if(any(is.null(chr),is.null(left),is.null(right))){
    tmplocus <- association %>% arrange(p) %>% slice(1) 
    chr <- tmplocus$Locus
    left <- tmplocus$Site - Distance
    right <- tmplocus$Site + Distance
  }
  # limmit to 100k
  
  if((right - left)  > Distance_max){
    middle <- (right-left)/2
    left <- middle - Distance_max/2
    right <- middle + Distance_max/2
  }
  
  ## redefine the region if the association region is mall
  if((max(association$Site)-min(association$Site)) < 2*Distance){
    left <- min(association$Site)
    right <- max(association$Site)
  }
  
  chromosome_association <- association[association$Locus == chr, ]
  transcript_min <- left
  transcript_max <- right
  transcript_association <- chromosome_association[chromosome_association$Site >= transcript_min & chromosome_association$Site <= transcript_max, ]
  transcript_association <- transcript_association[order(transcript_association$Site),]
  
  if (dim(transcript_association)[1] < 2) {
    stop("Less than 2 markers, can not compute LD")
  } else {
    
    ## tabix gtf file
    #if(genome_build == "GRCh37") { gtf_tabix_file <- paste0(gtf_tabix_folder,'/gencode.v19.annotation.gtf.gz')}
    #if(genome_build == "GRCh38") { gtf_tabix_file <- paste0(gtf_tabix_folder,'/gencode.v37.annotation.gtf.gz')}
    regionfile=paste0(chr,":",left,"-",right,'.gtf_tmp')
    #cmd_ztw=paste0('tabix ',gtf_tabix_file,' ',chr,":",left,"-",right,' > ',regionfile)
    cmd_ztw=paste0('tabix ',gtf_tabix_file,' ',chr,":",left,"-",right,' -D > ',regionfile)
    system(cmd_ztw)
    gtf <- read_delim(regionfile,delim = '\t',col_names = FALSE)
    colnames(gtf) <- paste0('V',1:9)
    unlink(x = regionfile,force = TRUE)
    
    R2 <- Site <- Site2 <- V4 <- V5 <- V9 <- group <- p <- NULL
    x <- xend <- y <- yend <- aggregate <- NULL
    pvalue_range <- pretty(-log10(transcript_association$p))
    fold <- ((transcript_max - transcript_min) * 2/3)/max(pvalue_range)
    n_pvalue_range <- length(pvalue_range)
    gene_list <- gtf[gtf$V1 == chr & grepl("exon",gtf$V3,ignore.case = TRUE),]
    gene_list$V9 <- gsub("\"|;", "", gene_list$V9)
    gene_list$V9 <- sub(".*gene_name (\\S+) .+", "\\1", gene_list$V9)
    gene_list_start <- aggregate(V4 ~ V9, data = gene_list, FUN = min)
    gene_list_end <- aggregate(V5 ~ V9, data = gene_list, FUN = max)
    gene_list$V4 <- NULL
    gene_list$V5 <- NULL
    gene_list <- merge(gene_list, gene_list_start, by = "V9")
    gene_list <- merge(gene_list, gene_list_end, by = "V9")
    gene_list <- gene_list[!duplicated(gene_list$V9), ]
    gene_list <- gene_list[gene_list$V4 >= transcript_min & gene_list$V5 <= 
                             transcript_max, ]
    gene_for <- gene_list[gene_list$V7 == "+", ]
    gene_rev <- gene_list[gene_list$V7 == "-", ]
    if (nrow(gene_for) >= 1) {
      # plot forward strand gene
      gene_for_seg <- list(geom_segment(data = gene_for, aes(x = V4, y = -(transcript_max -  transcript_min)/30, xend = V5, yend = -(transcript_max - transcript_min)/30), arrow = arrow(length = unit(0.1, "cm"))))
    } else {
      gene_for_seg <- NULL
    }
    if (nrow(gene_rev) >= 1) {
      # plot forward strand gene
      gene_rev_seg <- list(geom_segment(data = gene_rev, aes(x = V5, y = -(transcript_max -  transcript_min)/15, xend = V4, yend = -(transcript_max - transcript_min)/15),   arrow = arrow(length = unit(0.1, "cm"))))
    } else {
      gene_rev_seg <- NULL
    }
    
    if (isTRUE(label_gene_name) & nrow(gene_for) >= 1) {
      gene_for_seg_name <- list(geom_text_repel(data = gene_for, aes(x = V4,  y = -(transcript_max - transcript_min)/25, label = V9), size = 2, angle = 25))
      
    } else {
      gene_for_seg_name <- NULL
    }
    
    if (isTRUE(label_gene_name) & nrow(gene_rev) >= 1) {
      gene_rev_seg_name <- list(geom_text_repel(data = gene_rev, aes(x = V4,   y = -(transcript_max - transcript_min)/12, label = V9), size = 2,   angle = 25))
    } else {
      gene_rev_seg_name <- NULL
    }
    if (any(nrow(gene_rev) >= 1, nrow(gene_for) >= 1)) {
      gene_box <- list(
        geom_segment(aes(x = rep(transcript_min, 2), xend = rep(transcript_max,   2), y = c(-(transcript_max - transcript_min)/12.5, -(transcript_max -  transcript_min)/45), yend = c(-(transcript_max - transcript_min)/12.5,  -(transcript_max - transcript_min)/45))),
        geom_segment(aes(x = c(transcript_min,  transcript_max), xend = c(transcript_min, transcript_max), y = c(-(transcript_max -  transcript_min)/12.5, -(transcript_max - transcript_min)/12.5),   yend = c(-(transcript_max - transcript_min)/45, -(transcript_max -  transcript_min)/45))),
        geom_segment(aes(x = c(transcript_min,  transcript_max), xend = c(transcript_min, transcript_max), y = c(-(transcript_max -  transcript_min)/12.5, -(transcript_max - transcript_min)/12.5),  yend = c(-(transcript_max - transcript_min)/11.5, -(transcript_max - transcript_min)/11.5))),
        geom_text(aes(x = c(transcript_min,  transcript_max), y = rep(-(transcript_max - transcript_min)/10.2,   2)), label = c(transcript_min, transcript_max)))
    } else {
      gene_box <- list(NULL)
    }
    # decide whether to rotate x axis
    scale_x <- list(
      scale_x_continuous(limits = c(transcript_min - (transcript_max - transcript_min)/6, transcript_max), breaks = seq(transcript_min, transcript_max, (transcript_max - transcript_min)/4)))
    # label the yaxis
    scale_y_line <- list(
      geom_segment(aes(x = transcript_min - (transcript_max - transcript_min)/30, y = min(pvalue_range), xend = transcript_min -  (transcript_max - transcript_min)/30, yend = max(pvalue_range) * fold)))
    scale_y_ticks <- list(
      geom_segment(aes(x = rep(transcript_min - (transcript_max - transcript_min)/15, n_pvalue_range), y = pvalue_range * fold, xend = rep(transcript_min -  (transcript_max - transcript_min)/30, n_pvalue_range), yend = pvalue_range * fold)))
    scale_y_text <- list(
      geom_text(aes(x = rep(transcript_min - (transcript_max -  transcript_min)/12, n_pvalue_range), y = pvalue_range * fold, label = pvalue_range)))
    # add threshold line
    if (is.null(threshold)) {
      threshold_line <- list(NULL)
    }
    if (all(length(threshold) > 0, threshold <= max(pvalue_range))) {
      threshold_line <- list(geom_segment(aes(x = transcript_min, xend = transcript_max, 
                                              y = threshold * fold, yend = threshold * fold), linetype = "longdash", 
                                          colour = "gray"))
    }
    if (all(length(threshold) > 0, threshold > max(pvalue_range))) {
      threshold_line <- list(NULL)
      print("no -log10(p) pass the threshold, will not draw threshold line")
    }
    
    ld_leadsnp_colour <- list(NULL)
    bottom_trianglLD <- list(NULL)
    # link association and LD for the significant loci link between LD and genic
    # structure
    if (any(isTRUE(leadsnpLD), isTRUE(triangleLD))) {
      
      # ld.matrix <- read_delim(LDfile,delim = '\t',col_names = F,col_types = cols('X1'='c')) %>% rename(chr=X1,pos=X2,rsnum=X3,ref=X4,alt=X5)
      # ld.info <- ld.matrix %>% dplyr::select(chr,pos,rsnum,ref,alt) %>% mutate(Seq=seq_along(chr))
      # ld.matrix <- ld.matrix %>% dplyr::select(-c(chr,pos,rsnum,ref,alt)) %>% as.matrix
      # rownames(ld.matrix) <- ld.info$Seq
      # colnames(ld.matrix) <- ld.info$Seq
      
      ld.info <- ld.info[ld.info$chr==chr & ld.info$pos >= transcript_min & ld.info$pos <= transcript_max,]
      ld.matrix <- ld.matrix[ld.info$Seq,ld.info$Seq]
      
      ## reindex 
      ld.info <- ld.info %>% dplyr::select(chr,pos,rsnum,ref,alt) %>% mutate(Seq=seq_along(chr))
      rownames(ld.matrix) <- ld.info$Seq
      colnames(ld.matrix) <- ld.info$Seq
      
      ld <- ld.matrix
      
      if (ldstatistics == "rsquare") 
        ld <- ld^2
      names(ld) <- ld.info$rsnum
      ld <- reshape2::melt(ld)
      marker_info <- data.frame(index =ld.info$Seq, marker_name = ld.info$rsnum)
      ld$Var1 <- marker_info$marker[match(ld$Var1, marker_info$index)]
      ld$Var2 <- marker_info$marker[match(ld$Var2, marker_info$index)]
      if (ldstatistics == "rsquare") {
        lengend_name = expression(italic(r)^2)
      } else if (ldstatistics == "dprime") {
        lengend_name = expression(D * {
          "'"
        })
      }
      ld <- ld[!is.na(ld$value), ]
      ld <- data.frame(Var1 = c(as.character(ld$Var1), as.character(ld$Var2)), 
                       Var2 = c(as.character(ld$Var2), as.character(ld$Var1)), 
                       value = rep(ld$value, 2), stringsAsFactors = FALSE)
      ld0 <- ld
      marker_pos <- transcript_association[, c("Marker", "Site")]
      ld$Site1 <- marker_pos$Site[match(ld$Var1, marker_pos$Marker)]
      ld$Site2 <- marker_pos$Site[match(ld$Var2, marker_pos$Marker)]
      
      if (isTRUE(leadsnpLD)) {
        if (!is.null(leadsnp)) {
          leadsnp <- leadsnp
        }
        if (is.null(leadsnp)) {
          leadsnp <- as.character(transcript_association[which.min(transcript_association$p), "Marker"])
        }
        ld_leadsnp <- ld[ld$Var1 == leadsnp, ]
        ld_leadsnp <- merge(ld_leadsnp, transcript_association, by.x = "Var2", by.y = "Marker")
        if (length(which(ld_leadsnp$Var1 == leadsnp & ld_leadsnp$Var2 ==  leadsnp)) >= 1) {
          ld_leadsnp <- ld_leadsnp[!(ld_leadsnp$Var1 == leadsnp & ld_leadsnp$Var2 ==  leadsnp), ]
        }
        ld_leadsnp$R2 <- 0.2 * (ld_leadsnp$value%/%0.2 + as.logical(ld_leadsnp$value%/%0.2))
        ld_leadsnp$R2 <- as.character(ld_leadsnp$R2)
        ld_leadsnp$R2[ld_leadsnp$R2 == "0"] = "0.2"
        ld_leadsnp$R2[ld_leadsnp$R2 == "1.2"] = "1"
        
        
        if(is.null(association_file)){
          ld_leadsnp_colour <- list(scale_fill_manual(values = c(`0.2` = colour02, `0.4` = colour04, `0.6` = colour06, `0.8` = colour08, `1` = colour10), labels = c("0-0.2","0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"), name = lengend_name))
        }else{
          ld_leadsnp_colour <- list(
            geom_point(data = ld_leadsnp, aes(Site2, -log10(p) * fold, fill = R2), shape = 21, colour = "black",size=2), 
            scale_fill_manual(values = c(`0.2` = colour02, `0.4` = colour04, `0.6` = colour06, `0.8` = colour08, `1` = colour10), labels = c("0-0.2","0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"), name = lengend_name))
        }
        
      }
      if (!isTRUE(leadsnpLD)) {
        ld_leadsnp_colour <- list(NULL)
      }
      
      marker_number = dim(ld.info)[1]
      length = (transcript_max - transcript_min)
      distance = 0.5 * length/(marker_number - 1)
      
      ld <- ld0
      locib <- rep(1:(marker_number - 1), (marker_number - 1):1)
      locia <- sequence((marker_number - 1):1)
      marker_pair <- 1:length(locia)
      center_x <- distance * (locia + locia + locib - 2)
      center_y <- -locib * distance
      upper_center_x <- center_x
      upper_center_y <- center_y + distance
      lower_center_x <- center_x
      lower_center_y <- center_y - distance
      left_center_x <- center_x - distance
      left_center_y <- center_y
      right_center_x <- center_x + distance
      right_center_y <- center_y
      poly_data <- data.frame(
        group = rep(marker_pair, 4), 
        x = c(upper_center_x,right_center_x, lower_center_x, left_center_x) + transcript_min, 
        y = c(upper_center_y, right_center_y, lower_center_y, left_center_y) - 4 * max(pvalue_range) * fold/30, 
        label = rep(c(1, 2, 3, 4), each = length(upper_center_x)),stringsAsFactors = FALSE)
      
      
      poly_data$marker1 <- rep(locia, 4)
      poly_data$marker2 <- rep(locia + locib, 4)
      
      
      marker_index <- data.frame(rs=ld.info$rsnum, marker_number =ld.info$Seq)
      
      poly_data$Var1 <- marker_index$rs[match(poly_data$marker1, marker_index$marker_number)]
      poly_data$Var2 <- marker_index$rs[match(poly_data$marker2, marker_index$marker_number)]
      
      
      poly_data$value <- ld$value[match(paste0(poly_data$Var1, "/", poly_data$Var2), 
                                        paste0(ld$Var1, "/", ld$Var2))]
      
      poly_data$R2 <- 0.2 * (poly_data$value%/%0.2 + as.logical(poly_data$value%/%0.2))
      poly_data$R2 <- as.character(poly_data$R2)
      poly_data$R2[poly_data$R2 == "0"] = "0.2"
      poly_data$R2[poly_data$R2 == "1.2"] = "1"
      
      if (!isTRUE(leadsnpLD)) {
        bottom_trianglLD = list(
          geom_polygon(data = poly_data, aes(group = group, x = x, y = y - (transcript_max - transcript_min)/50, fill = R2)),  
          scale_fill_manual(
            values = c(`0.2` = colour02, `0.4` = colour04,  `0.6` = colour06, `0.8` = colour08, `1` = colour10),
            labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"), name = lengend_name))
      }
      if (isTRUE(leadsnpLD)) {
        bottom_trianglLD = list(geom_polygon(data = poly_data, aes(group = group, x = x, y = y - (transcript_max - transcript_min)/50, fill = R2)))
      }
    }
    if (!isTRUE(triangleLD)) {
      bottom_trianglLD <- list(NULL)
    }
  }
  # link line from significant loci to the strucuture
  if (!is.null(link2gene) & !is.null(threshold)) {
    link_association_structure <- transcript_association[transcript_association$Marker %in% link2gene$rs, ]
    if (dim(link_association_structure)[1] == 0) {
      print("no matched locis, will not draw linking line")
      link_asso_gene <- list(NULL)
      threshold_line <- list(NULL)
    }
    if (dim(link_association_structure)[1] > 0) {
      link_number <- dim(link_association_structure)[1]
      link_asso_gene <- list(geom_segment(data = link_association_structure, aes(x = Site, xend = Site, y = rep(-max(pvalue_range) * fold/30,  link_number), yend = -log10(p) * fold), linetype = "longdash",  colour = "red"))
    }
  }
  if (is.null(link2gene) & is.null(threshold)) {
    print("threshold acquired")
    link_asso_gene <- list(NULL)
    link_LD_genic_structure <- list(NULL)
  }
  if (is.null(link2gene) & !is.null(threshold)) {
    link_association_structure <- transcript_association[-log10(transcript_association$p) >=  threshold, ]
    link_association_structure <- link_association_structure[!duplicated(link_association_structure$p), 
    ]
    if (dim(link_association_structure)[1] == 0) {
      print("no -log10(p) pass the threshold, will not draw link")
      link_asso_gene <- list(NULL)
      threshold_line <- list(NULL)
    }
    if (dim(link_association_structure)[1] > 0) {
      link_association_structure <- transcript_association[-log10(transcript_association$p) >=  threshold, ]
      link_association_structure <- link_association_structure[!duplicated(link_association_structure$p), ]
      link_number <- dim(link_association_structure)[1]
      link_asso_gene <- list(geom_segment(data = link_association_structure,  aes(x = Site, xend = Site, y = rep(-max(pvalue_range) * fold/30,  link_number), yend = -log10(p) * fold), linetype = "longdash", colour = "red"))
    }
  }
  # add linking line to link gene structure and LD matrix
  if (isTRUE(triangleLD)) {
    if (is.null(link2gene) & is.null(link2LD)) {
      if(!is.null(threshold)) {
        link_association_structure <- transcript_association[-log10(transcript_association$p) >=  threshold, ]
        link_association_structure <- link_association_structure[!duplicated(link_association_structure$p), ]
        link_number <- dim(link_association_structure)[1]
        
        link_asso_gene <- list(
          geom_segment(
            data = link_association_structure, 
            aes(x = Site, xend = Site, y = rep(-max(pvalue_range) * fold/30, link_number), yend = -log10(p) * fold),
            linetype = "longdash", colour = "red"))
        marker_axis_LD_x <- transcript_min + (seq(1:marker_number) - 1) * 2 * distance
        marker_axis_genic_x <- ld.info$pos
        marker_axis_LD_y <- rep(-5 * max(pvalue_range) * fold/30, marker_number)
        marker_axis_genic_y <- rep(-max(pvalue_range) * fold/30, marker_number)
        link_ld_data <- data.frame(x = marker_axis_LD_x, xend = marker_axis_genic_x, 
                                   y = marker_axis_LD_y, yend = marker_axis_genic_y)
        link_ld_data <- link_ld_data[link_ld_data$xend %in% link_association_structure$Site, ]
        link_LD_genic_structure <- geom_segment(data = link_ld_data, aes(x = x,  xend = xend, y = y, yend = yend), colour = "red", linetype = "longdash")
      }
    }
    if (!is.null(link2gene) & !is.null(link2LD)) {
      link_association_structure <- transcript_association[transcript_association$Marker %in%  link2LD$rs, ]
      link_number <- dim(link_association_structure)[1]
      link_asso_gene <- list(geom_segment(data = link_association_structure,  aes(x = Site, xend = Site, y = rep(-max(pvalue_range) * fold/30, link_number), yend = -log10(p) * fold), linetype = "longdash",  colour = "red"))
      marker_axis_LD_x <- transcript_min + (seq(1:marker_number) - 1) * 2 * distance
      marker_axis_genic_x <- hapmap_ld$pos
      marker_axis_LD_y <- rep(-4.5 * max(pvalue_range) * fold/30, marker_number)
      marker_axis_genic_y <- rep(-max(pvalue_range) * fold/30, marker_number)
      link_ld_data <- data.frame(x = marker_axis_LD_x, xend = marker_axis_genic_x, 
                                 y = marker_axis_LD_y, yend = marker_axis_genic_y)
      link_ld_data <- link_ld_data[link_ld_data$xend %in% link_association_structure$Site, ]
      link_LD_genic_structure <- geom_segment(data = link_ld_data, aes(x = x, xend = xend, y = y, yend = yend), colour = "red", linetype = "longdash")
    }
  }
  if (!is.null(link2gene) & is.null(link2LD)) {
    link_LD_genic_structure <- list(NULL)
  }
  if (is.null(link2gene) & !is.null(link2LD)) {
    link_LD_genic_structure <- list(NULL)
  }
  if (!isTRUE(triangleLD)) {
    link_LD_genic_structure <- list(NULL)
  }
  
  if(is.null(trait)){
    labelx = "atop(-log[10]*italic(P))"
  }else{
    #labelx = paste0("atop(-log[10]*italic(P)), Trait=",trait)
    labelx = "atop(-log[10]*italic(P))"
  }
  y_axis_text <- list(geom_text(aes(x = transcript_min - (transcript_max - transcript_min)/6.5, y = mean(pvalue_range) * fold), label = labelx,  parse = TRUE, angle = 90))
  #"atop(-log[10]*italic(P)[observed])"
  if (isTRUE(triangleLD)) {
    xtext <- list(geom_text(aes(x = (transcript_max + transcript_min)/2,  y = min(poly_data$y) - 10 * distance, label = paste0("Position on chr.",  chr, " (bp)"))))
  } else {
    xtext <- list(geom_text(aes(x = (transcript_max + transcript_min)/2,  y = -(transcript_max - transcript_min)/10, label = paste0("Position on chr.", chr, " (bp)"))))
  }
  # add shape for the points of leadsnp
  leadsnp2highlight <- transcript_association[transcript_association$Marker == leadsnp, ]
  leadsnp2highlight_list <- list(geom_point(data = leadsnp2highlight, aes(x = Site,  y = -log10(p) * fold), shape = leadsnp_shape, colour = leadsnp_colour, fill = leadsnp_fill, size = leadsnp_size))
  # change the shape,size,colour, and fill for highlighted marker
  if (is.null(marker2highlight)) {
    marker2highlight_list = list(NULL)
  } else {
    marker2highlight <- merge(marker2highlight, transcript_association, by.x = "rs", by.y = "Marker")
    # marker2highlight_list = list(geom_point(data=marker2highlight, aes(Site,
    # -log10(p) * fold, shape=factor(shape), colour=factor(colour),
    # fill=factor(fill), size=factor(size))))
    
    marker2highlight_list = list(
      annotate("point", 
               x = marker2highlight$Site,
               y = -log10(marker2highlight$p) * fold, 
               shape = marker2highlight$shape, 
               colour = marker2highlight$colour, size = marker2highlight$size, 
               fill = marker2highlight$fill))
  }
  if (!is.null(marker2label)) {
    marker2label <- merge(marker2label, transcript_association, by.x = "rs", by.y = "Marker")
    # marker2label_list <-
    # list(annotate('text',x=marker2label$Site,y=-log10(marker2label$p) *
    # fold,label=marker2label$rs,angle=marker2label_angle))
    marker2label_list <- list(
      geom_text_repel(aes(x = marker2label$Site, y = -log10(marker2label$p) * fold, label = marker2label$rs), angle = marker2label_angle,  size = marker2label_size))
  } else {
    marker2label_list <- list(NULL)
  }
  # plot the reduced point if highlighted marker exited
  if (is.null(marker2highlight)) {
    transcript_association = transcript_association
  } else {
    transcript_association = transcript_association[!(transcript_association$Marker %in%  marker2highlight$rs), ]
  }
  
  if(is.null(association_file)){
    plot <- ggplot() + 
      link_asso_gene+
      link_LD_genic_structure + 
      # geom_point(data = transcript_association, aes(Site, -log10(p) * fold),
      #            shape = 21, colour = "black", fill = "black",size=3) +
      ld_leadsnp_colour+
      gene_box + 
      bottom_trianglLD +
      gene_for_seg_name + 
      gene_rev_seg_name + 
      gene_for_seg + 
      gene_rev_seg +
      scale_x + 
      xtext + 
      leadsnp2highlight_list + 
      marker2highlight_list + 
      marker2label_list +
      theme_bw() + 
      theme(legend.key = element_rect(colour = "black"), axis.ticks = element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), text = element_text(size = 15, face = "bold"))
    
  } else{
    plot <- ggplot() + 
      threshold_line +
      link_asso_gene + 
      link_LD_genic_structure + 
      # geom_point(data = transcript_association, aes(Site, -log10(p) * fold),
      #            shape = 21, colour = "black", fill = "black",size=3) +
      ld_leadsnp_colour + 
      gene_box + 
      bottom_trianglLD +
      gene_for_seg_name + 
      gene_rev_seg_name + 
      gene_for_seg + 
      gene_rev_seg +
      scale_x + 
      scale_y_line +
      scale_y_ticks + 
      scale_y_text + 
      y_axis_text + 
      xtext + 
      leadsnp2highlight_list + 
      marker2highlight_list + 
      marker2label_list +
      theme_bw() + 
      theme(legend.key = element_rect(colour = "black"), axis.ticks = element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), text = element_text(size = 15, face = "bold"))
  }
  
  if(is.null(output_file)){
    return(plot)
  }else{
    if(!str_detect(output_file,'png$')) {ggsave(filename = output_file,plot = plot,width = 16,height = 16) }
    if(str_detect(output_file,'png$')) {ggsave(filename = output_file,plot = plot,width = 16,height = 16,dpi = "retina")}
    if(!str_detect(output_file,'pdf$')) {
      output_file <- str_replace(output_file,"\\....$",".pdf")
      ggsave(filename = output_file,plot = plot,width = 16,height = 16)}
  }
  
}

locus_quantification_qtl <- function(gdata,qdata,gdata_queryid=NULL,qdata_queryid=NULL,output_plot=NULL,log2=TRUE){
  #  queryid <- '21:42642038:A:G'
  #queryid2 <- "ENSG00000226496.1"
  
  require(ggstatsplot)

  # box plot
  gdata_samples <- colnames(gdata)[-c(1:4)]
  qdata_samples <- colnames(qdata)[-c(1:4)]
  
  osamples <- intersect(gdata_samples,qdata_samples)
  gdata <- gdata %>% select(1:4,all_of(osamples))
  qdata <- qdata %>% select(1:4,all_of(osamples))
  
  if(!is.null(gdata_queryid)){
    gdata_tmp <- gdata %>% 
      mutate(id=paste(chr,pos,ref,alt,sep=':')) %>%
      filter(id == gdata_queryid) %>% 
      pivot_longer(cols = -c(chr,pos,ref,alt,),names_to = 'Sample',values_to = 'Genotype') %>% 
      select(-c(chr:alt))
  }else{
    gdata_tmp <- gdata %>% 
      mutate(id=paste(chr,pos,ref,alt,sep=':')) %>%
      dplyr::slice(1) %>% 
      pivot_longer(cols = -c(chr,pos,ref,alt,),names_to = 'Sample',values_to = 'Genotype') %>% 
      select(-c(chr:alt))
  }
  
  
  
  if(!is.null(qdata_queryid)){
    qdata_tmp <- qdata %>%
      filter(gene_id == qdata_queryid) %>% 
      pivot_longer(cols = -c(chr,start,end,gene_id),names_to = 'Sample',values_to = 'Quantification') %>% 
      select(-c(chr:end))
  }else{
    qdata_tmp <- qdata %>%
      dplyr::slice(1) %>% 
      pivot_longer(cols = -c(chr,start,end,gene_id),names_to = 'Sample',values_to = 'Quantification') %>% 
      select(-c(chr:end))
  }
  
  if(log2){
    qdata_tmp <- qdata_tmp %>% mutate(Quantification=log2(Quantification+0.1))
  }
  
  tmpdata <- gdata_tmp %>% left_join(qdata_tmp)
  
  p <- ggbetweenstats(
    data = tmpdata,
    x = Genotype,
    y = Quantification,
    #ggtheme = hrbrthemes::theme_ipsum_pub(),
  )
  
  if(is.null(output_plot)){
    return(p)
  }else{
    ggsave(filename = output_plot,plot = p,width = 6,height = 8)
  }
  
}




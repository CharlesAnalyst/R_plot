getIntersection <- function(df_ase, df_snpAnn, snpAnnType="GWAS", snpAnnCol2Sel="p"){
  if (snpAnnType == "GWAS") {
    df_snpAnn <- updatepvaleq0 (df_snpAnn)
    df_snpAnn <- getMinPval ( df_snpAnn)
  }
  df_inter <- getHetOverlapRef(df_ase, df_snpAnn, snpAnnCol2Sel )
  cat("\n No. of SNPs common =", nrow(df_inter ), "\n")
  
  return(df_inter)
}

#' Get the intersection for MAE -of the ASE calling, SNP annotation1
#'  and SNP annotation2
#'
#'The function takes as input 2 data frames, one with the
#'intersection of the ASE calling dataset and the SNP annotation1
#'and the second the intersection of the ASE calling and SNP
#'annotation2 dataset. These can be obtained by calling the
#'getIntersection(). Each of the data frame objects must have
#'a column named 'cmp.col', e.g. containing
#'the rsid, chr:pos values, that is used to merge the data.
#'
#' @param df_aseSnpAnn1 dataframe object containing the intersection
#' of the ASE calling dataset and SNP annotation1 -the output of
#' getIntersection()
#' @param df_aseSnpAnn2 dataframe object containing the intersection
#' of the ASE calling dataset and SNP annotation2 -the output of
#' getIntersection()
#' @return a dataframe object with SNPs common to df_aseSnpAnn1 and
#' df_aseSnpAnn2
#' @export
#'
#' @examples
getIntersectionMae <- function(df_aseSnpAnn1, df_aseSnpAnn2){
  return( merge( df_aseSnpAnn1, df_aseSnpAnn2, by="cmp.col" ) )
}


updatepvaleq0 <- function(df) {
  #if pval==0, assign the pval as (next min pval)*1e-100
  if( length(which(df$p==0) ) >0 ) {
    next.min <- min( df[ which(df$p != 0),"p" ] )
    new.min <- next.min*1e-100
    df[ which(df$p==0), "p"] <- new.min
  }
  return(df)
}

getMinPval<- function(df) {
  #function returns a dataframe with the min pval per SNP, if a SNP has > 1 pval
  recs <-  data.frame( p =  with( df,( tapply( p, cmp.col, function(z) { min(z, na.rm=T ) })) ) , stringsAsFactors=F)
  recs$cmp.col <- rownames(recs)
  recs$p <- as.numeric( recs$p )
  rownames(recs) <- NULL
  return(recs)
}

getHetOverlapRef <- function( df, p.df, field2sel) {
  #the function returns a dataframe with the overlapping SNPs i.e the SNPs in df present in p.df,
  #and the field2sel col e.g pval from the pval.df using the cmp.col
  # In the df returned the fields will be same as df with an additional column of pval (p)
  #p.df <- p.df[ , c("cmp.col", "p") ]
  p.df <- p.df[ , c("cmp.col", field2sel) ]
  overlap  <- merge( df , p.df, by="cmp.col")   #
  df <- overlap
  ##with( overlap,( tapply( p, cmp.col, function(z) { length(z) })) ) # to check nPvals present per SNP from the overlap
  return(df)
}


# load("testDataASE.rda")
# load("testDataGWAS.rda")
# asegwas <- getIntersection(dfAseData, dfGwasData)
# load("testDataASE.rda")
# load("testDataGWAS.rda")
# asegwas <- getIntersection(dfAseData, dfGwasData)
library(ERASE)
library(testthat)
library(optparse)
# setwd("/Charles/project/ASm6A/ASm6A/ASm6A_link_ALS/ERASE/")

option_list = list(
  make_option(c("-a", "--asm6A"), action="store", default=NA, type='character',
              help="ASm6A file"),
  make_option(c("-g", "--gwas"), action="store", default=NA, type='character',
              help="GWAS annotation file"))
opt = parse_args(OptionParser(option_list=option_list))

asm6a <- opt$a # "brain_IP_3.ASm6A.txt"
df_asm6a = read.table(asm6a, header = TRUE)
sample <- basename(asm6a)
# snpAnn = "/Charles/project/ASm6A/data/GWAS_summary/ALS/data_1/summary_statistics/total/als.sumstats.meta.total.GWASanno.txt"
# snpAnn = "/Charles/project/ASm6A/data/GWAS_summary/ALS/data_2/als.total.GWASanno.txt"
# snpAnn = "/Charles/project/ASm6A/data/GWAS_summary/UKB/wes/BMI/BMI.GWASanno.txt"
snpAnn = opt$g # "/Charles/project/ASm6A/data/GWAS_summary/WHR/CombinedSexesAllAncestry/WHRadjBMI.GWASanno.txt"
df_snpAnn = read.table(snpAnn, header = TRUE)
###
asegwas <- getIntersection(df_asm6a, df_snpAnn)
###
asegwas$neglogpval <- -log10(asegwas$p)
pval <- randomization(df_sigASE_SNPann = asegwas[asegwas$qvalue < 0.05,], # min.FDR
                      df_nonASE_SNPann = asegwas[asegwas$qvalue >= 0.05,], # min.FDR
                      colname_rankSNPann = "neglogpval",
                      colname_chk4distr = "depth", # avgReads
                      outFilePrefix = sample, # "SAE",
                      nIterations = 5 ,
                      binwidth = 100, # 10
                      mode="SAE",
                      seedValue=12345)

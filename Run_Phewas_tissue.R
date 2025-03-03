library(tidyr)
library(ggplot2)
library(parallel)
library(PheWAS)
library(stringr)
library(dplyr)


# read covariate
cov_all<-read.table("GE_DEMO_01102023_20230110.csv", header=T,sep="|",colClasses=c("character","character","integer","character","character","integer","integer","integer"))

#cov_all
covColCount<-ncol(cov_all)
cov_all$AGE_LAST_VISIT <- cov_all$YEAR_LAST_VISIT - cov_all$YEAR_OF_BIRTH

ncol(cov_all)
nrow(cov_all)
head(cov_all)

# generate phewas table
phecode_2023 <- read.table("GE_PHECODE_01102023_MEGA_20230110.csv",sep="|",header=T,colClasses=c("character","character","integer","character","character"))
#must have correct type for "PHECODE" as "character", otherwise it will fail on making phe_table
colnames(phecode_2023)<-c("GRID","PHECODE","COUNTS","FIRST_DX","LAST_DX")

genotype<- read.table("imputed_gene.txt",header = TRUE, sep = "\t")
name<-colnames(genotype)
lastSnpIndex<-ncol(genotype)

print(name) # check how many SNPs

ncol(genotype)
nrow(genotype)
head(genotype)

#find white people's id
white_cohort = read.csv("white_cohort.csv",encoding="UTF-8")
id <- subset(white_cohort,select="FID")
names(id)[1] <- "GRID"

ncol(id)
nrow(id)
head(id)

all <- merge(phecode_2023,id)

ncol(all)
nrow(all)
head(all)

all <- phecode_2023 %>% inner_join(id) %>% dplyr::select("GRID","PHECODE","COUNTS")%>%glimpse()
ncol(all)
nrow(all)
head(all)

phe_table<-createPhewasTable(all, min.code.count = 2, add.exclusions=T, translate=F) ### Note, if your phenotype is icd9 code, please use translate=T

ncol(phe_table)
nrow(phe_table)
head(phe_table)


# Read the data from the file
PC <- read.table("mega_white_common.pcs", header = FALSE)
# Keep only the 2-5 columns
PC <- PC[, 2:5]
# Rename the columns
names(PC) <- c("GRID", "PC1", "PC2", "PC3")

# merge cov_all and PC
cov_all <- merge(cov_all, PC, by = "GRID")

ncol(cov_all)
nrow(cov_all)
head(cov_all)



# Begin the loop
for (i in 3:lastSnpIndex) {
  # Read the genotype file only once if it's the same for every iteration (move outside the loop if appropriate)
  genotype <- read.table("imputed_gene.txt",header = TRUE, sep = "\t")
  names(genotype)[1] <- "GRID"
  name <- colnames(genotype)
  print(print(paste("Processing SNP:", name[i])))
  
  snpname <- name[i]
  
  # Perform operations for the valid SNP
  # Example: Create a directory for the SNP if it doesn't exist
  if (!file.exists(snpname)) {
    dir.create(snpname)
  }
  
  # Subset genotype data, merge with other tables, perform analysis, etc.
  # This section can include your data processing and analysis code.
  genotype<-subset(genotype,select=c(GRID, i))
  
  #genotype
  print(print(paste("N of rows in genotype:", nrow(genotype))))
  print(print(paste("N of columns in genotype:", ncol(genotype))))
  
  # merge tables
  all_data<-merge(phe_table, genotype)

  all_data<-merge(all_data, cov_all)
  
  all_data_w <- all_data %>% dplyr::filter(RACE == 'W') %>% glimpse
  
  
  # read merged table
  allColCount<-ncol(all_data)
  print(print(paste("allColCount:", allColCount)))
  
  #count how many phecode
  phecodeCount<-allColCount-covColCount
  print(print(paste("phecodeCount:", phecodeCount)))
  
  phenotypes<-names(all_data_w[2:ncol(phe_table)])
  predictors<-name[i]
  covariates<-c("AGE_LAST_VISIT","GENDER","EHR_LENGTH","PC1","PC2","PC3")
  output<-phewas(phenotypes, predictors, data=all_data_w, covariates)
  # now lets annotate our phewas
  annotated<-addPhecodeInfo(output)
  
  tabletitle=str_c(snpname,"/TG_w_phewas_",snpname,".txt")
  write.table(output, file=tabletitle,quote=F,row.names=F, sep="\t")
  
  annotatedtitle=str_c(snpname,"/TG_w_phewas_",snpname,"_annotated.txt")
  write.table(annotated, file=annotatedtitle,quote=F,row.names=F, sep="\t")
    
}





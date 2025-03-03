# Load necessary libraries
library(PheWAS)
library(tidyr)
library(ggplot2)
library(parallel)
library(stringr)
library(dplyr)
library(GBJ)

# Step 1: Read covariate data
cov_all <- read.table("GE_DEMO_01102023_20230110.csv", header = TRUE, sep = "|", 
                      colClasses = c("character", "character", "integer", 
                                     "character", "character", "integer", "integer", "integer"))

cov_all$AGE_LAST_VISIT <- cov_all$YEAR_LAST_VISIT - cov_all$YEAR_OF_BIRTH

# Step 2: Read white cohort IDs
white_cohort <- read.csv("white_cohort.csv", encoding = "UTF-8")
id <- subset(white_cohort, select = "FID")
names(id)[1] <- "GRID"

# Step 3: Read principal components (PC)
PC <- read.table("mega_white_common.pcs", header = FALSE)
PC <- PC[, 2:5]  # Keep only columns 2-5
names(PC) <- c("GRID", "PC1", "PC2", "PC3")

# Merge covariates and PCs
cov_all <- merge(cov_all, PC, by = "GRID")
cov_all <- cov_all %>%
  select(GRID, AGE_LAST_VISIT, GENDER, EHR_LENGTH, PC1, PC2, PC3)

cov <- cov_all %>%
  filter(GRID %in% id$GRID)

# Step 4: Read genotype data
genotype <- read.table("imputed_gene.txt", header = TRUE, sep = "\t")
names(genotype)[1] <- "GRID"
genotype <- genotype %>%
  select(GRID, chr19_updated_combo_Artery_Tibial,chr19_updated_combo_Kidney_Cortex, chr19_updated_combo_Heart_Atrial_Appendage) %>%
  filter(GRID %in% id$GRID)

# Step 5: Read phecode table
phe_table <- read.csv("phetable_W.csv", header = TRUE)

# List of phecodes to analyze
phecodes <- c("X272.1", "X272.11", "X272.12", "X272.13", "X401.1", "X411.3", 
              "X411.4", "X411.8", "X250.2", "X274.1", "X290.1", "X290.11", 
              "X290.16", "X290.2", "X339", "X340", "X386.9", "X577", 
              "X577.1", "X577.2", "X591")

# Initialize a dataframe to store results
results_summary <- data.frame(
  Phecode = character(), 
  Test_Stat_1 = numeric(), 
  Test_Stat_2 = numeric(), 
  Test_Stat_3 = numeric(),
  GBJ_Statistic = numeric(), 
  GBJ_Pvalue = numeric(), 
  Error_Code = numeric(),
  stringsAsFactors = FALSE
)

# Step 6: Loop through each phecode
for (phecode in phecodes) {
  # Filter phe_table for the current phecode
  phe_temp <- phe_table %>%
    select(GRID, !!sym(phecode)) %>%  # Dynamically select the phecode column
    filter(GRID %in% id$GRID) %>%
    filter(!is.na(!!sym(phecode)))    # Remove NA values
  
  # Rename the column to a standard name for GLM
  names(phe_temp)[2] <- "phecode_val"
  
  # Find common GRID values
  common_GRID <- Reduce(intersect, list(cov$GRID, genotype$GRID, phe_temp$GRID))
  
  # Subset datasets
  cov_sub <- cov[cov$GRID %in% common_GRID, ]
  geno_sub <- genotype[genotype$GRID %in% common_GRID, ]
  phe_sub <- phe_temp[phe_temp$GRID %in% common_GRID, ]
  
  # Merge all data
  demo_table <- merge(cov_sub, phe_sub, by = "GRID")
  full_data <- merge(demo_table, geno_sub, by = "GRID")
  
  # Convert GENDER to numeric
  full_data <- full_data[full_data$GENDER %in% c("F", "M"), ]
  full_data$GENDER <- ifelse(full_data$GENDER == "F", 1, 0)
  
  # Logistic regression
  null_mod <- glm(phecode_val ~ AGE_LAST_VISIT + GENDER + EHR_LENGTH + PC1 + PC2 + PC3, 
                  data = full_data, family = binomial(link = "logit"))
  
  # Prepare genotype data
  genotype_data <- as.matrix(full_data[, c("chr19_updated_combo_Artery_Tibial", "chr19_updated_combo_Kidney_Cortex", " chr19_updated_combo_Heart_Atrial_Appendage")])
  
  # Calculate GBJ statistics
  log_reg_stats <- calc_score_stats(null_model = null_mod, factor_matrix = genotype_data, link_function = "logit")
  cor_Z <- log_reg_stats$cor_mat
  score_stats <- log_reg_stats$test_stats
  
  result <- GBJ(test_stats = score_stats, cor_mat = cor_Z)
  
  # Extract test stats values
  test_stat_1 <- score_stats[1]  # First test statistic
  test_stat_2 <- score_stats[2]  # Second test statistic
  test_stat_3 <- score_stats[3]
  
  # Store results
  results_summary <- rbind(results_summary, 
                           data.frame(
                             Phecode = phecode, 
                             Test_Stat_1 = test_stat_1, 
                             Test_Stat_2 = test_stat_2, 
                             Test_Stat_3 = test_stat_3,
                             GBJ_Statistic = result$GBJ, 
                             GBJ_Pvalue = result$GBJ_pvalue, 
                             Error_Code = result$err_code
                           ))
  
}

# Step 7: Save the summary table to a file
write.csv(results_summary, "GBJ_results_ANGPTL4.csv", row.names = FALSE)

# Print the results summary
print(results_summary)



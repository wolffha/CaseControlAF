# Helper functions

# Function to calculate the p-value of a wald statistic assuming normal distribution
calc_waldP <- function(beta, SE) {
  z <- abs(beta/SE)
  p_val <- 2*(1-pnorm(z))
  return(p_val)
}

# function to calculate OR, SE, and p-value using a 2x2 table
# assumed to be using the AF of trait 1, AF of trait 2, and sample sizes
# DefFac is used to scale the allele counts calculated based on sample overlap
calc_ORStat <- function(AF1, N1, AF2, N2, DefFac) {
  # first calculate the 2x2 tables of allele counts
  a = AF1 * 2 * N1 * DefFac[1]
  b = (1-AF1) * 2 * N1 * DefFac[1]
  c = AF2 * 2 * N2 * DefFac[2]
  d = (1-AF2) * 2 * N2 * DefFac[2]
  
  OR <- (a*d)/(b*c)
  SE <- sqrt(1/a + 1/b + 1/c + 1/d)
  p <- calc_waldP(log(OR), SE)
  
  return(data.frame(OR = OR, SE = SE, p = p))
}

# main helper function for CC-GWAS
# uses the case and control AF for each of the two traits as well as the case and control
# sample sizes of each trait
calc_ccgwas <- function(trait1_AF_case, trait2_AF_case, trait1_N_case, trait2_N_case,
                       trait1_AF_control, trait2_AF_control, trait1_N_control, trait2_N_control) {
  # print(head(trait1_AF_case))
  # print(head(trait1_AF_control))
  # print(trait1_N_case)
  # first calculate the DefControl and DefCase - these are all 1 if there's no sample overlap
  # for now we will assume no sample overlap
  # TODO: implement sample overlap correction
  N_caseInCase <- matrix(data = c(trait1_N_case + trait2_N_case, 0, 0, trait1_N_case + trait2_N_case), nrow = 2)
  N_controlInControl <- matrix(data = c(trait1_N_control + trait2_N_control, 0, 0, 
                                        trait1_N_control + trait2_N_control), nrow = 2)
  N_caseInControl <- matrix(data = c(0, 0, 0, 0), nrow = 2)

  # now we can calculate the default factor for trait 1 and 2 in cases and controls
  DefCase <- c((N_caseInCase[1,1]/(N_caseInCase[1,1] + N_caseInCase[1,2])),
               (N_caseInCase[2,2]/(N_caseInCase[2,2] + N_caseInCase[2,1])))
  DefControl <- c((N_controlInControl[1,1]/(N_controlInControl[1,1] + N_controlInControl[1,2])),
                  (N_controlInControl[2,2]/(N_controlInControl[2,2] + N_controlInControl[2,1])))
   
  # now we calculate the OR, SE, and p for cases and controls separately
  # controls are used for analysis of stratification 
  # overall final CC-GWAS result OR is exp(log(OR_case) - log(OR_control))
  OR_stat_control <- calc_ORStat(trait1_AF_control, trait1_N_control, trait2_AF_control, trait2_N_control, DefControl)
  OR_stat_case <- calc_ORStat(trait1_AF_case, trait1_N_case, trait2_AF_case, trait2_N_case, DefCase)
  
  OR <- exp(log(OR_stat_case$OR)-log(OR_stat_control$OR))
  SE <- OR_stat_case$SE
  p <- calc_waldP(log(OR), SE)
  
  return(data.frame(OR = OR, SE = SE, p = p,
                    OR_control = OR_stat_control$OR,
                    SE_control = OR_stat_control$SE))
}


#' @title CC_GWAS
#' @description
#' This is a function to run a case-case GWAS using the framework introduced in the ReACt package. The method will first filter out any SNPs with NA in the case or control AF. It will then remove any SNPs that are not present in both datasets. Finally, it will remove duplicate SNPs, by defualt keeping the first record. Then the CC-GWAS will be run. 
#' 
#' @param trait1_data a dataframe containing the chromosome, position, case AF, and control AF for each SNP for trait 1
#' @param trait2_data a dataframe containing the chromosome, position, case AF, and control AF for each SNP for trait 2
#' @param trait1_N_case the number of cases in trait 1 dataset
#' @param trait1_N_control the number of controls in trait 1 dataset
#' @param trait2_N_case the number of cases in trait 2 dataset
#' @param trait2_N_control the number of controls in trait 2 dataset
#' @param AF_case_colname a string with the name for the column containing the case AF in both trait 1 and trait 2 datasets; default: "AF_case"
#' @param AF_control_colnamea a string with the name for the column containing the control AF in both trait 1 and trait 2 datasets; default: "AF_control"
#' @param chromosome_colname a string with the name for the column containing the SNP chromosome number in both trait 1 and trait 2 datasets; default: "chr"
#' @param position_colname a string with the name for the column containing the SNP position in both trait 1 and trait 2 datasets; default: "pos"
#' @param id_colname a string containing the id for each SNP in both trait 1 and trait 2 datasets; default: NA. If NA, then an ID column will be added and used for data filtering
#' @note
#'  trait1_data and trait2_data need to have the same column names for chromosome, position, and case and control AFs.
#' @return returns a dataframe with the original data (where each dataset's original columns with have .trait1 or .trait2 appended to its column name) and three additional columns for the CC-GWAS results: OR, SE, p
#' equal to the number of variants
#' @import tidyverse
#' @export 
CC_GWAS <- function(trait1_data, trait2_data, 
                    AF_case_colname = "AF_case", AF_control_colname = "AF_control",
                    trait1_N_case, trait1_N_control, trait2_N_case, trait2_N_control,
                    chromosome_colname = "chr", position_colname = "pos", id_colname = NA,
                    A1_colname = "A1", A2_colname = "A2", rsid_colname = "RSID") {
  require(tidyverse)
  # First need to make sure there are no SNPs with NA is the case or control AFs
  NA_case_1 = length(which(is.na(trait1_data[,c(AF_case_colname)])))
  NA_case_2 = length(which(is.na(trait2_data[,c(AF_case_colname)])))
  NA_control_1 = length(which(is.na(trait1_data[,c(AF_control_colname)])))
  NA_control_2 = length(which(is.na(trait2_data[,c(AF_control_colname)])))
  
  NA_rm1 = 0
  NA_rm2 = 0
  
  if(NA_case_1 > 0 | NA_control_1 > 0) {
    trait1_data_filt <- trait1_data[!is.na(trait1_data[,c(AF_case_colname)]) &
                                      !is.na(trait1_data[,c(AF_control_colname)]),]
    NA_rm1 = nrow(trait1_data)-nrow(trait1_data_filt)
    trait1_data <- trait1_data_filt
    rm(trait1_data_filt)
  }
  if(NA_case_2 > 0 | NA_control_2 > 0) {
    trait2_data_filt <- trait2_data[!is.na(trait2_data[,c(AF_case_colname)]) &
                                      !is.na(trait2_data[,c(AF_control_colname)]),]
    NA_rm2 = nrow(trait2_data)-nrow(trait2_data_filt)
    trait2_data <- trait2_data_filt
    rm(trait2_data_filt)
  }
  rm(NA_case_1)
  rm(NA_case_2)
  rm(NA_control_1)
  rm(NA_control_2)

  # now make sure that both datasets have the same variants
  if(is.na(id_colname)) {
    trait1_data$id <- paste0(trait1_data[,c(chromosome_colname)], ":", 
                          trait1_data[,c(position_colname)])
    trait2_data$id <- paste0(trait2_data[,c(chromosome_colname)], ":", 
                          trait2_data[,c(position_colname)])
    if(nrow(trait1_data) < nrow(trait2_data)) {
      overlap <- trait1_data$id[which(trait1_data$id %in% trait2_data$id)]
    } else {
      overlap <- trait2_data$id[which(trait2_data$id %in% trait1_data$id)]
    }
    
    id_colname = "id"
  } else {
    colnames(trait1_data)[which(colnames(trait1_data) == id_colname)] = "id"
    colnames(trait2_data)[which(colnames(trait2_data) == id_colname)] = "id"
    if(nrow(trait1_data) < nrow(trait2_data)) {
      overlap <- trait1_data$id[which(trait1_data$id %in% trait2_data$id)]
    } else {
      overlap <- trait2_data$id[which(trait2_data$id %in% trait1_data$id)]
    }
  }
  overlap_rm1 = nrow(trait1_data)-length(overlap)
  overlap_rm2 = nrow(trait2_data)-length(overlap)
  trait1_data <- trait1_data %>% filter(id %in% overlap)
  trait2_data <- trait2_data %>% filter(id %in% overlap)
  
  # now keep only distinct 
  trait1_data_filt <- trait1_data %>% distinct(id, .keep_all = T)
  trait2_data_filt <- trait2_data %>% distinct(id, .keep_all = T)
  distinct_rm1 = nrow(trait1_data) - nrow(trait1_data_filt)
  distinct_rm2 = nrow(trait2_data) - nrow(trait2_data_filt)
  trait1_data = trait1_data_filt
  rm(trait1_data_filt)
  trait2_data = trait2_data_filt
  rm(trait2_data_filt)

  nRemoved <- data.frame(dataset = c(1, 2),
                         NA_removed = c(NA_rm1, NA_rm1),
                         overlap_removed = c(overlap_rm1, overlap_rm2),
                         distinct_removed = c(distinct_rm1, distinct_rm2))
  print(("Removed the following number of variants due to NAs in AFs, SNPs not in both datasets, or duplicate SNPs:"))
  print(nRemoved)

      # rename columns so I can use dplyr to select desired ones for output
  colnames(trait1_data)[which(colnames(trait1_data) == chromosome_colname)] = "chr"
  colnames(trait1_data)[which(colnames(trait1_data) == position_colname)] = "pos"
  colnames(trait1_data)[which(colnames(trait1_data) == rsid_colname)] = "RSID"
  colnames(trait2_data)[which(colnames(trait2_data) == chromosome_colname)] = "chr"
  colnames(trait2_data)[which(colnames(trait2_data) == position_colname)] = "pos"
  colnames(trait2_data)[which(colnames(trait2_data) == rsid_colname)] = "RSID"
  
  if(AF_case_colname != "AF_case") {
    trait1_data$AF_case <- trait1_data[,c(AF_case_colname)]
    trait2_data$AF_case <- trait2_data[,c(AF_case_colname)]
  }
  if(AF_control_colname != "AF_control") {
    trait1_data$AF_control <- trait1_data[,c(AF_control_colname)]
    trait2_data$AF_control <- trait2_data[,c(AF_control_colname)]
  }
  if(A1_colname != "A1") {
    trait1_data$A1 <- trait1_data[,c(A1_colname)]
  }
  if(A2_colname != "A2") {
    trait1_data$A2 <- trait1_data[,c(A2_colname)]
  }
  
  # make sure trait1 and trait2 are in the same order 
  if(nrow(trait1_data) != nrow(trait2_data)) {
    stop("Number of SNPs doesn't match")
  }
  combined <- trait1_data %>% inner_join(trait2_data, by = c("id" = "id"),
                                         suffix = c(".trait1", ".trait2"))
  
  # run the actual CC-GWAS using the cleaned datasets
  results <- calc_ccgwas(trait1_AF_case = combined$AF_case.trait1,
                         trait2_AF_case = combined$AF_case.trait2,
                         trait1_AF_control = combined$AF_control.trait1,
                         trait2_AF_control = combined$AF_control.trait2,
                         trait1_N_case = trait1_N_case,
                         trait2_N_case = trait2_N_case,
                         trait1_N_control = trait1_N_control,
                         trait2_N_control = trait2_N_control)
  
  toreturn <- combined %>% select(chr.trait1, pos.trait1, RSID.trait1, A1.trait1, A2.trait1,
                                  AF_case.trait1, AF_case.trait2, AF_control.trait1,
                                  AF_control.trait2)
  colnames(toreturn)[1:5] <- c("chr", "pos", "RSID", "A1", "A2")
  toreturn <- toreturn %>% cbind(results)
  return(toreturn)
}

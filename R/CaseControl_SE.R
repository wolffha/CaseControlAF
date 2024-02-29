# Helper function for ReACt method to derive case/control AF using SE
# NOTE: this method exists as a C implementation at <https://github.com/Paschou-Lab/ReAct>
GroupFreq <- function(se, nCase, nControl, OR, freq = 0.0) {
  # using w, x, y, z to maintain consistency with the derivation from ReACt manuscript
  x = 2*nCase
  y = 2*nControl
  w = se^2
  z = OR
  
  
  #inflate w(se) by 1.001, for at max 49 times (usually can be done within 5 iterations. 
  #maximum inflation is 0.050195, ~5%), if tmp2 still < 0 then give up
  for(i in (1:50)) {
    w = w * (1.001^i)
    tmp2 <- (((2*z*y*(1-z)) - (w*x*y*z))^2) - 4*(w*x*z + (1-z)^2)*(y*z*(x + y*z))
    if (tmp2 >= 0) {
      break
    }
  }
  if (tmp2 < 0) { # if discriminate is <0 cannot solve and return 0s
    resFrq <- data.frame(MAF_case = 0.0,
                         MAF_control = 0.0,
                         MAF_pop = 0.0)
  } else { # actually solve the allele frequency
    tmp1 <- w * x * y * z - (2 * y * z * (1 - z))
    tmp2 = tmp2^0.5
    tmp3 = 2*(w*x*z + (1-z)^2)
    
    d1 = (tmp1 - tmp2) / tmp3
    c1 = y - d1
    b1 = x * d1 / ((z * y) - (z * d1) + d1)
    a1 = x - b1
    frq1 = c1 / (c1 + d1)
    
    d2 = (tmp1 + tmp2)/ tmp3
    c2 = y - d2
    b2 = x * d2 / ((z * y) - (z * d2) + d2)
    a2 = x - b2
    frq2 = c2 / (c2 + d2)
    
    flag1 = 0 # indicator to use root d1
    flag2 = 0 # indicator to use root d2

    # need to select the root that allows a, b, c, d (allele counts of the 2x2 table) to all be positive
    if (a1 > 0 & b1 > 0 & c1 > 0 & d1 > 0) { 
      flag1 = 1
      tmpRes1 <- data.frame(res11 = a1,
                            res12 = b1,
                            res21 = c1,
                            res22 = d1)
    }
    if (a2 > 0 & b2 > 0 & c2 > 0 & d2 > 0) {
      flag2 = 1
      tmpRes2 <- data.frame(res11 = a2,
                            res12 = b2,
                            res21 = c2,
                            res22 = d2)
    }
    
    if (flag1 == 1 & flag2 == 0) { # use d1 if d1 is all positive and d2 is NOT
      res = tmpRes1
    } else if (flag1 == 0 & flag2 == 1) { # use d2 if all positive and d1 is NOT
      res = tmpRes2
    } else if (flag1 == 1 & flag2 == 1) { # if both are all positive then
      if (abs(freq - frq1) < abs(freq - frq2)) { #use the larger value as it minimizes the fxn of interest
        res = tmpRes1
      } else {
        res = tmpRes2
      }
    } else { # if it cannot solve due to negatives return all 0s
      res =  data.frame(res11 = 0.0,
                        res12 = 0.0,
                        res21 = 0.0,
                        res22 = 0.0)
    }
    
    # solve for the minor allele frequency using the solved for allele counts
    MAF_case = res$res11 / (res$res11 + res$res12)
    MAF_control = res$res21 / (res$res21 + res$res22)
    MAF_pop = (res$res11 + res$res21) / (res$res11 + res$res12 + res$res21 + res$res22)
    
    # if any evaluates to inf replace with 0
    MAF_case <- ifelse(is.infinite(MAF_case), 0.0, MAF_case)
    MAF_control <- ifelse(is.infinite(MAF_control), 0.0, MAF_control)
    MAF_pop <- ifelse(is.infinite(MAF_pop), 0.0, MAF_pop)
    
    resFrq <- data.frame(MAF_case, MAF_control, MAF_pop)
  }
  return(resFrq)
}

#' @title CaseControl_SE
#' @description This is a function to derive the case and control AFs from GWAS summary statistics when 
#' the user has access to the whole sample AF, the sample sizes, and the OR (or beta).
#' If user has SE instead of sample AF use [CaseControlAF::CaseControl_SE()]
#' This code uses the GroupFreq function adapted from C from <https://github.com/Paschou-Lab/ReAct/blob/main/GrpPRS_src/CountConstruct.c>
#' @param OR a vector of odds ratios
#' @param se a vector of standard errors for the OR calculation *make sure the indices match between the two vectors
#' @param nCase an integer of the number of Case individuals
#' @param nControl an integer of the number of Control individuals
#' @return a dataframe with 3 columns: pCase, pControl, pPop for the estimated AFs for each variant which are the rows
#' @export 
CaseControl_SE <- function(OR, SE, N_case, N_control) {
  res <- data.frame(MAF_case = rep(0, length(OR)),
                    MAF_control = rep(0, length(OR)),
                    MAF_pop = rep(0, length(OR)))
  for(i in 1:length(OR)) {
    res[i,] <- GroupFreq(SE[i], nCase = N_case, nControl = N_control, OR[i])
  }
  
  return(res)
}

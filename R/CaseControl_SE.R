#' @title CaseControl_SE
#' @description This is a function to derive the case and control AFs from GWAS summary statistics when 
#' the user has access to the whole sample AF, the sample sizes, and the OR (or beta).
#' If user has SE instead of sample AF use [CaseControlAF::CaseControl_SE()]
#' This code uses the GroupFreq function adapted from C from <https://github.com/Paschou-Lab/ReAct/blob/main/GrpPRS_src/CountConstruct.c>
#' @param OR a vector of odds ratios
#' @param se a vector of standard errors for the OR calculation *make sure the indices match between the two vectors
#' @param nCase an integer of the number of Case individuals
#' @param nControl an integer of the number of Control individuals
#' @return a dataframe with 3 columns: MAF_case, MAF_control, MAF_pop for the estimated MAFs for each variant 
#' @import genpwr
#' @export 
CaseControl_SE <- function(OR, SE, N_case, N_control) {
  # genpwr package is used to solve for the roots of the quadratic equation
  require(genpwr)

  # this function uses w, x, y, z as the derivation in the ReACt paper does - see their supplement
  w = SE^2
  x = 2*N_case
  y = 2*N_control
  z = OR

  AC_control = rep(0, length(OR))
  AC_case = rep(0, length(OR))

  for(i in 1:length(OR)) {
    # need to make sure the discriminate will be positibe or it won't solve
  #inflate w(se) by 1.001, for at max 49 times (usually can be done within 5 iterations.
  #maximum inflation is 0.050195, ~5%), if disc still < 0 then give up
    for(j in 0:99) {
      w[i] = w[i] * (1.001^j)
      disc = (((2*z[i]*y*(1-z[i])) - (w[i]*x*y*z[i]))^2) - 4*(w[i]*x*z[i] + (1-z[i])^2)*(y*z[i]*(x + y*z[i]))
      if (disc >= 0) {
        break
      } 
    }
    if (disc < 0) { # if discriminate is <0 cannot solve and return 0s
      AC_control[i] = 0
      AC_case[i] = 0
    } else { # now actually solve the quadratic for allele counts (AC)
        # solve for the a, b, and c of the quadratic equation
        # this quadratic is solving for the allele count (AC) of the controls
        # overall their derivation relies on AC (rather than AF) and then calculates AF
        a = (w[i]*x*z[i]) + (1-z[i])^2
        b = 2*y*z[i]*(1-z[i]) - w[i]*x*y*z[i]
        c = y*z[i]*(x + y*z[i])
        
         #find roots of quadratic equation
        AF_control_opts =  quad_roots(a, b, c)
        # in order to select which root, we need to use each option (d1 and d2) to calculate a, b, c of the 2x2 table of allele counts
        d1 = AF_control_opts[1]
        c1 = y - d1
        b1 = (x*d1)/(y*z[i] - z[i]*d1 + d1)
        a1 = x - b1
        
        d2 = AF_control_opts[2]
        c2 = y - d2
        b2 = (x*d2)/(y*z[i] - z[i]*d2 + d2)
        a2 = x - b2

        vec1 = c(a1, b1, c1, d1) # vector of a,b,c,d using root 1
        vec2 = c(a2, b2, c2, d2) # vector of a,b,c,d using root 2

        # if both roots allow for all values to be positive, choose the larger
        if(!any(vec1 < 0) & !(any(vec2 < 0))) {
          if(d1 > d2) { # if d1 is the larger root, then the AC_control = c1 and AC_case = a1
            AC_control[i] = c1
            AC_case[i] = a1
          } else {
            AC_control[i] = c2
            AC_case[i] = a2
          }
        } else if(!any(vec1 < 0)) { # if d1 allows all values of 2x2 to be positive but NOT d2, use c1 and a1
          AC_control[i] = c1
          AC_case[i] = a1
        } else { # if d2 allows all values to be positive but NOT d1, use c2 and a2
          AC_control[i] = c2
          AC_case[i] = a2
        }
      }
    }

  # calculate and return the MAF
  MAF_case = AC_case/(2*N_case)
  MAF_control = AC_control/(2*N_control)
  MAF_pop = (MAF_case*N_case + MAF_control*N_control)/(N_case + N_control)

  return(data.frame(MAF_case = MAF_case,
                    MAF_control = MAF_control,
                    MAF_pop = MAF_pop))
}

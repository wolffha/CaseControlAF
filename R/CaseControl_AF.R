#' @title CaseControl_AF
#' @description
#' This is a function to derive the case and control AFs from GWAS summary statistics when 
#' the user has access to the whole sample AF, the sample sizes, and the OR (or beta).
#' If user has SE instead of sample AF use [CaseControlAF::CaseControl_SE()]
#' 
#' @param N_case the number of cases in the sample
#' @param N_control the number of controls in the sample
#' @param AF_population a vector with the observed sample AFs
#' @param OR a vector with the ORs 
#' @note
#'  Make sure the vectors are in the same order for your variants
#' @return returns a dataframe with two columns (AF_case, AF_control) and rows 
#' equal to the number of variants
#' @import genpwr
#' @export 

### COMMENT ###

CaseControl_AF <- function(N_case, N_control, AF_population, OR){
  
  require(genpwr)
  
  #calculate total sample size
  N_total = N_control+N_case
  
  #set a, b, c of quadratic equation derived in pdf
  a = (N_control/N_case)*(OR-1)
  b = (OR-((N_total/N_case)*AF_population*OR))+((N_control/N_case)+(N_total*AF_population/N_case))
  c = -(N_total/N_case)*AF_population
  
  AF_control <- rep(0, length(a))
  
  for(i in 1:length(a)) {
    #find roots of quadratic eq and choose root between 0 and 1 as AF_control
    AF_control_opts =  quad_roots(a[i], b[i], c[i])
    if(AF_control_opts[1]>1 | AF_control_opts[1]<0){
      AF_control[i] = AF_control_opts[2]
    }else{
      AF_control[i] = AF_control_opts[1]
    } 
  }
  
  
  #calculate AF_case with known relationship shown in pdf
  AF_case = (N_total/N_case)*AF_population - (N_control/N_case)*AF_control
  
  #Output shows case AF first, then control AF
  return(data.frame(AF_case = AF_case,
                    AF_control = AF_control))
}

#' @title CaseControl_AF
#' @description
#' This is a function to derive the case and control AFs from GWAS summary statistics when
#' the user has access to the whole sample AF, the sample sizes, and the OR (or beta).
#' If user has SE instead of sample AF use [CaseControlAF::CaseControl_SE()]
#'
#' @param N_case the number of cases in the sample
#' @param N_control the number of controls in the sample
#' @param OR a vector with the ORs
#' @param AF_population a vector with the observed sample AFs
#'
#' @note
#'  Make sure the vectors are in the same order for your variants
#'
#' @return returns a dataframe with two columns (AF_case, AF_control) and rows
#' equal to the number of variants
#'
#' @author Hayley Wolff (Stoneman), \email{hayley.wolff@cuanschutz.edu}
#'
#' @references https://github.com/wolffha/CaseControlAF
#'
#' @seealso \url{https://github.com/wolffha/CaseControlAF} for further documentation
#'
#' @examples
#' data("sampleDat")
#'
#' nCase_sample = 16550
#' nControl_sample = 403923
#'
#' # get the estimated case and control AFs
#' af_method_results <- CaseControl_AF(N_case = nCase_sample,
#'                                     N_control = nControl_sample,
#'                                     OR = sampleDat$OR,
#'                                     AF_population = sampleDat$true_maf_pop)
#'
#' head(af_method_results)
#'
#' @import genpwr
#' @export
CaseControl_AF <- function(N_case, N_control, OR, AF_population){

  # do input checking
  if(length(OR) != length(AF_population)) {
    stop("ERROR: 'OR' and 'AF_population' vectors need to be the same length")
  }
  if(typeof(OR) != "double") {
    stop("ERROR: 'OR' vector needs to be a vector of numbers (hint: check for NAs)")
  }
  if(typeof(AF_population) != "double") {
    stop("ERROR: 'AF_population' vector needs to be a vector of numbers (hint: check for NAs)")
  }
  if(N_case <= 0) {
    stop("ERROR: 'N_case' needs to be a number > 0")
  }
  if(N_control <= 0) {
    stop("ERROR: 'N_control' needs to be a number > 0")
  }
  if(any(AF_population < 0)) {
    stop("ERROR: 'AF_population' vector cannot contain negative AFs")
  }
  if(any(AF_population > 1)) {
    stop("ERROR: 'AF_population' vector cannot contain values > 1")
  }

  # uses the genpwr package for the quad_roots function to solve for the roots of quadratic
  #require(genpwr)

  #calculate total sample size
  N_total = N_control+N_case

  #set a, b, c of quadratic equation derived in manuscript
  a = (N_control/N_case)*(OR-1)
  b = (OR-((N_total/N_case)*AF_population*OR))+((N_control/N_case)+(N_total*AF_population/N_case))
  c = -(N_total/N_case)*AF_population

  AF_control <- rep(0, length(a))

  for(i in 1:length(a)) {
    #find roots of quadratic eq and choose root between 0 and 1 as AF_control
    AF_control_opts =  genpwr::quad_roots(a[i], b[i], c[i])
    if(AF_control_opts[1]>1 | AF_control_opts[1]<0){
      AF_control[i] = AF_control_opts[2]
    }else{
      AF_control[i] = AF_control_opts[1]
    }
  }


  #calculate AF_case with known relationship shown in manuscript
  AF_case = (N_total/N_case)*AF_population - (N_control/N_case)*AF_control

  #Output shows case AF first, then control AF
  return(data.frame(AF_case = AF_case,
                    AF_control = AF_control))
}

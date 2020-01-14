#' dimsumms_errormodel_random_effect
#'
#' estimate variant specific random error terms
#'
#' @param score_rx replicate fitness scores (#score_rx is GxR matrix with G number of genotypes, R number of replicates)
#' @param sigma_rx replicate count-based errors
#'
#' @return Nothing
#' @export
#' @import data.table
dimsumms_errormodel_random_effect = function(
  score_rx,
  sigma_rx) {
  
  ### random-effect model
  score = rowSums(score_rx * sigma_rx^-2,na.rm=T) / rowSums(sigma_rx^-2,na.rm=T)
  sigma_s2 = 1/(rowSums(!is.na(score_rx))-1) * rowSums((score - score_rx)^2,na.rm=T) #already squared
  #fisher scoring iterations
  for (i in 1:50) {
    weights = 1 / (sigma_rx^2 + sigma_s2)
    score = rowSums(score_rx * weights,na.rm=T) / rowSums(weights,na.rm=T)
    sigma_s2= sigma_s2 * rowSums(weights^2 * (score - score_rx)^2,na.rm=T) / 
      (rowSums(weights,na.rm=T) - rowSums(weights^2,na.rm=T)/rowSums(weights,na.rm=T))
  }
  sigma_s2[is.na(sigma_s2)] = 0
  return(sigma_s2)
}

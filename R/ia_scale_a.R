#' Title
#'
#' @param input_array
#' @param target
#' @param method
#' @param reduce_target
#'
#' @return
#' @export
#'
#' @examples
ia_scale <- function(input_array, target,
                     method = "sum", reduce_target = FALSE){

  if(!(method %in% c("sum", "product"))){
    stop("method must either be 'sum' or 'product'")
  }

  #Check input sizes
  if(length(target) > 1 & !reduce_target){
    warning("Length of `target` is greater than 1. Only first value will be used. Set reduce_target = TRUE to sum over target.")
    tar <- as.numeric(target[1])
  } else if(length(target) > 1 & reduce_target){
    tar <- case_when(
      method == "sum" ~ sum(target, na.rm=TRUE),
      method == "product" ~ prod(target, na.rm=TRUE)
    )
  } else {
    tar <- target
  }

  scale_length <- length(input_array)

  scale_fct <- case_when(
    method == "sum" ~ (tar / sum(input_array)),
    method == "product" ~ (tar^(1/scale_length)) / (prod(input_array)^(1/scale_length))
  )

  return(scale_fct * input_array)

}

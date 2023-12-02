
#' Title
#'
#' @param datatable
#' @param target
#' @param series_start
#' @param series_target
#' @param series_type
#'
#' @return
#' @export
#'
#' @examples
ip_miss <- function(datatable, target, series_start = "value", series_target = "value", series_type = "tar") {

  names(datatable)[names(datatable) == series_start] <- "value"
  names(target)[names(target) == series_target] <- "target_value"

  error <- datatable %>%
    left_join(target) %>%
    do(
      if(series_type != "subtl"){
        group_by_(., .dots = as.list(c(names(target))))
      } else {
        group_by_(., .dots = as.list("target_value"))
      }
    ) %>%
    summarize(sum = sum(value, na.rm=T)) %>%
    ungroup() %>%
    mutate(error = abs(sum - target_value))
}

#' Compare data frame to target subtotals and calculate absolute error.
#'
#' @param datatable
#' @param target_series
#' @param series_start
#' @param series_target
#' @param series_type

ip_miss_a <- function(datatable, target_series, series_start = "value", series_target = "tar1", series_type = "tar") {

  names(datatable)[names(datatable) == series_target] <- "target"
  target_series[target_series == series_target] <- "target"

  error <- datatable %>%
    do(
      if(series_type != "subtl"){
        group_by_(., .dots = as.list(unique(c(target_series, "target"))))
      } else {
        group_by_(., .dots = as.list("target"))
      }
    ) %>%
    summarize(sum = sum(value, na.rm=T)) %>%
    ungroup() %>%
    mutate(error = abs(sum - target))

  names(error)[names(error) == "target"] <- series_target

  return(error)
}

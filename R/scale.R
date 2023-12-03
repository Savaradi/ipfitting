
#' Title
#'
#' @param datatable
#' @param target
#' @param series_start
#' @param series_target
#'
#' @return
#' @export
#'
#' @examples
ip_scale <- function(datatable, target, series_start = "value", series_target = "value") {

  names(datatable)[names(datatable) == series_start] <- "value"
  names(target)[names(target) == series_target] <- "target_value"

  datatable <- datatable %>%
    left_join(target) %>%
    group_by_(.dots = as.list(names(target)[!(names(target) %in% c(series_target))])) %>%
    mutate(shr = ifelse(value == 0 | is.na(value), 0, value/sum(value, na.rm=T))) %>%
    ungroup() %>%
    mutate(value = shr * target_value) %>%
    select(-shr, -target_value)

  return(datatable)
}

#' Single target scale function using single data frame.
#'
#' @param datatable
#' @param target_series
#' @param series_start
#' @param series_target
#' @param series_type


ip_scale_a <- function(datatable, target_series, series_start = "value", series_target = "tar1", series_type = "tar") {

  names(datatable)[names(datatable) == series_start] <- "value_temp"
  names(datatable)[names(datatable) == series_target] <- "tar_temp"
  target_series[target_series == series_target] <- "tar_temp"

  datatable <- datatable %>%
    ungroup() %>%
    do(
      if(series_type == "subtl") {
        group_by_(., .dots = as.list("tar_temp"))
      } else {
        group_by_(., .dots = as.list(target_series[!(target_series %in% c(series_target))]))
      }
    ) %>%
    # group_by_(.dots = as.list(target_series[!(target_series %in% c(series_target))])) %>%
    mutate(shr = ifelse(value_temp == 0 | is.na(value_temp), 0, value_temp/sum(value_temp, na.rm=T))) %>%
    ungroup() %>%
    #Discount factor for any missing targets
    #By however much the existing targets raise, lower the other values by that much
    mutate(
      val_to_be_scaled = ifelse(is.na(tar_temp), NA, value_temp),
      proportional_movement =  sum(shr * tar_temp, na.rm=T) / sum(val_to_be_scaled, na.rm=T)
    ) %>%
    mutate(value_temp = ifelse(is.na(tar_temp), value_temp * 1/proportional_movement, shr * tar_temp)) %>%
    select(-shr, -val_to_be_scaled,-proportional_movement)

  names(datatable)[names(datatable) == "value_temp"] <- series_start
  names(datatable)[names(datatable) == "tar_temp"] <- series_target

  return(datatable)
}

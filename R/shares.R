
#' Title
#'
#' @param datatable
#' @param target
#' @param series_start
#' @param series_target
#' @param series_share_of
#' @param sep
#'
#' @return
#' @export
#'
#' @examples
ip_shares_transform <- function(datatable, target, series_start = "value", series_target = "share",
                                series_share_of = "share__of", sep = " + ") {

  names(datatable)[names(datatable) == series_start] <- "value"
  names(target)[names(target) == series_target] <- "share"
  names(target)[names(target) == series_share_of] <- "share__of"

  tars <- target
  #Add column for grouping - unique sets of series included
  tars$type <- apply(tars, 1, function(x){
    series <- names(target)[!is.na(x)]
    series <- series[!(series  %in% c("share", "share__of"))]

    share__of <- x[names(x) == "share__of"]

    return(paste0(share__of, " | ", paste(series, collapse = " ")))
  })

  tars.list <- lapply(unique(tars$type), function(x){
    df <- tars %>% filter(type == x) %>% select(-type)
    df <- df[, colSums(is.na(df)) == 0]
    return(df)
  })

  tars.values <- lapply(seq_along(tars.list), function(i){
    x <- tars.list[[i]] %>%
      select(-share__of)

    share__of <- strsplit(tars.list[[i]]$share__of[1], split=sep)[[1]]

    dat <- datatable  %>%
      group_by_(.dots = as.list(names(x %>% select(-share)))) %>%
      summarize(value_start = sum(value, na.rm=T)) %>%
      ungroup() %>%
      left_join(x, by = names(x %>% select(-share))) %>%
      group_by_(.dots = as.list(names(.)[!(names(.) %in% c("share", "value_start", share__of))])) %>%
      filter(!all(is.na(share))) %>%
      mutate(value_new = share * sum(value_start, na.rm=T)) %>%
      mutate(value = ifelse(!is.na(share), value_new,
                            (sum(value_start, na.rm=T) - sum(value_new, na.rm=T)) * value_start / sum(value_start * is.na(share), na.rm=T))) %>%
      select(-share, -dplyr::starts_with("value_")) #%>%
    #rename_(.dots = setNames("value", paste0("tar_shares__", i)))

  })

  tars.values <- bind_rows(tars.values)

  return(tars.values)

}

#' Scale data frame values to shares over specified series
#'
#' @param datatable A data frame of values.
#' @param series_start The name of the series in \code{datatable} to be converted to shares.
#' @param groups A character vector of series names to group data, excluding from share calculations.
#' @return A summarized data frame with the same dimensionality as \code{datatable}, with values grouped by \code{groups}.
#' @export
ip_shares_calc <- function(datatable, series_start = "value", groups = NULL) {

  names(datatable)[names(datatable) == series_start] <- "value"

  dat <- datatable %>%
    group_by_(.dots = as.list(groups)) %>%
    mutate(share = value / sum(value, na.rm=T)) %>%
    ungroup()

  return(dat)

}

#' Title
#'
#' @param datatable
#' @param series_start
#' @param groups
#'
#' @return
#' @export
#'
#' @examples
ip_shares_tot <- function(datatable, series_start = "value", groups = NULL) {

  names(datatable)[names(datatable) == series_start] <- "value"

  dat <- datatable[, !(names(datatable) %in% c(groups, "value"))]

  dat <- dat %>%
    distinct() %>%
    mutate(value = 1)

  return(dat)

}


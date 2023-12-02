#' Title
#'
#' @param .data
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
group_by_not <- function(.data, ...){
  dont_group <- gsub("~", "", as.character(rlang::quos(...)))

  .data %>%
    dplyr::group_by_(.dots = names(.data)[!(names(.data) %in% dont_group)])
}

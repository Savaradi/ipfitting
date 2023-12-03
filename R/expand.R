
#' Title
#'
#' @param targets
#' @param target.value.names
#' @param names.exclude
#' @param value.set
#' @param value.name
#' @param override.warning
#' @param max.error
#' @param max.iterations
#' @param freeze_cells
#' @param freeze_cells.value.name
#' @param freeze_slice
#' @param freeze_slice.value.names
#' @param minmax_cells
#' @param minmax_cells.value.names
#' @param minmax_slice
#' @param minmax_slice.value.names
#' @param minmax.smash.param
#' @param save.tars
#' @param show.messages
#'
#' @return
#' @export
#'
#' @examples
ip_expand <- function(targets, target.value.names = "value",
                      names.exclude = c("value"), value.set = 1, value.name = "value",
                      override.warning = FALSE,
                      max.error = 0.01, max.iterations = 25,
                      freeze_cells = NULL, freeze_cells.value.name = "value",
                      freeze_slice = NULL, freeze_slice.value.names = "value",
                      minmax_cells = NULL, minmax_cells.value.names = c("value_min", "value_max"),
                      minmax_slice = NULL, minmax_slice.value.names = c("value_min", "value_max"),
                      minmax.smash.param = 1/3,
                      save.tars = FALSE, show.messages = TRUE) {

  df <- ip_create_seed(targets, names.exclude = names.exclude,
                       value.set = value.set, value.name = value.name,
                       override.warning = override.warning) %>%
    ip_fitting(targets, datatable.value.name = value.name, target.value.names = target.value.names,
           max.error = max.error, max.iterations = max.iterations,
           freeze_cells = freeze_cells, freeze_cells.value.name = freeze_cells.value.name,
           freeze_slice = freeze_slice, freeze_slice.value.names = freeze_slice.value.names,
           minmax_cells = minmax_cells, minmax_cells.value.names = minmax_cells.value.names,
           minmax_slice = minmax_slice, minmax_slice.value.names = minmax_slice.value.names,
           minmax.smash.param = minmax.smash.param,
           save.tars = save.tars, show.messages = show.messages)

  return(df)

}

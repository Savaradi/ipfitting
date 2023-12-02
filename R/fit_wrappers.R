#' Title
#'
#' @param datatable
#' @param targets
#' @param assumptions
#' @param datatable.value.name
#' @param target.value.names
#' @param assumption.value.names
#' @param assumption.drop.names
#' @param max.error
#' @param max.iterations
#' @param minmax.smash.param
#' @param save.tars
#' @param show.messages
#'
#' @return
#' @export
#'
#' @examples
ip_fit_sl <- function(datatable, targets, assumptions,
                      datatable.value.name = "value", target.value.names = "value",
                      assumption.value.names = c("value", "value_min", "value_max"),
                      assumption.drop.names = c("Notes"),
                      max.error = 0.01, max.iterations = 25,
                      minmax.smash.param = 1/3,
                      save.tars = FALSE, show.messages = TRUE) {

  assumptions_pro <- ip_load_assumptions(assumptions,
                                         freeze.name = assumption.value.names[1], minmax.name = assumption.value.names[2:3],
                                         drop.names = assumption.drop.names)

  df <- datatable %>%
    ip_fit(targets,
           datatable.value.name = datatable.value.name, target.value.names = target.value.names,
           max.error = max.error, max.iterations = max.iterations,
           freeze_cells = assumptions_pro[["freeze_cells"]],
           freeze_slice = assumptions_pro[["freeze_slice"]],
           minmax_cells = assumptions_pro[["minmax_cells"]],
           minmax_slice = assumptions_pro[["minmax_slice"]],
           minmax.smash.param = minmax.smash.param,
           save.tars = save.tars, show.messages = show.messages)

  return(as.data.frame(df))

}

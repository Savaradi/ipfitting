
#' Title
#'
#' @param targets
#' @param seed
#' @param target.value.names
#' @param seed.value.name
#' @param max.error
#'
#' @return
#' @export
#'
#' @examples
check_targets <- function(targets, seed = NULL,
                          target.value.names = "value", seed.value.name = "value",
                          max.error = 0.01) {

  #Warnings
  if(is.null(targets) | !is.list(targets) | !is.data.frame(targets[[1]])) {
    stop("Targets must be a list of data frames.")
  }

  #Initialize
  if(length(targets) == 1 & is.null(seed)) {
    stop("Check_targets() looks for compatibility issues between targets and the seed. If only one target is supplied, you must also supply a seed.
         Otherwise, no need to check. :)")
  } else if(length(targets) == 1) {
    message("Beginning check... one target and seed provided")
  } else {
    message(paste0("Beginning check... ", length(targets), " targets provided",
                   if(is.null(seed)){"."} else {", as well as a seed."}))
  }

  ###
  # Check the targets against each other ----
  ###
  num.tars <- length(targets)
  tar.list <- targets

  names(tar.list) <- paste0("Tar", 1:num.tars)

  if(length(tar.list) <= 1){
    target.checks.op <- list()
  } else {
    tar.combo <- data.frame(TarA = c(), TarB = c(), stringsAsFactors = FALSE)
    for(i in num.tars:2){
      for(j in (i-1):1){
        tar.combo <- tar.combo %>% bind_rows(
          data.frame(TarA = paste0("Tar", j), TarB = paste0("Tar", i), stringsAsFactors = FALSE)
        )
      }
    }
    tar.combo <- tar.combo %>%
      arrange(TarA, TarB)

    combine_tars_a <- function(tara, tarb, list_of_tars){
      TarA <- list_of_tars[[tara]]
      TarB <- list_of_tars[[tarb]]

      common.dims <- names(TarA)[names(TarA) %in% names(TarB)]
      common.dims <- common.dims[!(common.dims %in% target.value.names)]

      names(TarA)[names(TarA) %in% target.value.names] <- ".value"
      names(TarB)[names(TarB) %in% target.value.names] <- ".value"

      combo.tar <- TarA %>%
        mutate(.dftype = "Checker") %>% #This cleans up joins for targets with no common series
        group_by_(.dots = as.list(c(".dftype", common.dims))) %>%
        summarize(TarA = sum(.value, na.rm=TRUE)) %>%
        ungroup() %>%
        full_join(
          TarB %>%
            mutate(.dftype = "Checker")  %>%
            group_by_(.dots = as.list(c(".dftype", common.dims))) %>%
            summarize(TarB = sum(.value, na.rm=TRUE)) %>%
            ungroup() ,
          by = c(".dftype", common.dims)
        ) %>%
        mutate(Check_value = (TarA - TarB),
               Check_trigger = abs(Check_value) > max.error,
               Check_trigger = ifelse(is.na(Check_trigger), FALSE, Check_trigger),
               Mismatch_trigger = is.na(TarA) | is.na(TarB)) %>%
        select(-.dftype)

      names(combo.tar)[names(combo.tar) == "TarA"] <- tara
      names(combo.tar)[names(combo.tar) == "TarB"] <- tarb

      return(combo.tar)
    }

    target.checks <- purrr::pmap(list(a = tar.combo$TarA, b = tar.combo$TarB),
                                 function(a, b){combine_tars_a(a, b, tar.list)})

    target.checks.names <- tar.combo %>% mutate(.name = paste(TarA, " & ", TarB)) %>% pull
    names(target.checks) <- target.checks.names

    #Only keep dfs with violations
    target.checks.op <- target.checks[purrr::map_lgl(target.checks, function(x){any(x$Check_trigger, na.rm=TRUE)})]

    #Look for targets with values in some but not others
    target.mismatch.op <- target.checks[purrr::map_lgl(target.checks, function(x){any(x$Mismatch_trigger, na.rm=TRUE)})]

    if(length(target.checks.op) == 0 && length(target.mismatch.op) == 0) {
      message("\nTargets are good! No issues here.\n===================================")
    } else {
      message("\nAt least one violation has been found within the targets. See output.\n===================================")
    }
  } #End target check


  ####
  # Check the targets against the seed ----
  ####

  if(is.null(seed)){ seed.checks.op <- list()} else {

    message("Checking each target against the seed... This will look for 0 or NA values over seed subtotals and compare them with the targets.\nIf the seed has a 0 subtotal, then the matching target should also be 0 (IPF cannot scale zero to a non-zero).")

    check_seed_a <- function(TarA, SeedA){
      dims.in.tar <- names(TarA)
      dims.in.tar <- dims.in.tar[!(dims.in.tar %in% c(target.value.names))]

      names(SeedA)[names(SeedA) == seed.value.name] <- ".value"
      names(TarA)[names(TarA) %in% target.value.names] <- ".target"

      seed.collapse <- SeedA %>%
        mutate_if(is.factor, as.character) %>%
        group_by_(.dots = as.list(dims.in.tar)) %>%
        summarize(.seed = sum(.value, na.rm = TRUE)) %>%
        ungroup() %>%
        full_join(TarA %>% mutate_if(is.factor, as.character), by = dims.in.tar) %>%
        mutate(.seed = round(.seed, 3),
               .target = round(.target, 3)) %>%
        mutate(Check_trigger = (.seed == 0 | is.na(.seed)) & (.target > 0 | !is.na(.target)),
               Mismatch_trigger = is.na(.seed) | is.na(.target))

    }

    seed.checks <- purrr::map(tar.list, function(x){check_seed_a(x, seed)})

    names(seed.checks) <- paste(names(seed.checks), "& Seed")

    #Only keep dfs with violations
    seed.checks.op <- seed.checks[purrr::map_lgl(seed.checks, function(x){any(x$Check_trigger)})]
    seed.mismatch.op <- seed.checks[purrr::map_lgl(seed.checks, function(x){any(x$Mismatch_trigger)})]

    #Output message
    if(length(seed.checks.op) == 0 && length(seed.mismatch.op) == 0){
      message("\nThe seed and targets line up! No issues here.\n===================================")
    } else {
      message("\n!!! Zero subtotals found in seed where targets have values. IPF will not converge. See output.\n===================================")
    }

  } #End seed check

  check.op <- purrr::map(c(target.checks.op, target.mismatch.op, seed.checks.op), function(x){
    x %>% filter(Check_trigger | Mismatch_trigger) %>% select(-Check_trigger, -Mismatch_trigger)
  })

  return(check.op)

}



#' Title
#'
#' @param datatable
#' @param targets
#' @param datatable.value.name
#' @param target.value.names
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
#' @param growth_targets
#' @param save.tars
#' @param show.messages
#'
#' @return
#' @export
#'
#' @examples
ip_fitting <- function(datatable, targets,
                   datatable.value.name = "value", target.value.names = "value",
                   max.error = 0.01, max.iterations = 25,
                   freeze_cells = NULL, freeze_cells.value.name = "value",
                   freeze_slice = NULL, freeze_slice.value.names = "value",
                   minmax_cells = NULL, minmax_cells.value.names = c("value_min", "value_max"),
                   minmax_slice = NULL, minmax_slice.value.names = c("value_min", "value_max"),
                   minmax.smash.param = 1/3,
                   growth_targets = NULL,
                   save.tars = FALSE, show.messages = TRUE) {

  #Warnings
  if(is.null(targets) | !is.list(targets) | !is.data.frame(targets[[1]])) {
    stop("Targets must be a list of data frames.")
  }

  #Single Target
  if(length(targets) == 1) {
    if(show.messages) {
      message(paste("Initializing IPF...",
                    "Only 1 target has been given. Data will be scaled to single target with minmax conditions."))
    }
    targets[[2]] <- targets[[1]]
  } else {
    message(paste("Initializing IPF...", length(targets), "targets supplied."))
  }

  #Check targets
  target.checker <- vector(length = length(targets)) #Use Targets, not tar.list
  for(i in seq_along(targets)){
    x <- targets[[i]]
    target.checker[i] <- sum(x[, target.value.names], na.rm=T)
  }
  if(length(unique(round(target.checker, 4))) > 1){
    message("Warning: Given targets do not have the same totals. IPF may not converge.")
  }

  #Set initial conditions
  current.error     <- 10^9
  current.iteration <- 1
  minmax_cells.oob <- TRUE
  minmax_slice.oob <- TRUE

  tar.list <- targets

  #New for 0.0.0.9005 - Slice targets now given as single data frame
  #Freeze 2- Ice Slices
  if(!is.null(freeze_slice)){
    if(!is.data.frame(freeze_slice)) {stop("Parameter freeze_slice must be a data frame containing partial targets or subtotals.")}
    if(show.messages) {
      message(paste(nrow(freeze_slice), "Frozen slices (partial) targets supplied."))
    }

    freeze_slice_list <- ip_load_slice_a(freeze_slice, slice.value.name = freeze_slice.value.names)

    tar.list <- c(tar.list, freeze_slice_list)
  } else {
    for(i in seq_along(tar.list)){
      names(tar.list)[i] <- paste0("tar__", i)
    }
  }

  #initialize
  df0 <- datatable
  names(df0)[names(df0) == datatable.value.name] <- "value"

  for(i in seq_along(tar.list)){
    if(names(tar.list)[i] == "" | is.null(names(tar.list))){
      names(tar.list)[i] <- paste0("tar__", i)
    }
  }

  #Freeze 1 - Ice Cells
  if(!is.null(freeze_cells)){
    if(show.messages) {
      message(paste(nrow(freeze_cells), "Frozen cell values will be hit."))
    }

    names(freeze_cells)[names(freeze_cells) == freeze_cells.value.name] <- "frz__c"

    df0 <- df0 %>%
      left_join(freeze_cells, by = names(freeze_cells %>% select(-frz__c))) %>%
      #Iced cells are not unknown, so we don't need to include them in the IPF
      mutate(value = ifelse(is.na(frz__c), value, 0))

    #Since we 0'd out the seed, we should also remove all the iced values from the targets
    tar.list2 <- lapply(seq_along(tar.list), function(i){
      x <- tar.list[[i]]
      names(x)[names(x) %in% c(target.value.names, freeze_slice.value.names)] <- "value"

      nm_tar <- names(tar.list)[i]

      ice_target <- freeze_cells %>%
        group_by_(.dots = as.list(names(x)[!(names(x) == "value")])) %>%
        summarize(iced = sum(frz__c, na.rm=T)) %>%
        ungroup()

      df <- x %>%
        left_join(ice_target, by = names(ice_target %>% select(-iced))) %>%
        mutate(iced = ifelse(is.na(iced), 0, iced)) %>%
        mutate(value = value - iced) %>%
        select(-iced)

      if(any(df$value < 0)){stop("Frozen cell values exceed given targets. Unable to rationalize.")}

      #For Subtotal slices, need to set them all equal to the minimum
      if(grepl("frz__subtl", nm_tar, fixed = TRUE)){
        df <- df %>%
          mutate(value = min(value, na.rm=T))
      }

      return(df)
    })

    names(tar.list2) <- names(tar.list)
    tar.list <- tar.list2
  } #End ice cells

  #Freeze 3 - Slush Cells. Similar to Ice, but cells have a min/max bound, not specific value
  if(!is.null(minmax_cells)) {
    if(show.messages) {
      message(paste(nrow(minmax_cells), "Min/max cell values supplied."))
    }

    names(minmax_cells)[names(minmax_cells) == minmax_cells.value.names[1]] <- "minmax__c_min"
    names(minmax_cells)[names(minmax_cells) == minmax_cells.value.names[2]] <- "minmax__c_max"

    df0 <- df0 %>%
      left_join(minmax_cells, by = names(minmax_cells %>% select(-starts_with("minmax__"))))
  }

  #Freeze 4 - Slush Slices. Partial (or complete?) targets with min/max bound
  if(!is.null(minmax_slice)) {
    if(!is.data.frame(minmax_slice)) {
      stop("Parameter minmax_slice must be a data frame containing partial min/max targets.")
    }

    slush.list <- ip_load_slice_a(minmax_slice, slice.value.name = minmax_slice.value.names,  prefix = "mm")

    if(show.messages) {
      message(paste(nrow(slush.list), "Min/max slice targets given."))
    }

    #Freeze 1 inside Freeze 4 - Ice Cells
    if(!is.null(freeze_cells)){

      #Like the normal targets, need to subtract the frozen values for the min/max
      slush.list <- lapply(slush.list, function(x){
        names(x)[names(x) == minmax_slice.value.names[1]] <- "value_min"
        names(x)[names(x) == minmax_slice.value.names[2]] <- "value_max"

        ice_slush <- freeze_cells %>%
          group_by_(.dots = as.list(names(x)[!(names(x) %in% c("value_min", "value_max"))])) %>%
          summarize(iced = sum(frz__c, na.rm=T)) %>%
          ungroup()

        df <- x %>%
          left_join(ice_slush, by = names(ice_slush %>% select(-iced))) %>%
          mutate(iced = ifelse(is.na(iced), 0, iced)) %>%
          mutate(value_min = value_min - iced,
                 value_max = value_max - iced) %>%
          select(-iced)

        if(any(df$value_min < 0)){stop("Frozen cell values exceed supplied minmax_slice targets. Unable to rationalize.")}

        names(x)[names(x) == "value_min"] <- minmax_slice.value.names[1]
        names(x)[names(x) == "value_max"] <- minmax_slice.value.names[2]

        return(df)
      })
    } #End ice cells

    #Format each slush slice df - give each value unique names
    for( i in seq_along(slush.list)){
      x <- slush.list[[i]]

      names(x)[names(x) == minmax_slice.value.names[1]] <- paste0("minmax__s_min", i)
      names(x)[names(x) == minmax_slice.value.names[2]] <- paste0("minmax__s_max", i)

      names(slush.list[[i]]) <- names(x)
    }

    #Add each of those slice targets to the master data frame
    for(x in slush.list){
      df0 <- df0 %>%
        left_join(x, by = names(x %>% select(-starts_with("minmax__"))))
    }

  } #End Freeze 4 setup

  #Format targets - give each value unique names
  for( i in seq_along(tar.list)){
    x <- tar.list[[i]]

    nm <- names(tar.list)[i]

    #If it's a target, give it the prefix tar__
    if(i <= length(targets)) {
      names(x)[names(x) %in% target.value.names] <- paste0("tar__", i)
      names(tar.list[[i]]) <- names(x)
    } else if(grepl("__slice", nm, fixed=T)) { #Else it's a freeze slice
      names(x)[names(x) %in% target.value.names] <- paste0("tar__frz__slice", i)
      names(tar.list[[i]]) <- names(x)
    } else {
      names(x)[names(x) %in% target.value.names] <- paste0("tar__frz__subtl", i)
      names(tar.list[[i]]) <- names(x)
    }

  }

  #Format input seed
  for(i in seq_along(tar.list)){
    x <- tar.list[[i]]

    df0 <- df0 %>%
      left_join(x, by = names(x %>% select(-dplyr::contains("__"))))

    #If it's a target and has NAs, replace them with 0
    if(i <= length(targets)) {
      df0[, paste0("tar__", i)][is.na(df0[, paste0("tar__", i)])] <- 0
    }
  }

  # Growth targets
  # These targets are for structure. The actual values come later
  if(!is.null(growth_targets)){
    if(show.messages) {
      message(paste(nrow(growth_targets), "Growth rate cell/slice targets supplied."))
    }
  }

  #IPF Loop
  while((current.error > max.error | minmax_cells.oob == TRUE | minmax_slice.oob == TRUE) &
        current.iteration <= max.iterations ) {

    if(show.messages) {
      message(paste("Iteration Summary:", current.iteration))
    }

    #Reset Freeze 3+4 - Slush Cells/Slice - to off
    minmax_cells.oob <- FALSE
    minmax_slice.oob <- FALSE

    ##
    # Scaling
    ##

    df1 <- df0
    for(i in seq_along(tar.list)){
      x <- tar.list[[i]]
      nm <- names(tar.list)[i]

      #If it's a target
      if(i <= length(targets)) {
        df1 <- df1 %>%
          ip_scale_a(target_series = names(x)[!(names(x) %in% c("value"))], series_target = paste0("tar__", i))
      } else if(grepl("__slice", nm, fixed=T)) { #Else it's a freeze slice
        df1 <- df1 %>%
          ip_scale_a(target_series = names(x)[!(names(x) %in% c("value"))],
                     series_target = paste0("tar__frz__slice", i), series_type = "slice")
      } else {
        df1 <- df1 %>%
          ip_scale_a(target_series = names(x)[!(names(x) %in% c("value"))],
                     series_target = paste0("tar__frz__subtl", i), series_type = "subtl")
      }
    }

    ###
    # Growth targets
    ###
    if(!is.null(growth_targets)){
      df1 <- df1[, !grepl("tar__growth__", names(df1), fixed = T)] #Remove previous yy value targets

      #Create list of target dataframes, with growth rates converted to levels
      growth_tars_current_df <- df1

      #Need to add back in frozen cells for Growth calc
      if(!is.null(freeze_cells)){
        growth_tars_current_df <- growth_tars_current_df %>%
          mutate(value = ifelse(is.na(frz__c), value, frz__c))
      }
      growth_tars_current <- growth_tars_current_df %>% ip_growth_transform(growth_targets)

      for(i in seq_along(growth_tars_current)) {
        names(growth_tars_current)[i] <- paste0("tar__growth__", i)
      }

      #Scale each growth value target
      for(i in seq_along(growth_tars_current)){
        x <- growth_tars_current[[i]]

        df1 <- df1 %>%
          left_join(x, by = names(x)[!(names(x) %in% paste0("tar__growth__", i))]) %>%
          ip_scale_a(target_series = names(x)[!(names(x) %in% paste0("tar__growth__", i))],
                     series_target = paste0("tar__growth__", i))
      }

    } else {
      growth_tars_current <- NULL
    }

    ###
    #Error
    ###

    #For each target (except final one), calculate the error over each element in the target,
    # Save results to a list of dfs
    #err.list <- lapply(seq_along(head(tar.list, -1)), function(i){
    err.list <- lapply(seq_along(tar.list), function(i){
      x <- tar.list[[i]]
      nm <- names(tar.list)[i]

      #If it's a target
      if(i <= length(targets) | grepl("__growth", nm, fixed=T)) {
        df <- df1 %>%
          ip_miss_a(names(x)[!(names(x) %in% c("value"))], series_target = paste0("tar__", i))
      } else if(grepl("__slice", nm, fixed=T)) { #Else it's a freeze slice
        df <- df1 %>%
          ip_miss_a(names(x)[!(names(x) %in% c("value"))],
                    series_target = paste0("tar__frz__slice", i), series_type = "slice")
        # } else if(grepl("__subtl", nm, fixed=T)){
      } else {
        df <- df1 %>%
          ip_miss_a(names(x)[!(names(x) %in% c("value"))],
                    series_target = paste0("tar__frz__subtl", i), series_type = "subtl")
      }
    })

    #Growth error
    if(!is.null(growth_targets)){

      err.list.growth <- lapply(seq_along(growth_tars_current), function(i){
        x <- growth_tars_current[[i]]
        nm <- names(growth_tars_current)[i]

        df <- df1 %>%
          ip_miss_a(names(x)[!(names(x) %in% c("value"))], series_target = paste0("tar__growth__", i))
      })

      err.list <- c(err.list, err.list.growth)

    }

    #For each error df, sum the error column, then sum those sums for total abs error
    current.error <- sum(sapply(err.list, function(x){
      return(sum(abs(x$error), na.rm=T))
    }), na.rm=T)

    if(show.messages) {
      message(paste("    Iteration Error: ", round(current.error, 5)))
    }

    ###
    # Freeze 4 - Min/Max Slice
    ###
    #Similar to Ice, but cells have a min/max bound, not specific value
    #Calculate AFTER error, as this is its own deal breaker. If everything is in bounds, no problem
    if(!is.null(minmax_slice)) {

      #1) For each minmax_slice target, check to see there are any OOB
      for(i in seq_along(slush.list)) {
        x <- slush.list[[i]]

        if(grepl("__subtl", names(slush.list)[i], fixed = T)){
          type <- "subtl"
        } else{
          type <- "slice"
        }

        #Roll-up version of df1 to check slush targets. Use generic target names for easier programming
        names(df1)[names(df1) == paste0("minmax__s_min", i)] <- "minmax__s_min"
        names(df1)[names(df1) == paste0("minmax__s_max", i)] <- "minmax__s_max"

        df1_slushs <- df1
        names(x)[names(x) == paste0("minmax__s_min", i)] <- "minmax__s_min"
        names(x)[names(x) == paste0("minmax__s_max", i)] <- "minmax__s_max"

        df1_slushs <- df1_slushs %>%
          #Group by names of the target, INCLUDING the target value names, since they are the same for each group element
          do(
            if(type == "subtl"){
              group_by_(., .dots = as.list(c("minmax__s_min", "minmax__s_max")))
            } else {
              group_by_(., .dots = as.list(names(x)))
            }
          ) %>%
          summarize(value = sum(value, na.rm = T)) %>%
          ungroup() %>%
          mutate(check_slush_min = value < minmax__s_min,
                 check_slush_max = value > minmax__s_max) %>%
          rowwise() %>%
          mutate(check_slush_s = any(c(check_slush_max, check_slush_min), na.rm=T)) %>% #Look for OOB
          ungroup() %>%
          select(-check_slush_min, -check_slush_max)

        #If any of the bounded cells are OOB,
        # Then proceed to edit the df1 values to smash them into bounds
        # Then adjust the other values in the opposite direction

        if(any(df1_slushs$check_slush_s, na.rm = T)){
          df1_slushs$check_slush_s <- NULL

          if(show.messages) {
            message(paste("    Out of bounds conditions present for minmax_slice: Target", i))
          }

          #First update the trigger to True - There needs to be another iterations
          minmax_slice.oob <- TRUE

          #Then figure out by how much to move the values
          df1_slushs <- df1_slushs %>%
            #Allow for no lower bound
            mutate(minmax__s_min = ifelse(is.na(minmax__s_min) & !is.na(minmax__s_max), 0 , minmax__s_min)) %>%
            mutate(
              value_o_slush = ifelse(is.na(minmax__s_min) & is.na(minmax__s_max), NA, ifelse((value >= minmax__s_min) | is.na(minmax__s_min), value,
                                                                                             #Allow for no upper bound
                                                                                             ifelse(!is.na(minmax__s_max), minmax__s_min + (minmax.smash.param)*(minmax__s_max - minmax__s_min),
                                                                                                    (1+minmax.smash.param)*minmax__s_min))),
              value_o_slush = ifelse((value_o_slush <= minmax__s_max) | is.na(minmax__s_max), value_o_slush, minmax__s_max - (minmax.smash.param)*(minmax__s_max - minmax__s_min))
            ) %>%
            #How much the remaining values have to move to compensate for that shift
            mutate(value_o_remainder = ifelse(is.na(value_o_slush), value, NA),
                   ratio_remainder = (sum(value, na.rm = T) - sum(value_o_slush, na.rm=T))/sum(value_o_remainder, na.rm=T)) %>%
            #Combine into a single ratio to apply to data cube
            mutate(ratio = ifelse(is.na(value_o_slush), ratio_remainder, value_o_slush/value)) %>%
            do(
              if(type == "subtl"){
                select(., -value, -value_o_slush, -value_o_remainder, -ratio_remainder)
              } else {
                select(., -value, -minmax__s_min, -minmax__s_max, -value_o_slush, -value_o_remainder, -ratio_remainder)
              }
            )

          #Apply this df of ratios to our main data frame
          df1 <- df1 %>%
            left_join(df1_slushs, by = names(df1_slushs %>% select(-ratio))) %>%
            mutate(ratio = ifelse(is.na(ratio), 1, ratio)) %>%
            mutate(value = value * ratio) %>%
            select(-ratio)

        }

        names(df1)[names(df1) == "minmax__s_min"] <- paste0("minmax__s_min", i)
        names(df1)[names(df1) == "minmax__s_max"] <- paste0("minmax__s_max", i)

        df1_slushs <- NULL
      }

    } #End Freeze 4

    ###
    # Freeze 3 - Min/Max Cells
    ###
    #Similar to Ice, but cells have a min/max bound, not specific value
    #Calculate AFTER error, as this is its own deal breaker. If everything is in bounds, no problem
    #Do Freeze 3 AFTER Freeze 4 because this is more restrictive

    if(!is.null(minmax_cells)) {

      df1 <- df1 %>%
        mutate(check_slush_min = value < minmax__c_min,
               check_slush_max = value > minmax__c_max) %>%
        rowwise() %>%
        mutate(check_slush_c = any(c(check_slush_max, check_slush_min), na.rm=T)) %>% #Look for OOB
        ungroup() %>%
        select(-check_slush_min, -check_slush_max)

      #If any of the bounded cells are OOB, first update the trigger to True
      minmax_cells.oob <- any(df1$check_slush_c, na.rm = T)
      df1$check_slush_c <- NULL

      #Then replace those values with 1/3 closer above/below that missed bound
      if(minmax_cells.oob == TRUE) {
        if(show.messages) {
          message("    Out of bounds conditions present for minmax_cells.")
        }

        df1 <- df1 %>%
          #Allow for no lower bound
          mutate(minmax__c_min = ifelse(is.na(minmax__c_min) & !is.na(minmax__c_max), 0 , minmax__c_min)) %>%
          mutate(value = ifelse((value >= minmax__c_min) | is.na(minmax__c_min), value,
                                #Allow for no upper bound
                                ifelse(!is.na(minmax__c_max), minmax__c_min + (minmax.smash.param)*(minmax__c_max - minmax__c_min),
                                       (1+minmax.smash.param)*minmax__c_min)),
                 value = ifelse((value <= minmax__c_max) | is.na(minmax__c_max), value, minmax__c_max - (minmax.smash.param)*(minmax__c_max - minmax__c_min))
          )
      }
    } #End Freeze 3


    ###
    # Prepare for next loop
    ###

    current.iteration <- current.iteration + 1

    df0 <- df1

  }

  #Freeze 1 - Iced Cells
  #Add back Iced Cells
  if(!is.null(freeze_cells)){
    df0 <- df0 %>%
      mutate(value = ifelse(is.na(frz__c), value, frz__c))
  }

  #Save/print all the assumption columns?
  # Currently, tars are reduced by Freeze tars.
  # TODO: Add frozen values back into tars... or just replace tars. That's easier
  if(!save.tars){
    df0 <- df0 %>%
      select(-dplyr::contains("__"))
  }

  return(df0)

}

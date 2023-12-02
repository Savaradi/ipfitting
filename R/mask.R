
#' Title
#'
#' @param mask_raw
#' @param mask.value.name
#' @param drop.names
#'
#' @return
#' @export
#'
#' @examples
ip_load_mask_a <- function(mask_raw, mask.value.name = "mask",
                           drop.names = c("Notes")) {
  mask <- mask_raw
  df <- mask_raw
  df <- df[, names(df)[!(names(df) %in% c(mask.value.name, drop.names))]]

  mask$tar__names <- lapply(1:nrow(df), function(i){
    x <- df[i, ]
    dat <- names(df)[!is.na(x) & !(names(df) %in% mask.value.name)]

    return(paste(c(dat), collapse = " "))
  })

  mask_list <- lapply(unlist(unique(mask$tar__names)), function(x){
    dat <- mask %>%
      filter(tar__names == x) %>%
      select(-tar__names)

    i <- which(unlist(unique(mask$tar__names)) == x)
    names(dat)[names(dat) == mask.value.name] <- paste0("mask_", i)

    dat <- dat[, colSums(is.na(dat)) < nrow(dat) | names(dat) %in% mask.value.name]
    dat <- dat[, !(names(dat) %in% c(drop.names))]

    return(dat)

  })

  for(i in seq_along(mask_list)){
    names(mask_list)[i] <- paste0("mask", "__", i)
  }

  return(mask_list)
}

#' Apply a layer of multiplicative masks to a seed.
#'
#' @param seed
#' @param mask
#' @param mask.value.name
#' @param drop.names


ip_mask_seed <- function(seed, mask, mask.value.name = "mask",
                         drop.names = c("Notes")){

  seed_mask <- ip_load_mask_a(mask, mask.value.name = mask.value.name,
                              drop.names = drop.names)

  seed_list <- vector("list", length = (length(seed_mask) + 1))
  seed_list[[1]] <- seed
  seed_list[2:length(seed_list)] <- seed_mask

  seed_list %>%
    reduce(left_join) %>%
    mutate_at(vars(dplyr::starts_with("mask_")), funs(ifelse(is.na(.), 1, .))) %>%
    gather(value__type, value__mask, value, dplyr::starts_with("mask_")) %>%
    group_by_(.dots = as.list(names(.)[!(names(.) %in% c("value__type", "value__mask"))])) %>%
    summarize(value = prod(value__mask, na.rm=TRUE)) %>%
    ungroup()

}


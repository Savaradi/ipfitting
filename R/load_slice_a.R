
#' Title
#'
#' @param slice
#' @param slice.value.name
#' @param drop.names
#' @param prefix
#'
#' @return
#' @export
#'
#' @examples
ip_load_slice_a <- function(slice, slice.value.name = "value", drop.names = c("Notes"), prefix = "frz") {

  df <- slice
  df <- df[, names(df)[!(names(df) %in% slice.value.name)]]

  slice$tar__names <- lapply(1:nrow(df), function(i){
    x <- df[i, ]
    dat <- names(df)[!is.na(x) & !(names(df) %in% slice.value.name)]

    if(any(grepl(" + ", x, fixed = T))) { subtotal.checker <- paste0("SUBTOTAL_PRESENT", i)} else {subtotal.checker <- NULL}

    return(paste(c(dat, subtotal.checker), collapse = " "))
  })

  slice_list <- lapply(unique(slice$tar__names), function(x){
    dat <- slice %>%
      filter(tar__names == x) %>%
      select(-tar__names)

    dat <- dat[, colSums(is.na(dat))<nrow(dat) | names(dat) %in% slice.value.name]

    if(grepl("SUBTOTAL_PRESENT", x, fixed = T)) {
      daty <- as.data.frame(expand.grid(strsplit(as.character(dat[1, ]), " + ", fixed = T), stringsAsFactors = F))
      names(daty) = names(dat)
      daty[, slice.value.name[1]] <- as.numeric(daty[, slice.value.name[1]])
      if(length(slice.value.name)==2){
        daty[, slice.value.name[2]] <- as.numeric(daty[, slice.value.name[2]])
      }
      dat <- as.data.frame(daty)
    }

    dat <- dat[, colSums(is.na(dat))<nrow(dat) | names(dat) %in% slice.value.name]

    return(dat)

  })

  for(i in seq_along(slice_list)){
    if(grepl("SUBTOTAL_PRESENT", unique(slice$tar__names)[i], fixed = T)) {
      names(slice_list)[i] <- paste0(prefix, "__subtl_", i)
    } else {
      names(slice_list)[i] <- paste0(prefix, "__slice_", i)
    }

  }

  return(slice_list)
}

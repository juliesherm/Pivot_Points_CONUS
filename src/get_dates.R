## ---------------------------
##
## Script name: get_dates.R
##
## Purpose of script: gets dates from either water balance, gridmet, or modis data
##
## Author: Julie Sherman
##
## ---------------------------
##
## Notes:
##  - function that extracts dates from names of spatraster depending on source
##  - valid sources are "MODIS", "GRIDMET", "WB"
##
## ---------------------------

get_dates <- function(source, stk) {
  if (source == "WB") {
    dates <- lapply(stk, time)

    # if wb_stk has no time attribute, get date from file name
    if (any(unlist(lapply(dates, is.na)))) {
      dates <- lapply(stk, function(var_stk) {
        names <- sources(var_stk)
        names <- paste(str_split_i(names, "_", 5), str_split_i(names, "_", 13), "01", sep = "-")
        return(as.Date(names))
      })
    }
  } else if (source == "GRIDMET") {
    gm_vars <- names(stk)
    dates <- sapply(gm_vars, function(x) {
      list(as.Date(
        as.numeric(str_sub(
          names(stk[[x]]),
          -5, -1
        )),
        origin = as.Date("1900-01-01")
      ))
    })
  } else if (source == "MODIS") {
    #vi_vars <- names(stk)

    #dates <- sapply(vi_vars, function(x) {
    #  m_dates <- vi_stk[[x]]@cpp[["names"]] %>%
    #    str_split(., "_", simplify = TRUE) %>%
    #    as.data.frame() %>%
    #    select(V3, V4) %>%
    #    mutate_all(as.numeric) %>%
    #    rename(year = V3, doy = V4)

    #  return(list(as.Date(m_dates$doy, origin = as.Date(paste0(m_dates$year, "-01-01")) - 1)))
    #})

    #names(dates) <- vi_vars

    fnames <- list.files(here("raw_data", "MODIS", "CONUS", "VI_16Days_250m_v61", "EVI"),
                          pattern = "ydays",
                          full.names = TRUE)
    
    dates <- sapply(fnames, function(fname) {
      yr <- str_sub(basename(fname), 6, 9)
      ydays <- as.vector(read.csv(fname)[[1]])

      days <- as.Date(ydays, origin = as.Date(paste0(yr,"-01-01")))
    })

    dates <- unlist(dates)
  } else {
    stop("unknown date source")
  }

  if (length(dates) > 1) {
    # check if the dates are the same across all variables
    all.identical.list <- function(l) identical(unname(l[-length(l)]), unname(l[-1]))

    if (all.identical.list(dates)) {
      dates <- dates[[1]]
      print("all variables have  the same dates. returning single vector of  dates")
    } else {
      warning("variables have unequal dates! returning list of dates for each")
    }
  }

  return(dates)
}

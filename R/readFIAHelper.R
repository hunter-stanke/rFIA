readFIAHelper1 <- function(x, dir, ...){
  # Read in and append each file to a list
  file <- data.table::fread(paste(dir, x, sep = ""), showProgress = FALSE, logical01 = FALSE, integer64 = 'numeric', ...)
  # suppressMessages({suppressWarnings({
  # f <- read_csv(paste(dir, x, sep = ""), progress = FALSE, col_type = cols(.default = col_character()), ...)
  # f <- type_convert(f)
  # })})
  file <- as.data.frame(file)

  file
}

getFIAHelper <- function(x, dir, ...){
  # Download and append each file to a list
  file <- data.table::fread(x, showProgress = FALSE, logical01 = FALSE, integer64 = 'numeric', ...)

  # Write the data out the directory they've chosen
  if(!is.null(dir)){
    data.table::fwrite(x = file, file = paste0(dir, str_sub(x, 43, -1)), showProgress = FALSE)
  }

  file <- as.data.frame(file)

  file
}


readFIAHelper2 <- function(x, tables){
  subList <- tables[str_sub(x, 4) == unique(str_sub(x,4))]
  mergedTable <- bind_rows(subList, .id = NULL)

  gc()
  mergedTable
}


#   if (anyDuplicated(str_sub(x, 4)) != 0){
#
#     #name <- unique(str_sub(x, 4, -5))
#
#     for (i in 1:length(unique(str_sub(tableNames, 4)))){
#       subList <- fileList[str_sub(tableNames, 4) == unique(str_sub(tableNames,4))[i]]
#       name <- unique(str_sub(tableNames, 4, -5))[i]
#       # Give a ton of warnings about factors and characters, don't do that
#       outList[[name]] <- bind_rows(subList, .id = NULL)
#     }
#   } else {
#     outList <- fileList
#     names(outList) <- unique(str_sub(tableNames, 4, -5))
#   }
# }

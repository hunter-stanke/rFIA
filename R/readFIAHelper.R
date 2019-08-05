readFIAHelper1 <- function(x, dir, ...){
  # Read in and append each file to a list
  file <- read.csv(paste(dir, x, sep = ""), ...)

  gc()
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

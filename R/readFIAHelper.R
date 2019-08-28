readFIAHelper1 <- function(x, dir, ...){
  # Read in and append each file to a list
  file <- fread(paste(dir, x, sep = ""), showProgress = FALSE, logical01 = FALSE, integer64 = 'double', ...)

  # We don't want data.table formats
  file <- as.data.frame(file)

  ## Use a loop to avoid using bit64 package (integer64 columns still appear for some reason)
  classes <- sapply(file, class)
  for (i in 1:ncol(file)){
    if (classes[i] == 'integer64'){
      file[,i] <- as.double(file[,i])
    }
  }


  file
}

getFIAHelper <- function(x, dir, ...){
  # Download and append each file to a list
  file <- fread(x, showProgress = FALSE, logical01 = FALSE, integer64 = 'double', ...)

  # Write the data out the directory they've chosen
  if(!is.null(dir)){
    fwrite(x = file, file = paste0(dir, str_sub(x, 43, -1)), showProgress = FALSE)
  }

  # We don't want data.table formats
  file <- as.data.frame(file)

  ## Use a loop to avoid using bit64 package (integer64 columns still appear for some reason)
  classes <- sapply(file, class)
  for (i in 1:ncol(file)){
    if (classes[i] == 'integer64'){
      file[,i] <- as.double(file[,i])
    }
  }

  file
}


# writeFIAHelper <- function(x, db, dir, ...){
#   # Write the data out the directory they've chosen
#   fwrite(x = db[[x]], file = paste0(x, '.csv'), showProgress = FALSE)
# }



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

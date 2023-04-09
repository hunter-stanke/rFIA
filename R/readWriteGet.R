
## Sometimes these read as dates, and that breaks readFIA when mutliple states are read
## Maybe also think about dropping uncommon columns from TREE. Could save a lot of memory
dropTheseCols <- function() {

  dropThese <- c("CREATED_BY", "CREATED_IN_INSTANCE", "CREATED_DATE", "MODIFIED_BY", "MODIFIED_IN_INSTANCE", "MODIFIED_DATE")
  dropThese <- c(dropThese, tolower(dropThese))
  return(dropThese)
}

# Read in FIA database files (.csv) from local directory
#' @export
readFIA <- function(dir = NULL,
                    con = NULL,
                    schema = NULL,
                    common = TRUE,
                    tables = NULL,
                    states = NULL,
                    inMemory = TRUE,
                    nCores = 1,
                    ...){

  ## Must specify dir or con
  if (is.null(dir) & is.null(con)) {
    stop('Must specify `dir` or `con`. If reading from csv files, set `dir` to the directory where your FIA data lives.')
  }

  ## Need schema if reading from a database
  if (!is.null(con) & is.null(schema)) {
    stop('Must specify `schema` when reading data from a database connection (`con`).')
  }

  ## Attempt to read from database
  if (!is.null(con)) {
    db <- readFIA.db(con, schema, states, tables, common, inMemory)
    return(db)
  }


  # Add a slash to end of directory name if missing
  if (str_sub(dir,-1) != '/') dir <- paste(dir, '/', sep = "")
  # Grab all the file names in directory
  files <- list.files(dir)

  ## methods for reading the full database into memory
  if (inMemory){

    inTables <- list()

    # Some warnings
    if(!dir.exists(dir)) {
      stop(paste('Directory', dir, 'does not exist.'))
    }
    if(length(files[str_sub(files,-4, -1) == '.csv']) < 1){
      stop(paste('Directory', dir, 'contains no .csv files.'))
    }

    ## Some warnings up front
    ## Do not try to merge ENTIRE with other states
    if (length(states) > 1 & any(str_detect(str_to_upper(states), 'ENTIRE'))){
      stop('Cannot merge ENITRE with other state tables. ENTIRE includes all state tables combined. Do you only need data for a particular region?')
    }

    # Only csvs
    files <- files[str_sub(files,-4,-1) == '.csv']

    ## All states present in the directory, adding reference if present
    allStates <- unique(str_to_upper(str_sub(files, 1, 2)))
    allStates <- allStates[allStates %in% intData$stateNames$STATEAB]
    if (any(str_to_upper(str_sub(files, 1, 3)) == 'REF')) allStates <- c(allStates, 'REF')

    ## If states are given, then subset files now,
    ## otherwise, bring them all in
    if (!is.null(states)){

      ### I'm not very smart and like specify the name twice sometimes,
      ### --> making the function smarter than me
      states <- str_to_upper(unique(states))

      ## Make sure the states exist to begin with
      if (any(states %in% allStates == FALSE)){
        missStates <- states[states %in% allStates == FALSE]
        stop(paste('Data unavailable for: ', paste(as.character(missStates),collapse = ', '), '. States not found in specified directory.'))
      }

      state.files <- files[str_to_upper(str_sub(files, 1, 2)) %in% states]
      if ('REF' %in% states) {
        ref.files <- files[str_to_upper(str_sub(files, 1,3)) == 'REF']
        files <- c(state.files, ref.files)
      } else {
        files <- c(state.files)
      }

    } else {
      ## Checking if state files and merged state files are mixed in the directory.
      states <- unique(str_to_upper(str_sub(files, 1, 3)))
      ## Dropping lichen spp summary and project reference tables
      trueStates <- states[str_sub(states, 3,3) == '_']
      if ("REF" %in% states) {
        trueStates <- c(trueStates, 'REF')
      }

      ## If length is zero, then all merged states - great
      ## If length states is the same as true States, then all state files - great
      ## Otherwise, they're probably mixed. Throw a warning and read only the states
      if (length(trueStates) > 0 & length(trueStates) < length(states)) {
        warning('Found data from merged states and individual states in same directory. Reading only individual states files.')
        state.files <- files[str_to_upper(str_sub(files, 3, 3)) == '_']
        ref.files <- files[str_to_upper(str_sub(files, 1,3)) == 'REF']
        files <- c(state.files, ref.files)
      }
    }

    ## Common file names to consider
    ## All reference tables should be common
    cFiles <- c('COND', 'COND_DWM_CALC', 'INVASIVE_SUBPLOT_SPP', 'PLOT', 'POP_ESTN_UNIT',
                'POP_EVAL', 'POP_EVAL_GRP', 'POP_EVAL_TYP', 'POP_PLOT_STRATUM_ASSGN', 'POP_STRATUM',
                'SUBPLOT', 'TREE', 'TREE_GRM_COMPONENT', 'TREE_GRM_MIDPT', 'TREE_GRM_BEGIN', 'SUBP_COND_CHNG_MTRX',
                'SEEDLING', 'SURVEY', 'SUBP_COND', 'P2VEG_SUBP_STRUCTURE')

    if ('REF' %in% allStates) {
      cFiles <- c(cFiles,
                  "CITATION", "DIFFERENCE_TEST_PER_ACRE", "DIFFERENCE_TEST_TOTALS", "FIADB_VERSION",
                  "FOREST_TYPE_GROUP", "FOREST_TYPE", "GRM_TYPE", "HABTYP_DESCRIPTION", "HABTYP_PUBLICATION",
                  "INVASIVE_SPECIES", "LICHEN_SPECIES", "LICHEN_SPP_COMMENTS", "NVCS_HIERARCHY_STRCT",
                  "NVCS_LEVEL_1_CODES", "NVCS_LEVEL_2_CODES", "NVCS_LEVEL_3_CODES", "NVCS_LEVEL_4_CODES",
                  "NVCS_LEVEL_5_CODES", "NVCS_LEVEL_6_CODES", "NVCS_LEVEL_7_CODES", "NVCS_LEVEL_8_CODES",
                  "OWNGRPCD", "PLANT_DICTIONARY", "POP_ATTRIBUTE", "POP_EVAL_TYP_DESCR", "RESEARCH_STATION",
                  "SPECIES_GROUP", "SPECIES", "STATE_ELEV", "UNIT")
    }

    ## Drop state/ref abbreviations and .csv
    abbs <- c(paste0('start_here_', intData$stateNames$STATEAB, '_'), 'start_here_REF_')
    file.sub <- paste0('start_here_', str_remove(files, '.csv'))
    ## only allow subset at the start at the beginning of string
    for ( i in abbs ) file.sub <- str_remove(file.sub, i)


    # Only read in the specified tables
    if (!is.null(tables)){
      files <- files[file.sub %in% tables]

    # Otherwise if common=TRUE, grab all common
    } else if (common) {
      files <- files[file.sub %in% cFiles]

    # Otherwise only known FIA tables
    } else {
      files <- files[file.sub %in% intData$fiaTableNames]
    }



    ## Now read them in
    inTables <- list()
    for (n in 1:length(files)){
      ## If MODIFIED_DATE is not present, will warn
      suppressWarnings({
        # Read in and append each file to a list
        file <- data.table::fread(paste(dir, files[n], sep = ""), showProgress = FALSE,
                                  integer64 = 'double', logical01 = FALSE, nThread = nCores,
                                  drop = dropTheseCols(), ...)
      })

      # We don't want data.table formats
      #file <- as.data.frame(file)
      fileName <- str_sub(files[n], 1, -5)

      # Skip over files that are empty
      if(nrow(file) > 0){
        inTables[[fileName]] <- file
      }
    }

    # Give them some names
    #names(inTables) <- files
    #inTables <- lapply(inTables, as.data.frame)

    outTables <- list()
    if (any(str_sub(names(inTables), 3, 3) == '_')){ ## STATE NAMING CONVENTION
      # Remove the state prefix
      names(inTables) <- str_sub(names(inTables), 4, -1)
    }
    uniqueNames <- unique(names(inTables))
    ## Works regardless of whether or not there are duplicate names (multiple states)
    for (i in 1:length(uniqueNames)){
      outTables[[uniqueNames[i]]] <-  data.table::rbindlist(inTables[names(inTables) == uniqueNames[i]], fill = TRUE)
    }


    # NEW CLASS NAME FOR FIA DATABASE OBJECTS
    out <- lapply(outTables, as.data.frame)

    ## Add REF to the reference tables, will have _ prefix
    names(out)[str_sub(names(out), 1, 1) == '_'] <- paste0('REF', names(out)[str_sub(names(out), 1, 1) == '_'])


    ### Methods for keeping data remote until they are needed
    ### Chunking up data into states
    ## inMemory = FALSE
  } else {

    ## Can't use remote methods w/ full database
    if (any(str_detect(str_to_upper(states), 'ENTIRE'))) {
      stop('Remote methods undefind for ENTIRE. Please use `inMemory=TRUE`. Resulting FIA.Database will be quite large (~50GB).')
    }

    ## IF states isn't given, default to all
    ## states in the directory
    if (is.null(states)){
      states <- unique(str_to_upper(str_sub(files, 1, 3)))
      states <- states[str_sub(states, 3,3) == '_']
      states <- str_sub(states, 1, 2)
      ## Only states where abbreviations make sense
      states <- states[states %in% intData$stateNames$STATEAB]
      ## Don't fail if states have been merged
      if (length(states) < 1) states <- 1
    }


    ## Saving the call to readFIA, for eval later
    out <- list(dir = dir,
                common = common,
                tables = tables,
                states = states,
                ... = ...)

    class(out) <- 'Remote.FIA.Database'
  }


  if (inMemory) {

    ## Names of all tables forced to uppercase, since FIA changed this in Oct 2022
    out <- lapply(out, FUN = function(x){dplyr::rename_with(x, toupper)})

    ## If WY updates the END_INVYR for 2018 & 2019 (i.e., not set it as 2020)
    ## then we can drop this. Until then, we have to strong arm END_INVYR
    out <- handleWY(out)

    class(out) <- 'FIA.Database'
  }

  return(out)

}




## Access FIA Database files from the FIA Datamart
#' @export
getFIA <- function(states,
                   dir = NULL,
                   common = TRUE,
                   tables = NULL,
                   load = TRUE,
                   nCores = 1){

  if (!is.null(dir)){
    # Add a slash to end of directory name if missing
    if (str_sub(dir,-1) != '/'){
      dir <- paste(dir, '/', sep = "")
    }
    # Check to see directory exists, if not, make it
    if(!dir.exists(dir)) {
      dir.create(dir)
      message(paste('Creating directory:', dir))
    }

    message(paste0('Saving to ', dir, '. NOTE: modifying FIA tables in Excel may corrupt csv files.'))
  }

  if (is.null(dir) & load == FALSE){
    stop('Must specify a directory ("dir") to save data when "load = FALSE".')
  }

  if ('REF' %in% states & length(unique(states)) > 1) {
    stop('Please download reference tables and state subsets seperately.')
  }

  ## All or nothing w/ new FIADB for some reason
  if (states == 'REF') {
    tables <- NULL
  }

  ## If dir is not specified, hold in a temporary directory
  #if (is.null(dir)){tempDir <- tempdir()}
  tempDir <- tempdir()
  ## Some warnings up front
  ## Do not try to merge ENTIRE with other states
  if (length(states) > 1 & any(str_detect(str_to_upper(states), 'ENTIRE'))){
    stop('Cannot merge ENITRE with other state tables. ENTIRE includes all state tables combined. Do you only need data for a particular region?')
  }
  ## Check to make sure states exist
  allStates <- c('AL', 'AK', 'AZ', 'AR', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'HI', 'ID',
                 'IL', 'IN', 'IA', 'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 'MI', 'MN', 'MS',
                 'MO', 'MT', 'NE', 'NV', 'NH', 'NJ', 'NM', 'NY', 'NC', 'ND', 'OH', 'OK',
                 'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 'WA', 'WV',
                 'WI', 'WY', 'AS', 'FM', 'GU', 'MP', 'PW', 'PR', 'VI', 'ENTIRE', 'REF')
  if (any(str_to_upper(states) %in% allStates == FALSE)){
    missStates <- states[str_to_upper(states) %in% allStates == FALSE]
    stop(paste('Data unavailable for: ', paste(as.character(missStates),collapse = ', '), '. Did you use state/terrority abbreviations? e.g. use states = "AL" (correct) instead of states = "ALABAMA".'))
  }

  ## Check to make sure tables exist
  allTables <- c("BOUNDARY", "COND_DWM_CALC", "COND","COUNTY","DWM_COARSE_WOODY_DEBRIS",
                 "DWM_DUFF_LITTER_FUEL","DWM_FINE_WOODY_DEBRIS","DWM_MICROPLOT_FUEL",
                 "DWM_RESIDUAL_PILE", "DWM_TRANSECT_SEGMENT", "DWM_VISIT","GRND_CVR",
                 "INVASIVE_SUBPLOT_SPP","LICHEN_LAB","LICHEN_PLOT_SUMMARY","LICHEN_VISIT",
                 "OZONE_BIOSITE_SUMMARY","OZONE_PLOT_SUMMARY","OZONE_PLOT","OZONE_SPECIES_SUMMARY",
                 "OZONE_VALIDATION","OZONE_VISIT", "P2VEG_SUBP_STRUCTURE","P2VEG_SUBPLOT_SPP",
                 "PLOT_REGEN","PLOT", "PLOTGEOM", "PLOTSNAP","POP_ESTN_UNIT","POP_EVAL_ATTRIBUTE",
                 "POP_EVAL_GRP","POP_EVAL_TYP","POP_EVAL","POP_PLOT_STRATUM_ASSGN","POP_STRATUM",
                 "SEEDLING_REGEN","SEEDLING","SITETREE","SOILS_EROSION","SOILS_LAB","SOILS_SAMPLE_LOC" ,
                 "SOILS_VISIT", "SUBP_COND_CHNG_MTRX","SUBP_COND","SUBPLOT_REGEN","SUBPLOT",
                 "SURVEY","TREE_GRM_BEGIN","TREE_GRM_COMPONENT","TREE_GRM_ESTN", "TREE_GRM_MIDPT",
                 "TREE_REGIONAL_BIOMASS", "TREE_WOODLAND_STEMS","TREE","VEG_PLOT_SPECIES",
                 "VEG_QUADRAT","VEG_SUBPLOT_SPP","VEG_SUBPLOT", "VEG_VISIT",
                 'CITATION', 'DIFFERENCE_TEST_PER_ACRE', 'DIFFERENCE_TEST_TOTALS',
                 'FIADB_VERSION', 'FOREST_TYPE', 'FOREST_TYPE_GROUP',
                 'GRM_TYPE', 'HABTYP_DESCRIPTION', 'HABTYP_PUBLICATION',
                 'INVASIVE_SPECIES', 'LICHEN_SPECIES', 'LICHEN_SPP_COMMENTS',
                 'NVCS_HEIRARCHY_STRCT', 'NVCS_LEVEL_1_CODES', 'NVCS_LEVEL_2_CODES',
                 'NVCS_LEVEL_3_CODES', 'NVCS_LEVEL_4_CODES', 'NVCS_LEVEL_5_CODES',
                 'NVCS_LEVEL_6_CODES', 'NVCS_LEVEL_7_CODES', 'NVCS_LEVEL_8_CODES',
                 'OWNGRPCD', 'PLANT_DICTIONARY', 'POP_ATTRIBUTE', 'POP_EVAL_TYP_DESCR',
                 'RESEARCH_STATION', 'SPECIES', 'SPECIES_GROUP', 'STATE_ELEV', 'UNIT')
  if (any(str_to_upper(tables) %in% allTables == FALSE)){
    missTables <- tables[str_to_upper(tables) %in% allTables == FALSE]
    stop(paste('Tables: ', paste(as.character(missTables),collapse = ', '), ' unavailble. Check the FIA Datamart at https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html for a list of available tables for each state. Alternatively, specify common = TRUE to download the most commonly used tables.

Did you accidentally include the state abbreviation in front of the table name? e.g. tables = "AL_PLOT" (wrong) instead of tables = "PLOT" (correct).'))
  }

  ### I'm not very smart and like specify the name twice sometimes,
  states <- unique(states)

  ## If individual tables are specified, then just grab those .csvs, otherwise download the .zip file, extract and read with fread. Should be quite a bit quicker.
  if (!is.null(tables)){
    ## Make a list of tables names to read in
    ## Append table names with state abbs and then add url link
    tables <- str_to_upper(tables)

    # Make sure state Abbs are in right format
    states <- str_to_upper(states)
    if ('ENTIRE' %in% states == FALSE) {
      states <- paste0(states, '_')
    } else {
      states <- ''
    }


    tblNames <- sapply(states, FUN = function(x, y){paste0(x,y,'.zip')}, y = tables)
    tblNames <- c(tblNames)
    urls <- paste0('https://apps.fs.usda.gov/fia/datamart/CSV/', tblNames)


    inTables = list()
    for (n in 1:length(urls)){

      newName <- paste0(str_sub(tblNames[n], 1, -5), '.csv')

      ## Download the zip to a temporary file
      temp <- tempfile()
      download.file(urls[n], temp, timeout = 3600)

      # Write the data out the directory they've chosen
      if(is.null(dir)){
        unzip(temp, exdir = tempDir)
        file <- fread(paste0(tempDir, '/', newName), showProgress = FALSE, logical01 = FALSE, integer64 = 'double', nThread = nCores)
        ## Unlinking the directory is bad news, so just delete the file from tempdir
        file.remove(paste0(tempDir, '/', newName))
      } else {
        #download.file(urls[n], paste0(dir, tblNames[n]))
        unzip(temp, exdir = str_sub(dir, 1, -2))
        if (load) file <- fread(paste0(dir, newName), showProgress = FALSE, logical01 = FALSE, integer64 = 'double', nThread = nCores)
      }

      unlink(temp)

      if (load){
        # We don't want data.table formats
        file <- as.data.frame(file)
        inTables[[str_sub(urls[n], 43, -5)]] <- file
      }

    }

    if (load){
      # Check for corresponding tables (multiple states)
      # If they exist, loop through and merge corresponding tables
      ## Check if the directory has the entire US naming convention or state naming convention
      tableNames <- names(inTables)
      outTables <- list()
      if (any(str_sub(tableNames, 3, 3) == '_')){ ## STATE NAMING CONVENTION
        if (anyDuplicated(str_sub(tableNames, 4)) != 0){
          for (i in 1:length(unique(str_sub(tableNames, 4)))){
            subList <- inTables[str_sub(tableNames, 4) == unique(str_sub(tableNames,4))[i]]
            name <- unique(str_sub(tableNames, 4))[i]
            # Give a ton of warnings about factors and characters, don't do that
            outTables[[name]] <- suppressWarnings(do.call(rbind, subList))
          }
        } else {
          outTables <- inTables
          names(outTables) <- unique(str_sub(tableNames, 4))
        }
      } else{ ## ENTIRE NAMING CONVENTION
        if (anyDuplicated(tableNames) != 0){
          for (i in 1:length(unique(tableNames))){
            subList <- inTables[tableNames == unique(tableNames)[i]]
            name <- unique(str_sub(tableNames, 1))[i]
            # Give a ton of warnings about factors and characters, don't do that
            outTables[[name]] <- suppressWarnings(do.call(rbind, subList))
          }
        } else {
          outTables <- inTables
          names(outTables) <- unique(str_sub(tableNames, 1))
        }
      }
      # NEW CLASS NAME FOR FIA DATABASE OBJECTS
      #outTables <- lapply(outTables, as.data.frame)
      class(outTables) <- 'FIA.Database'
    }



    #### DOWNLOADING THE WHOLE ZIP FILE
  } else {
    ## Download to a temporary file location, then extract to the permanent if wanted. Read extracted tables into R w/ readFIA.

    # Make sure state Abbs are in right format
    states <- str_to_upper(states)

    ## Download each state and extract to directory
    for (i in 1:length(states)){
      # Temporary directory to download to
      temp <- paste0(tempDir, '/', states[i],'_CSV.zip') #tempfile()
      ## Make the URL
      url <- paste0('https://apps.fs.usda.gov/fia/datamart/CSV/', states[i],'_CSV.zip')
      if (states == 'REF') {url <- paste0('https://apps.fs.usda.gov/fia/datamart/CSV/FIADB_REFERENCE.zip')}
      #newName <- paste0(str_sub(url, 1, -4), 'csv')
      ## Download as temporary file
      download.file(url, temp, timeout = 3600)
      ## Extract
      if (is.null(dir)){
        unzip(temp, exdir = tempDir)
      } else {
        unzip(temp, exdir = str_sub(dir, 1, -2))
      }
      unlink(temp)
    }

    if (load){

      ## Read in the files w/ readFIA
      if (is.null(dir)){
        outTables <- readFIA(tempDir, nCores = nCores, common = common, states = states)
        #unlink(tempDir, recursive = TRUE)
      } else {
        outTables <- readFIA(dir, nCores = nCores, common = common, states = states)
      }
    }



  }


  if (is.null(dir)){
    tmp <- list.files(tempDir, full.names = TRUE, pattern = '.csv')
    invisible(file.remove(tmp))
  }

  if (load) return(outTables)
}


## Read FIA data from a database connection
readFIA.db <- function(con,
                       schema,
                       states,
                       tables,
                       common,
                       inMemory) {

  ## Load all data now, or wait until later (process state-by-state)?
  if (is.null(inMemory) | !is.logical(inMemory)) stop('`inMemory` must be TRUE/FALSE.')
  if (inMemory) {

    ## Check that specified states are valid -------------------------------------
    if (is.null(states)) states <- pull_all_states(con, schema)
    states <- unique(stringr::str_to_upper(states))
    stateNames <- intData$stateNames
    badStates <- states[!c(states %in% stateNames$STATEAB)]
    if (length(badStates) > 0) {
      if (any(stringr::str_length(badStates) > 2)) {
        stop (paste0('The following states are not supported: ',
                     paste(badStates, collapse = ', '),
                     '. Note that state abbreviations are required in place of full names (e.g., "MI" instead of "Michigan").'))
      } else {
        stop (paste0('The following states are not supported: ',
                     paste(badStates, collapse = ', '), '.'))
      }
    }


    ## Check that specified tables are valid -------------------------------------
    commonTables <- c('COND', 'COND_DWM_CALC', 'INVASIVE_SUBPLOT_SPP', 'PLOT',
                      'POP_ESTN_UNIT','POP_EVAL', 'POP_EVAL_GRP', 'POP_EVAL_TYP',
                      'POP_PLOT_STRATUM_ASSGN', 'POP_STRATUM', 'SUBPLOT', 'TREE',
                      'TREE_GRM_COMPONENT', 'TREE_GRM_MIDPT', 'TREE_GRM_BEGIN',
                      'SUBP_COND_CHNG_MTRX', 'SEEDLING', 'SURVEY', 'SUBP_COND',
                      'P2VEG_SUBP_STRUCTURE')
    allTables <- sort(c(commonTables,
                        'COUNTY', 'SUBP_COND', 'BOUNDARY', 'TREE_WOODLAND_STEMS',
                        'TREE_REGIONAL_BIOMASS', 'TREE_GRM_ESTN', 'SITETREE',
                        'P2VEG_SUBPLOT_SPP', 'GRND_CVR', 'DWM_VISIT',
                        'DWM_COARSE_WOODY_DEBRIS', 'DWM_DUFF_LITTER_FUEL',
                        'DWM_FINE_WOODY_DEBRIS', 'DWM_MICROPLOT_FUEL',
                        'DWM_RESIDUAL_PILE', 'DWM_TRANSECT_SEGMENT',
                        'COND_DWM_CALC', 'PLOT_REGEN', 'SUBPLOT_REGEN',
                        'SEEDLING_REGEN', 'POP_EVAL_ATTRIBUTE', 'PLOTGEOM',
                        'PLOTSNAP'))

    # If `tables` isn't specified, default to common or all
    if (is.null(common) | !is.logical(common)) stop('`common` must be TRUE/FALSE.')
    if (is.null(tables) & common) tables <- commonTables
    if (is.null(tables) & !common) tables <- allTables

    # Checking all tables are cool
    tables <- unique(stringr::str_to_upper(tables))
    badTables <- tables[!c(tables %in% allTables)]
    if (length(badTables) > 0) {
      stop (paste0('The following tables are not supported: ',
                   paste(badTables, collapse = ', '), '.'))
    }


    ## Load tables -------------------------------------------------------------
    out <- lapply(X = tables,
                  FUN = read_fia_table_from_db,
                  con = con,
                  schema = schema,
                  states = states)
    names(out) <- tables

    ## POP_EVAL_TYP doesn't have a STATECD variable, but we can filter it by
    ## EVAL_CN if available. Otherwise you get the whole thing.
    if ('POP_EVAL' %in% tables & 'POP_EVAL_TYP' %in% tables) {
      out$POP_EVAL_TYP <- out$POP_EVAL_TYP %>%
        dplyr::filter(EVAL_CN %in% out$POP_EVAL$CN)
    }

    ## Add package class and done
    class(out) <- 'FIA.Database'



    # Set up "remote" FIA Database
  } else {

    # If states not specified, grab all available from the plot table
    if (is.null(states)) {
      message('`states` not provided. Defaulting to all states listed in PLOT table.')
      allStates <- dplyr::tbl(
        con,
        dplyr::sql(
          paste0("SELECT DISTINCT statecd FROM ", schema, '.plot')
        )
      ) %>%
        dplyr::collect() %>%
        dplyr::rename_with(stringr::str_to_upper) %>%
        dplyr::left_join(intData$stateNames, by = 'STATECD')

      states <- allStates$STATEAB
    }

    ## Saving the call to readFIA, for eval later
    out <- list(con = con,
                schema = schema,
                common = common,
                tables = tables,
                states = states)
    class(out) <- 'Remote.FIA.Database'

  }

  return(out)

}

## Load individual tables from a connection to a mirror of the FIADB
read_fia_table_from_db <- function(table, con, schema, states) {

  ## Reformat states so we can use "IN"
  statecds <- intData$stateNames$STATECD[intData$stateNames$STATEAB %in% states]
  states <- paste0('(', paste(statecds, collapse = ', '), ')')

  ## All tables except POP_EVAL_TYP have a STATECD variable, so we can filter
  ## tables by state and then load. For POP_EVAL_TYP we have to load the whole
  ## thing each time, and then apply an adhoc filter from EVAL_CN if POP_EVAL
  ## is also available
  query <- paste0("SELECT * FROM ", schema, '.', table)
  if (table != 'POP_EVAL_TYP') query <- paste0(query, ' WHERE statecd IN ', states)

  ## Load the table
  fiaTable <- dplyr::tbl(
    con,
    dplyr::sql(query)
  ) %>%
    dplyr::collect() %>%
    dplyr::rename_with(stringr::str_to_upper)

  return(fiaTable)
}

# If states not specified, grab all available from the plot table
pull_all_states <- function(con, schema) {
  allStates <- dplyr::tbl(
    con,
    dplyr::sql(
      paste0("SELECT DISTINCT statecd FROM ", schema, '.plot')
    )
  ) %>%
    dplyr::collect() %>%
    dplyr::rename_with(stringr::str_to_upper) %>%
    dplyr::left_join(intData$stateNames, by = 'STATECD')
  return(allStates$STATEAB)
}



## Write out the raw FIA files
#' @export
writeFIA <- function(db,
                     dir,
                     byState = FALSE,
                     nCores = 1,
                     ...){
  if (class(db) == 'Remote.FIA.Database'){
    stop('Cannot write remote database.')
  }

  # if (byState & !c('SURVEY' %in% names(db))){
  #   stop('Need survey table for state abbreviations.')
  # }

  #cat(sys.call()$dir)
  if (!is.null(dir)){
    # Add a slash to end of directory name if missing
    if (str_sub(dir,-1) != '/'){
      dir <- paste(dir, '/', sep = "")
    }
    # Check to see directory exists, if not, make it
    if(!dir.exists(dir)) {
      dir.create(dir)
      message(paste('Creating directory:', dir))
    }
    message(paste0('Saving to ', dir, '. NOTE: modifying FIA tables in Excel may corrupt csv files.'))
  }

  ## Method to chunk up the database into states before writing it out
  if (byState){

    db$PLOT <- db$PLOT %>%
      select(-c(any_of('STATEAB'))) %>%
      left_join(intData$stateNames, by = 'STATECD')

    ## Unique state abbreviations
    states <- unique(db$PLOT$STATEAB)

    ## Chunk up plot
    pltList <- split(db$PLOT, as.factor(db$PLOT$STATEAB))

    ## Loop over states, do the writing
    for (s in 1:length(states)){
      db_clip <- db
      ## Overwriting plot with shortened table
      db_clip$PLOT <- pltList[[s]]
      ## Subsetting the remaining database based
      ## on plot table
      db_clip <- clipFIA(db_clip, mostRecent = FALSE)

      ## Write it all out
      tableNames <- names(db_clip)[names(db_clip) != 'mostRecent']
      tableNames <- paste(unique(db_clip$PLOT$STATEAB), tableNames, sep = '_')
      for (i in 1:length(tableNames)){
        if (is.data.frame(db_clip[[i]])){
          fwrite(x = db_clip[[i]], file = paste0(dir, tableNames[i], '.csv'), showProgress = FALSE, nThread = nCores)
        }
      }
    }


    ## Merge states together on writing
  } else {
    tableNames <- names(db)

    for (i in 1:length(tableNames)){
      if (is.data.frame(db[[i]])){
        fwrite(x = db[[i]], file = paste0(dir, tableNames[i], '.csv'), showProgress = FALSE, nThread = nCores)
      }
    }
  }

}

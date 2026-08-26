#------------------------------------------------------------------------------#
# script: CleanNames.R, formerly adapted for the LRP scallop project
#
# objectives: Take the column names of a dataframe and standardize them somewhat.
#   It does this by:
#     - removing leading or trailing white space, 
#     - replacing spaces with underscores, 
#     - converting symbols like % and # with abbreviations and
#     - making all lower case (including converting camelCase to camel_case).
#     - I added a line to the end which changes underscores to periods.
#
# author : William Doane 2019 
#   (https://www.r-bloggers.com/2019/07/clean-consistent-column-names/)
# adapted by : Andrew Harbicht
# date : 2022
# update by : LLandry
#
# input: 
#
# output: 
# - CleanNames function
#   which requires the attributes: .data, unique = FALSE (default)
#
# information/references:
#   to help to understand "regular expression" (regex)
#     see https://regex101.com/r/9HLin6/1 
#
# todo:
#
#------------------------------------------------------------------------------#
require(data.table)

CleanNames <- function(.data, unique = FALSE, separator = "dot") {
  
  n <- if (is.data.frame(.data)) {colnames(.data)} else {.data}
  
  # replace symbols with text equivalent #######################################
  # plus surrounded by underscores
  n <- gsub("%+", "_pct_", n)
  n <- gsub("\\$+", "_dollars_", n)
  n <- gsub("\\++", "_plus_", n)
  n <- gsub("-+", "_minus_", n)
  n <- gsub("\\*+", "_star_", n)
  n <- gsub("#+", "_cnt_", n)
  n <- gsub("&+", "_and_", n)
  n <- gsub("@+", "_at_", n)
  
  # replace residual symbols and whitespace(s) by underscore  ##################
  # replace a single non-alphanumeric/underscore/whitespace(s) character by one underscore
  n <- gsub("[^a-zA-Z0-9_]+", "_", n)
  
  # put an underscore before capital letter followed by a lowercase letter #####
  n <- gsub("([A-Z][a-z])", "_\\1", n)
  
  # remove leading and/or trailing whitespace and convert all to lower case ####
  n <- tolower(trimws(n))
  
  # remove leading and trailing underscore(s) ##################################
  n <- gsub("(^_+|_+$)", "", n)
  
  # replace double underscores with single #####################################
  n <- gsub("_+", "_", n)
  
  # make the variable names unique if the attribute unique == TRUE #############
  # by adding an underscore followed by a distinct number
  if (unique) n <- make.unique(n, sep = "_")
 
  # replace underscores with periods ###########################################
  if (separator == "dot"){
    n <- gsub("_", ".", n)
  }else{
    if (separator == "underscore"){
      cat("underscores are used as separator in column variables")
    }
  }
  
  # put the new clean names into the column headers ############################
  if (is.data.frame(.data)) {
    colnames(.data) <- n
    .data
  } else {
    n
  }
} # end of function CleanNames

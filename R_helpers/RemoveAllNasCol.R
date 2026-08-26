#'
#' RemoveAllNasCol : remove string col or vector of string col or
#'                   all NA value col.
#'
#' @param input_data 
#'
#' @returns
#' @export
#'
#' @examples
#'

RemoveAllNasCol <- function(input_data){
  
  for (col in names(input_data)){
    vals <- input_data[[col]]
    
    if (all(is.na(vals) | grepl("^\\s*NA(\\s*,\\s*NA)*\\s*$", vals))){
      input_data[[col]] <- NULL
      
    } 
  }
  return(input_data)
}
#' A function to load in sample data
#'
#' This function allows you to load in sample-data.
#' @param
#' @keywords
#' @export
#' @examples
#' load_data()

load_data <- function(){
  giturl = "https://github.com/yuanjing-ma/RioNorm2_simulation/blob/master/sample-data.RData?raw=true"
  load(url(giturl))
  data = sample
  return (data)
}


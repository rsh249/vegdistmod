#' Download distribution data from NEOTOMA
#' 
#' This function requests data from the NEOTOMA database for a single taxon.
#' 
#' @param taxon A string of the form 'genus species' or 'genus'.
#' @export
#' @examples \dontrun{
#' abies <- rNeotoma('Abies');
#' }

neotoma.sp <- function(taxon) {
  #require('jsonlite');
  #require('urltools')
  taxon = paste(taxon, "*", sep ='');
  taxon = urltools::url_encode(taxon)
  
  html_str = paste("http://api.neotomadb.org/v1/data/sampledata?taxonname=", taxon, sep = '');
  jsonget = jsonlite::fromJSON(html_str);
  return(jsonget$data)
}  

#' #' Download data set information from NEOTOMA
#' #' 
#' #' This function requests data from the NEOTOMA database for a given dataset.
#' #' 
#' #' @param taxon A string giving the dataset ID.
#' #' @export
#' #' @examples \dontrun{
#' #' abies <- rNeotoma('Abies');
#' #' }
#' neotoma.setinfo <- function(data_set_id) {
#'   require('jsonlite');
#'   require('urltools')
#'   data_set_id = url_encode(data_set_id)
#'   #http://api.neotomadb.org/v1/data/sites?sitename=**&format=xml
#'   site_list = "http://api.neotomadb.org/v1/data/sites?sitename=**&format=json";
#'   jsonget=fromJSON(site_list)
#'   ret = subset(jsonget$data, jsonget$data$SiteID == 1766); return(ret);
#'   html_str = paste("http://api.neotomadb.org/v1/data/sites?sitename=", sitename, "&format=json", sep = '');
#'   jsonget = fromJSON(html_str);
#'   # return(html_str)
#'   return(jsonget)
#' } 
#' 
#' #' Download data set IDs from NEOTOMA
#' #' 
#' #' This function requests ID values from the NEOTOMA database for all NEOTOMA data sets.
#' #' 
#' #' @param taxon A string of the form 'genus species' or 'genus'.
#' #' @export
#' #' @examples \dontrun{
#' #' abies <- rNeotoma('Abies');
#' #' }
#' neotoma.setdata <- function(data_set_id) {
#'   require('jsonlite');
#'   require('urltools')
#'   data_set_id = url_encode(data_set_id)
#'   #http://api.neotomadb.org/v1/data/sites?sitename=**&format=xml
#'   site_list = "http://api.neotomadb.org/v1/data/sites?sitename=**&format=json";
#'   jsonget=fromJSON(site_list)
#'   sitename = url_encode(subset(jsonget$data, jsonget$data$SiteID == 1766)$SiteName)
#'   html_str = paste("http://api.neotomadb.org/v1/data/sampledata?sitename=", sitename, "&format=json", sep = '');
#'   jsonget = fromJSON(html_str);
#'   # return(html_str)
#'   return(jsonget)
#' } 
#' 
#'  

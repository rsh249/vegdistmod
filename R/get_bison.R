#' Download distribution data from BISON (https://bison.usgs.gov)
#' 
#' This function requests data from the BISON database for a single species. 
#' Note that BISON requires exact name matching to binomial. Searching on a
#' genus, for example, will match only records that are enterred with that name
#' only, NOT all records in that genus identified by a binomial. Also note that
#' BISON has no name correction so mispelling and errors are likely.
#' 
#' @param taxon A string of the form 'genus species'.
#' @param maxrec Limit on number of records to download.
#' @export
#' @examples \dontrun{
#' abies <- get_bison('Abies fraseri', 10000);
#' }


get_bison <- function(taxon, maxrec = 10000) {
  #require('jsonlite');
  #require('urltools')
  taxon = urltools::url_encode(taxon)
  html_str = paste("https://bison.usgs.gov/api/search.json?species=", taxon, "&type=scientific_name&start=0&count=", maxrec, sep = '');
  #return(html_str); 
  jsonget = jsonlite::fromJSON(html_str); 
  #return(jsonget)
  dat = jsonget$data;
  dat = subset(dat, dat$geo != "No");
  if(class(dat)=='list'){
    stop('no data from bison\n')
    return(NA)
  }
  return(dat)
}  
  
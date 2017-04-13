#' Download distribution data from iNaturalist
#' 
#' This function requests data from the iNaturalist database for a single species that
#' wraps the data access function(s) from the rinat library.
#' 
#' @param taxon A string of the form 'genus species'.
#' @param maxrec Limit on number of records to download.
#' @export
#' @examples \dontrun{
#' abies <- inat('Abies fraseri', 10000);
#' }

inat <- function(taxon, maxrec = 10000){
  #require(rinat);
  di = rinat::get_inat_obs(taxon_name=taxon, maxresults=maxrec)
  return(di);
}

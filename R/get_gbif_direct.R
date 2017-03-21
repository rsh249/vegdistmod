#' Download distribution data directly from GBIF API
#' 
#' This function requests data from the GBIF database for a single taxon using the GBIF callback API.
#' 
#' @param taxon A string of the form 'genus species' or 'genus'.
#' @param maxrec Maximum number of records to download.
#' @export
#' @examples \dontrun{
#' abies <- gbif_get('Abies');
#' }

gbif_get <- function(taxon, maxrec=Inf) { 
  
 # require('jsonlite');
  #require('urltools')
  n=0;
  round = 0;
  hold = list();
  offset=0;
  while(n < 1){

    
    taxon = urltools::url_encode(taxon)
    html_str = paste("https://api.gbif.org/v1/occurrence/search?scientificName=", taxon, "&limit=300&offset=", offset, sep = ''); 
    jsonget = jsonlite::fromJSON(html_str);
    round = round+1;
    if(is.null(nrow(jsonget$results))){ print("ERR: 11");return(jsonget)} else {
      hold[[round]] = jsonget$results;
    }  
    
    if(jsonget$endOfRecords == TRUE){
      n=1;
    } else {
      offset = offset+300;
    }
    if(offset > maxrec){
      break;
    }
    
    


  }
  df = hold[[1]][,c('key', 'genus', 'specificEpithet', 'decimalLongitude', 'decimalLatitude')]
  for (n in 2:length(hold)){
    print(n);
    nex = hold[[n]][,c('key', 'genus', 'specificEpithet', 'decimalLongitude', 'decimalLatitude')]
    df = rbind(df, nex)
    
  }
  return(df)
}



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

gbif_get <- function(taxon, maxrec = 200000) {
  # require('jsonlite');
  #require('urltools')
  n = 0
  
  round = 0
  
  hold = list()
  
  offset = 0
  
  taxon = urltools::url_encode(taxon)
  
  while (n < 1) {
    html_str = paste(
      "https://api.gbif.org/v1/occurrence/search?scientificName=",
      taxon,
      "&limit=300&offset=",
      offset,
      sep = ''
    )
    
    jsonget = jsonlite::fromJSON(html_str)
    
    round = round + 1
    
    if (is.null(nrow(jsonget$results))) {
      print("ERR: 11")
      return(jsonget)
    } else {
      hold[[round]] = jsonget$results
      
    }
    
    if (jsonget$endOfRecords == TRUE) {
      n = 1
      
    } else {
      offset = offset + 300
      
    }
    if (offset > maxrec) {
      break
      
    }
    
    
    
    
  }
  df = hold[[1]][, c('key',
                     'genus',
                     'specificEpithet',
                     'decimalLatitude',
                     'decimalLongitude'
                     )]
  if (length(hold) > 1) {
    for (n in 2:length(hold)) {
      # print(n);
      nex = hold[[n]][, c('key',
                          'genus',
                          'specificEpithet',
                          'decimalLatitude',
                          'decimalLongitude'
                          )]
      df = rbind(df, nex)
      
    }
  }
  df[,2] = paste(df[,2], df[,3]);
  df = df[,-3];
  
  colnames(df) = c('ind_id', 'tax', 'lat', 'lon')
  return(df)
}

#' Download distribution data from BIEN, GBIF, Inaturalist,
#'
#' This function requests data from the GBIF database for a single taxon using the GBIF callback API.
#'
#' @param taxon A string of the form 'genus species' or 'genus'.
#' @param maxrec Maximum number of records to download.
#' @param local TRUE or FALSE. To use the GBIF API use FALSE
#' @export
#' @examples \dontrun{
#' abies <- gbif_dist_all('Abies', maxrec = 1000, local = FALSE);
#' }
get_dist_all <- function(taxon, maxrec = 19999, local = FALSE) {
  ###GET DATA
  #GET GBIF DATA direct
  library(vegdistmod)
  gbif = cbind(1,1,1,1);
  tryCatch({
    if(local == FALSE){
      gbif <- gbif_get(taxon, maxrec = 199999)
    } else {}
  },
  error = function(cond) {
    message(cond)
    return(NA)
  })
  ##Use .gbif_sql()
  ##set flag: if(nrow(ab.cgbd)>200000){print 'taxon' to file}
  
  
  
  #GET BIEN DATA
  bien = cbind(1,1,1,1);
  tryCatch({
  bien <-
    BIEN::BIEN_occurrence_species(
      species = taxon,
      native.status = TRUE,
      only.new.world = TRUE
    )
  },
  error = function(cond) {
    message(cond)
    return(NA)
  })
  
  #get bison data
  bison = cbind(1,1,1,1)
  tryCatch({
    bison <- vegdistmod::get_bison(taxon, maxrec = maxrec)
    
  },
  error = function(cond) {
    message(cond)
    return(NA)
  })
  
  #get inaturalist data
  source("~/Desktop/cracle_testing/development_files/rInat.R")
  inatr = cbind(1,1,1,1);
  tryCatch({
    inatr = vegdistmod::inat(taxon)
  }, 
  error = function(cond) {
    message(cond)
    return(NA)
  })
  
  
  cnames <- c('ind_id', 'tax', 'lat', 'lon')
  
  if (nrow(inatr) > 5) {
    inatr = inatr[, c('id', 'scientific_name', 'latitude', 'longitude')]
    colnames(inatr) = cnames
  } else {
    inatr = NA
  }
  if (nrow(bison) > 5) {
    bison = bison[, c('occurrenceID',
                      'name',
                      'decimalLatitude',
                      'decimalLongitude')]
    colnames(bison) = cnames
  } else {
    bison = NA
  }
  if (nrow(gbif) > 5) {
    #gbif = gbif[, c('key',
    #                   'genus',
    #                   'specificEpithet',
    #                   'decimalLatitude',
    #                   'decimalLongitude')]
    #gbif[, 2] = paste(gbif[, 2], gbif[, 3], sep = ' ')
    #gbif = gbif[,-3]
    #colnames(gbif) = cnames
  } else {
    gbif = NA
  }
  
  if (nrow(bien) > 5) {
    bien = bien[, c('datasource_id',
                    'scrubbed_species_binomial',
                    'latitude',
                    'longitude')]
    colnames(bien) = cnames
  } else {
    bien = NA
  }
  data <- rbind(inatr, bison, gbif, bien)
  
  data$lat <- as.numeric(as.character(data$lat))
  data$lon <- as.numeric(as.character(data$lon))
  data = subset(data, data$tax == taxon)
  data = stats::na.omit(data)
  return(data)
}

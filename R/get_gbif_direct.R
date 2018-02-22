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
  tori = taxon;
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
      return(NULL)
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
  cols = c('key',
  'genus',
  'specificEpithet',
  'decimalLatitude',
  'decimalLongitude');
  
  if(sum(cols %in% names(hold[[1]]))==length(cols)) {
    df = hold[[1]][, c(cols )]
    if (length(hold) > 1) {
      for (n in 2:length(hold)) {
        # print(n);
        if(sum(cols %in% colnames(hold[[n]]))==4){
      
        nex = hold[[n]][, c(cols)]
        df = rbind(df, nex)
      } else {next;}
      
    }

  } 
  
  df[,2] = paste(df[,2], df[,3]);
  df = df[,-3];
  
  colnames(df) = c('ind_id', 'tax', 'lat', 'lon')
  #df$tax = rep(tori, nrow(df));
  return(df)
  } else { 
    return(NULL); 
    }
}

#' Download distribution data from BIEN, GBIF, Inaturalist,
#'
#' This function requests data from the GBIF database for a single taxon using the GBIF callback API.
#'
#' @param taxon A string of the form 'genus species' or 'genus'.
#' @param maxrec Maximum number of records to download.
#' @param local TRUE or FALSE. To use the GBIF API use FALSE
#' @param db SQL database. ONLY FOR LOCAL GBIF INSALLATION
#' @param h SQL host
#' @param u SQL user
#' @param pass SQL password
#' @export
#' @examples \dontrun{
#' abies <- gbif_dist_all('Abies', maxrec = 1000, local = FALSE);
#' }
get_dist_all <- function(taxon, maxrec = 19999, local = FALSE, db = 0, h = 0, u = 0, pass = 0) {
  ###GET DATA
  #GET GBIF DATA direct
  library(vegdistmod)
  gbif = cbind(1,1,1,1);
  tryCatch({
    if(local == FALSE){
      gbif <- gbif_get(taxon, maxrec = 199999)
    } else {
      gbif <- .gbif_local(taxon, limit = maxrec, db = db, h = h, u = u, pass = pass)
    }
  },
  error = function(cond) {
    message(paste("GBIF", cond))
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
    message(paste("BIEN", cond))
    return(NA)
  })
  
  #get bison data
  bison = cbind(1,1,1,1)
  tryCatch({
    bison <- vegdistmod::get_bison(taxon, maxrec = maxrec)
    
  },
  error = function(cond) {
    message(paste("BISON", cond))
    return(NA)
  })
  
  #get inaturalist data
  #source("~/Desktop/cracle_testing/development_files/rInat.R")
  inatr = cbind(1,1,1,1);
  tryCatch({
    inatr = vegdistmod::inat(taxon, maxrec = maxrec)
  }, 
  error = function(cond) {
    message(paste("inat", cond))
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

#' Download distribution data, filter, and merge with climate or environmental
#'
#' getextr is a function that gets GBIF data and extracts climate or environmental 
#' data for each occurrence. This is a whole workflow for distribution 
#' data acquisition and value addition that draws on several other functions in vegdistmod
#' including gbif_get and extraction. Parallel option is useful for speeding up data collection for
#' many species when computational resources are available.
#' 
#' @param x A taxon name or list of taxon names. It is sometimes good to 
#' test these on the vegdistmod::get_gbif() function first.
#' @param maxrec Maximum number of records to download.
#' @param clim A raster object of climate or other environmental data to extract from.
#' @param schema To be passed to vegdistmod::extraction
#' @param rm.outlier To be passed to vegdistmod::extraction
#' @param factor To be passed to vegdistmod::extraction
#' @param alpha To be passed to vegdistmod::extraction
#' @param nmin To be passed to vegdistmod::extraction
#' @param parallel TRUE or FALSE. Should this be executed in parallel.
#' @param nclus If parallel == TRUE then how many cores should be used? Default is 4.
#' 
#' @export
#' @examples \dontrun{
#' abies <- getextr('Abies fraseri', 
#' clim = clim, maxrec=500, 
#' schema= 'flat', rm.outlier = TRUE, 
#' alpha=0.01, factor = 2, nmin = 5, parallel=FALSE, nclus = 4));
#' }
#' 
getextr = function(x, clim = clim, maxrec=500, schema= 'flat', 
                   rm.outlier = TRUE, alpha=0.01, 
                   factor = 2, nmin = 5, parallel=FALSE, nclus = 4){
  
  clim = clim;
  maxrec = maxrec;
  schema = schema;
  rm.outlier = rm.outlier;
  alpha = alpha;
  factor = factor;
  nmin = nmin;
  parallel = parallel;
  nclus = nclus;
  
  
  subfun = function(x){
    ex = list();
    for(i in 1:length(x)){
      print(x[i]);
      ex[[i]] = NULL;
      dat2 = vegdistmod::gbif_get(x[i], maxrec = maxrec)
      if(is.null(dat2)){ ex[[i]]=NULL; next; }
      dat2 = stats::na.omit(dat2);
      if(any(is.na(dat2))){ ex[[i]]=NULL; next;}
      if(nrow(dat2)<nmin){ ex[[i]]=NULL; next; }
      ex.hold = vegdistmod::extraction(dat2, clim, 
                                       schema = schema, 
                                       rm.outlier = rm.outlier, 
                                       alpha = alpha, 
                                       factor = factor, 
                                       nmin = nmin);
      if(length(ex.hold) == 0){ ex[[i]] = NA;next;} else {
        ex.hold$tax = rep(x[i], nrow(ex.hold))
        ex[[i]] = ex.hold;
      }
    }
    
    ex = stats::na.omit(ex);
    #	if(any(is.null(ex))){ return(NULL); }
    if(length(ex) == 0) { return(NULL); }
    
    ex2 = rbind(ex[[1]]);
    if(length(ex)>1){
      for(k in 2:length(ex)){
        ex2 = rbind(ex2, ex[[k]]);
      }
    } else { return(ex); }
    
    return(ex2);
  }
  
  
  if(parallel==FALSE){
    return(subfun(x));
  } else {
    clim = clim;
    maxrec = maxrec;
    schema = schema;
    rm.outlier = rm.outlier;
    alpha = alpha;
    factor = factor;
    nmin = nmin;
    parallel = parallel;
    nclus = nclus;
    
    cl = parallel::makeCluster(nclus, type = "SOCK", outfile = '')
    
    parallel::clusterExport(cl, varlist = c('clim',  'maxrec', 'nmin', 'schema', 'rm.outlier', 'alpha', 'factor' ), envir = environment())
    splits = parallel::clusterSplit(cl, x);
    extr = parallel::parLapply(cl, splits, subfun);
    parallel::stopCluster(cl);
    return(extr);
    
    ##code below here not executed and problematic::
    extall = rbind(extr[[1]][[1]]); ##Need to check that this object is OK as below.
    
    for(k in 2:length(extr)){
      
      if(is.null(extr[[k]])){} else {
        if(ncol(extr[[k]][[1]])==1){} else {
          extall=rbind(extall, extr[[k]][[1]]);
        }
      }
    }
    return(extall);
    
  }
}

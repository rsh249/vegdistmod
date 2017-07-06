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
  
  ###Run test query here to get total number of records so can adjust maxrecs
  html_test = paste("https://bison.usgs.gov/api/search.json?species=", 
                   taxon, "&type=scientific_name&start=0&count=1", sep = '');
  #print(html_test); 
  
  jsontest = jsonlite::fromJSON(html_test); 
  #return(jsontest)
  ntotal = sum(unlist(jsontest$occurrences$legend))
  if(maxrec>ntotal){maxrec=ntotal}
  dat = matrix(ncol=8, nrow=maxrec+100);
  dat = as.data.frame(dat);
  tdat = jsontest$data;
  colnames(dat)=colnames(tdat);
  if(maxrec<500){recs = maxrec}else{recs=500}
  if(maxrec>500){
    i = maxrec/500; print(i);
    for(z in 0:i){
      n=z*500;
      if(n>maxrec){break}
      if(n+500>maxrec){recs = maxrec-n}
      html_str = paste("https://bison.usgs.gov/api/search.json?species=", 
                       taxon, "&type=scientific_name&start=", n, "&count=", 
                       recs, sep = '');
     # print(html_str); 
      
      jsonget = jsonlite::fromJSON(html_str); 
      #return(jsonget)
      d = jsonget$data;
      #d2 = subset(d, dat$geo != "No");
      #print(d)
     # print(n)
     # print(n+nrow(d))
     # print(nrow(dat))
      dat[(n+1):(n+nrow(d)),]=d;
      if(class(dat)=='list' & n == 0){
        stop('no data from bison\n')
        return(NA)
      }
    }
    
  } else {
    n=0
    html_str = paste("https://bison.usgs.gov/api/search.json?species=", 
                   taxon, "&type=scientific_name&start=", n, "&count=", 
                   recs, sep = '');
   # print(html_str); 
  
    jsonget = jsonlite::fromJSON(html_str); 
    #return(jsonget)
    d = jsonget$data;
   #d2 = subset(d, dat$geo != "No");
   # print(d)
    dat[n+1:(n+nrow(d)),]=d;
    if(class(dat)=='list'){
      stop('no data from bison\n')
      return(NA)
    }
  
  }  
  
  ##Opportunity to provide filtering options her
  ##geo==TRUE
  
  ##fossil==FALSE
  
  
  return(dat)
}  
  
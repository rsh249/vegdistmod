


.gbif_local <- function(taxon, limit=1000000000, db, u, h, pass){
  #this will only work on Cuvier
  #if(!is.null(grep("[[:space:]]", taxon))){
    split <- strsplit( taxon, ' ');
    genus = split[[1]][1];
    species = split[[1]][2];
    if(is.na(species)){
      species = "%%";
    }
  #} 
  query = paste("SELECT div_id, genus, species, lat, lon from div_base where genus = \'", genus, "\' and species = \'", species, "\' and is_bad !=1 and cultivated != 1 and no_precision != 1 and hasGeospatialIssues != 1 ", " LIMIT ", sprintf("%i", limit), sep='');
  #return(query)
  con = DBI::dbConnect(RMySQL::MySQL(), dbname=db, username=u, host = h, password = pass);
  get = DBI::dbGetQuery(con, query);
  DBI::dbDisconnect(con);
  return(get)
}


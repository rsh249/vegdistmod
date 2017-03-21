


.gbif_local <- function(taxon, limit=10000000000000){
  #this will only work on Cuvier
  if(!is.null(grep(" ", taxon))){
    split <- strsplit( taxon, ' ');
    genus = split[[1]][1];
    species = split[[1]][2];
  } else {
    genus = taxon;
    species = '%%';
  }
  
  con = DBI::dbConnect(RMySQL::MySQL(), dbname='rh_div_amnh', username='rharbert', host = 'localhost');
  query = paste("SELECT div_id, genus, species, lat, lon from div_base where genus = \'", genus, "\' and species = \'", species, "\' LIMIT ", limit, sep='');
  get = DBI::dbGetQuery(con, query);
  return(get)
}

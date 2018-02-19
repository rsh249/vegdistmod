#' Download WorldClim data for modern or paleo models.
#' 
#' This function requests climate raster objects from www.worldclim.org, 
#' downloads the files, reads them into R and disposes of the original files.
#' 
#' @param period A string. Either 'cur', 'midholo', or 'lgm'.
#' @param model For paleo models. Which to use (i.e., 'ccsm4'). See http://worldclim.org/paleo-climate1 for options.
#' @param varset Either 'bio', 'tmean', 'tmin', 'tmax', or 'prec'.
#' @param version Either '1_4', or '2.0'
#' @param res What spatial resolution? '10', '2.5' arcmin, or '30' seconds. (options are '10', '2.5', '30').
#' @export
#' @examples \dontrun{
#' #get 10 minute mean monthly temperature grids.
#'
#' abies <- get_worldclim(period='cur', varset = 'tmean', res=10); 
#' }
get_worldclim <- function(period = 'cur', model = '', version = '1_4', varset = 'bio', res = 2.5) {
 
  # period = 'cur';
  bs = '';
  spacer='climate/';
  version=version;
  modelstub = '';
  if(period == 'cur'){
    overhead = 'cur';
    bs = '_bil';
    spacer = paste('climate/worldclim/', version, '/grid/', sep = '');
  }
  if(period == 'midholo'){
    if(varset=='bio'){varset='bi'}
    if(varset=='tmin'){varset='tn'}
    if(varset=='tmax'){varset='tx'}
    if(varset=='tmean'){print("ERR: Variable not available in this model");return(0);}
    if(varset=='prec'){varset='pr'}
    overhead = 'cmip5/mid';
    if(model == 'ccsm4'){
      modelstub = 'ccmid';
    }
    if(model == 'BCC-CSM1-1'){
      modelstub = 'bcmid';
    }
    if(model == 'CNRM-CM5'){
      modelstub = 'cnmid';
    }
    if(model == 'HadGEM2-CC'){
      modelstub = 'hgmid';
    }
    if(model == 'HadGEM2-ES'){
      modelstub = 'hemid';
    }
    if(model == 'IPSL-CM5A-LR'){
      modelstub = 'ipmid';
    }
    if(model == 'MIROC-ESM'){
      modelstub = 'mrmid';
    }
    if(model == 'MPI-ESM-P'){
      modelstub = 'memid';
    }
    if(model == 'MRI-CGCM3'){
      modelstub = 'mgmid';
    }
  } 
  
  if(period == 'lgm'){
    overhead = 'cmip5/lgm';
    if(varset=='bio'){varset='bi'}
    if(varset=='tmin'){varset='tn'}
    if(varset=='tmax'){varset='tx'}
    if(varset=='tmean'){print("ERR: Variable not available in this model");return(0);}
    if(varset=='prec'){varset='pr'}
    if(model == 'ccsm4'){
      modelstub = 'cclgm';
    }
    
    if(model == 'MIROC-ESM'){
      modelstub = 'mrlgm';
    }
    if(model == 'MPI-ESM-P'){
      modelstub = 'melgm';
    }
  } 
  #varset = 'bio'; #bio, tmean, prec, tmin, tmax
  # res = '10'; #10, 5, 2-5, 30s*caveat that the 30s bioclim variables come in to sets 1-9 and 10-19
  if(res == 30){
    res = paste(res, 's', sep ='');
    if(varset == 'bio'){
      varset = c('bio1-9', 'bio10-19');
    }
  } else {
    res=chartr('.', '-', res)
    res = paste(res, 'm', sep ='');
  }
  if(version != 2 | version != "2.0"){
    varset = paste(modelstub, varset, sep = '');
    res = paste(res, bs, sep = '');
    var = paste(varset, res, sep = "_");
    http_str = 
      paste("http://biogeo.ucdavis.edu/data/", spacer, overhead, "/", var, ".zip", sep = '');
    
  }
  
  if(version == 2 | version == "2.0"){
    spacer = 'worldclim/v2.0/tif/base/wc2.0_';
    overhead = res;
    overhead = gsub('_bil', '', overhead)
    
    var = paste("_", varset, sep = '');
    http_str = 
      paste("http://biogeo.ucdavis.edu/data/", spacer, overhead, var, ".zip", sep = '');
    
  }
  
  
  #Version 2.0
  #http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_10m_vapr.zip
  
  #return(http_str)
  
  ##NOTE: TO SAVE TIME ON LARGE DOWNLOADS CHECK CURRENT DIRECTORY FOR FILES
  
  ###
  
  
  for(zz in 1:length(http_str)){
    temp = tempfile()
    td = tempdir();
    utils::download.file(http_str[[zz]], temp)
    list <- utils::unzip(temp, exdir=td);
    #return(list)
    if(length(grep(list, pattern = '*.bil'))>0){
      r= raster::stack(list[grep(list, pattern='*.bil')])
    } else if (length(grep(list, pattern = '*.tif'))>0){
      r= raster::stack(list[grep(list, pattern='*.tif')])
    } else {
      r = raster::stack(list);
    }
    raster::rasterOptions(maxmemory=1500000000) ##Set max ram for raster to 15GB
    r=raster::brick(r); #should be able to brick to memory now
  
    #unlink(temp)
    system(paste('rm ', temp));
    for(n in 1:length(list)){
      system(paste('rm ', list[[n]]))
    }
    if(zz > 1){
      rr=raster::brick(rr,r);
    } else {
      rr = r;
    }
  }
    return(rr);
}


#' Download ENVIREM topographical data set(s)
#' 
#' This function requests raster layers from the envirem dataset (Bemmels, et al. 2017, Ecography)
#' downloads the zip archive, reads files into R and disposes of the originals because you should 
#' save the layers in an R ready raster object using raster::writeRaster().
#' @param period A string for time period to get. Either 'cur', 'midholo', or 'lgm'.
#' @param res Spatial Resolution in arcminutes. eg, 2.5
#' @export
#' @examples \dontrun{
#' #get 2.5 arcmin grid. 
#' envir <- get_envirem_elev(period='cur'); 
#' }
get_envirem_elev <- function(period = 'cur', res=2.5) {
  res=res; if(res == ''){print("ERR: No resolution chosen"); return(0);}
  period = period;
  ##res is assumed to be 2.5 arcmin. need to add this as an option and pick according url.
  http_str = '';
  if(res == '2.5'){
    if(period == 'cur'){
      http_str = 'https://deepblue.lib.umich.edu/data/downloads/05741r787?locale=en' #new world elev params
    }
    if(period == 'midholo'){
      http_str = 'https://deepblue.lib.umich.edu/data/downloads/s4655g713?locale=en'
    }
    if(period == 'lgm'){
      http_str = 'https://deepblue.lib.umich.edu/data/downloads/8910jt72z?locale=en'
    }
  } 
  if(res == '30'){
    if(period == 'cur'){
      http_str = 'https://deepblue.lib.umich.edu/data/downloads/44558d44j?locale=en' #new world elev params
    }
    if(period == 'midholo'){
      http_str = '' ##Not available
    }
    if(period == 'lgm'){
      http_str = '' ##Not available
    }
    
  }
  
  
    if(http_str == ''){print("ERR: Invalid http string. Check model and period selections."); return(0);}
    temp = tempfile()
    td = tempdir();
    utils::download.file(http_str[[1]], temp)
    list <- utils::unzip(temp, exdir=td);
    #return(list)
    if(length(grep(list, pattern = '*.bil'))>0){
      r = raster::stack(list[grep(list, pattern='*.bil$')])
    } else if (length(grep(list, pattern = '*.tif'))>0){
      r = raster::stack(list[grep(list, pattern='*.tif')])
    } else {
      r  = raster::stack(list);
    }
    raster::rasterOptions(maxmemory=1500000000) ##Set max ram for raster to 15GB
    #consider commenting out the next line, but it does make it temp.file safe.
    r=raster::brick(r); #should be able to brick to memory now
    
    #unlink(temp)
    system(paste('rm ', temp));
    for(n in 1:length(list)){
      system(paste('rm ', list[[n]]))
    }
    ##update layer names
    n.arr = strsplit(names(r),'arcmin_')
    for (i in 1:length(n.arr)){
      n.arr[[i]] = n.arr[[i]][2]
    }
    n.arr=unlist(n.arr);
    names(r) = n.arr;
  return(r);
}

#' Download ENVIREM climate data set(s)
#' 
#' This function requests raster layers from the envirem dataset (Bemmels, et al. 2017, Ecography)
#' downloads the zip archive, reads files into R and disposes of the originals because you should 
#' save the layers in an R ready raster object using raster::writeRaster().
#' @param period A string for time period to get. Either 'cur', 'midholo', or 'lgm'.
#' @param model A string for which model to use. Either 'ccsm4', 'miroc-esm', or 'mpi-esm'.
#' @export
#' @examples \dontrun{
#' #get 2.5 arcmin grid for North America (only option currently.
#' envir <- get_envirem_clim(period='cur', model='');
#' }
get_envirem_clim <- function(period='cur', model='') {
  model = model;
  period = period;
  ##res is assumed to be 2.5 arcmin. need to add this as an option and pick according url.
  http_str = '';
  if(period == 'cur'){
    http_str = 'https://deepblue.lib.umich.edu/data/downloads/12579s441?locale=en';
  }
  if(period == 'midholo'){
    if(model == 'ccsm4') {
      http_str = 'https://deepblue.lib.umich.edu/data/downloads/70795779p?locale=en'
    }
    if(model == 'miroc-esm'){
      http_str = 'https://deepblue.lib.umich.edu/data/downloads/v118rd682?locale=en'
    }
    if(model == 'mpi-esm'){
      http_str='https://deepblue.lib.umich.edu/data/downloads/z029p4825?locale=en'
    }
  }
  if(period == 'lgm'){
    if(model == 'ccsm4') {
      http_str = 'https://deepblue.lib.umich.edu/data/downloads/9593tv32k?locale=en'
    }
    if(model == 'miroc-esm'){
      http_str = 'https://deepblue.lib.umich.edu/data/downloads/tm70mv28r?locale=en'
    }
    if(model == 'mpi-esm'){
      http_str='https://deepblue.lib.umich.edu/data/downloads/6969z089x?locale=en'
    }
  }
  
  temp = tempfile()
  td = tempdir();
  utils::download.file(http_str[[1]], temp)
  list <- utils::unzip(temp, exdir=td);
  #return(list)
  if(length(grep(list, pattern = '*.bil'))>0){
    r = raster::stack(list[grep(list, pattern='*.bil$')])
  } else if (length(grep(list, pattern = '*.tif'))>0){
    r = raster::stack(list[grep(list, pattern='*.tif')])
  } else {
    r  = raster::stack(list);
  }
  raster::rasterOptions(maxmemory=1500000000) ##Set max ram for raster to 15GB
  r=raster::brick(r); #should be able to brick to memory now
  
  #unlink(temp)
  system(paste('rm ', temp));
  for(n in 1:length(list)){
    system(paste('rm ', list[[n]]))
  }
  ##update layer names
  n.arr = strsplit(names(r),'arcmin_')
  for (i in 1:length(n.arr)){
    n.arr[[i]] = n.arr[[i]][2]
  }
  n.arr=unlist(n.arr);
  names(r) = n.arr;
  return(r);
}


#' @import classInt
#' @import grDevices
#' @import doParallel
NULL

# A function for retrieving data from cloud.diversityoflife.org
# 
# This function requests data from a private service that
# mirrors GBIF distribution data at cloud.diversityoflife.org. Account required.
# 
# @param taxon A string of the form 'genus' or 'genus_species'.
# @export
# @examples
# abies <- get_gbif_cloud('abies');
# abies.fraseri <- get_gbif_cloud('abies_fraseri') 


.get_gbif_cloud <- function(taxon) {
  #curl_string = paste("curl --user ", user, ":", pw, " http://cloud.diversityoflife.org/cgi-div/tmp_mat_get.pl?taxon=", sep = '');
  curl_string = paste("curl http://cloud.diversityoflife.org/docs/tmp_mat_get.pl?taxon=", sep = '');
  curl_string = paste(curl_string, taxon, sep = '')
  s = system(curl_string, intern = TRUE)
  if(length(s) < 2){return();}
  splithdr <- strsplit(s[1], split = '\t')
  df <- matrix(ncol = length(splithdr[[1]]), nrow = length(s) - 1)
  colnames(df) = rbind(unlist(splithdr[[1]]))
  for (i in 2:length(s)) {
    split <- strsplit(s[i], split = '\t')
    n = length(split[[1]])
    new = rbind(unlist(split[[1]]))
    df[i - 1,] = new
  }
  df <- data.frame(df)
  df$lon <- as.numeric(as.character(df$lon))
  df$lat <- as.numeric(as.character(df$lat))
  return(df)
}


.latlonfromcell <- function(cells, clim){
  #assumes square (in degrees) cells
  #returns lat/lon coordinates from the center of each cell
  d <- dim(clim)
  e <- raster::extent(clim)
  uplon <- e[1]
  uplat <- e[4]
  nrow = d[1]
  ncol = d[2]
  xres = (e[2] - e[1])/ncol;
  yres = (e[4] - e[3])/nrow;
  

  r <- raster::res(clim)
  m <- matrix(nrow = length(cells), ncol = 4);
  for(i in 1:length(cells)){
   # print(cells[[i]])
    row = ceiling(cells[[i]]/ncol); #round up so it is 1 : nrow. not 0 : nrow-1
    col = cells[[i]] - ((row-1)*ncol);
    
    lat = uplat - (row * yres) + (0.5 * yres)
    lon = uplon + (col * xres) - (0.5 * xres)
   # print(row); print(col); print(lat); print(lon);
    m[i,] = c("1111", "bg", lat, lon)
    
  }
  df = data.frame(m)
  df[,1] = as.numeric(as.character(df[,1]))
  df[,3] = as.numeric(as.character(df[,3]))
  df[,4] = as.numeric(as.character(df[,4]))
  
  return(df)
  
}
##Hidden function that is called under several others to get background data points.
# call with .get_bg()

.get_bg <- function(clim){
  n = 1000

  if(class(clim) == 'list'){
    m = list();
    for(a in 1:length(clim)){
      r = clim[[a]];
      d = dim(r);
      l = d[1] * d[2]
      v = seq(1:l);
      m[[a]] = .latlonfromcell(v, r);
      colnames(m[[a]]) = c('ind_id', 'tax', 'lat','lon');
      #bg <- randomPoints(r[[a]], n)
   #   bg <- cbind(rep('0000', length(bg[,1])), rep('bg', length(bg[,1])), bg[,2], bg[,1])
  #    colnames(bg) <- c('ind_id', 'tax', 'lat', 'lon');
   #   bg = as.data.frame(bg);
    #  bg$lon <- as.numeric(as.character(bg$lon));
     # bg$lat <- as.numeric(as.character(bg$lat));
    #  bg.l[[a]] = bg;
    }
    
    ret = m;
    return(ret);
  } else {
    
    r = clim;
    d = dim(r);
    l = d[1] * d[2]
    v = seq(1:l);
    m = .latlonfromcell(v, r);
    colnames(m) = c('ind_id', 'tax', 'lat','lon');
    
    return(m)
    
    
 
  }
}

#.get_bg = compiler::cmpfun(.get_bg);


#' The multivariate likelihood (log-likelihood) for a given set of PDF climate functions and localities.
#'
#'  Returns a log-likelihood calculated from the climate data at one 
#'  or more localities referenced against a set of climate PDFs.
#' @param x A data.frame of climate values (i.e., extracted from the climate raster object)
#' @param clim A raster object of climate data (matching x)
#' @param dens A density object from vegdistmod::densform(), vegdistmod::and_fun(), or vegdistmod::or_fun();
#' @param type Designate either ".gauss" or ".kde".
#' @param w To weight or not to weight by coefficient of variation.
#'
#' @export
#' @examples \dontrun{
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='raw', 
#'  rm.outlier=FALSE);
#' dens.abies = densform(ext.abies, climondbioclim);
#' for(i in 1:length(ext.abies[,1])){
#'   m = multiv_likelihood(ext.abies[i,6:length(ext.abies[1,])],
#'      climondbioclim, dens.abies, type = '.kde');
#'    ext.abies[i, 'prob'] = m[[1]];
#' }
#' print(ext.abies$prob)
#' }

multiv_likelihood <- function(x, clim, dens, type, w = FALSE) {
  dens.ob1 <- dens;
  varlist <- names(dens.ob1);
  varlist <- (varlist[1:((length(varlist) - 1) / 6)]);
  varlist <- base::sub(".kde", "", varlist);
  p <- vector()
  weight = vector();
  x <- as.data.frame(x)
  sumit = vector();
  centerby = vector();
  sd = vector();
  nrecs = nrow(x);
  if(w == TRUE){
  for (n in 1:length(varlist)) {
    var = varlist[[n]]
    varx <- paste(var, "x", sep = ".")
    varsd <- paste(var, "sd", sep = ".")
    varmean <- paste(var, "mean", sep = ".")
    varlook <- paste(var, type, sep = '')
    to <- max(dens.ob1[[varx]])
    from <- min(dens.ob1[[varx]])
    num = length(dens.ob1[[varx]])
    by = (to - from) / num;
    #print(by);
   # sumit[[n]] = max(dens.ob1[[varx]]*by);
    centerby[[n]] = (max(dens.ob1[[varx]]) - min(dens.ob1[[varx]]));
    sd[[n]] = dens.ob1[[varsd]]
 #   sumit[[n]] = dens.ob1[[varsd]]/(centerby[[n]]);
    sumit[[n]] = sd[[n]] / centerby[[n]];
    #sumit[[n]] = 1;
  }
    
  #  return(sumit)
    weight = 1/sumit;  
    weight = weight-(0.9*min(weight));
    weight = weight/max(weight);  
  } else {
    weight = rep(1, length(varlist));
  }
  #return(weight)
  set = vector();
  for (z in 1:nrecs){
    for (j in 1:length(varlist)) {
      var = varlist[[j]];
      varx <- paste(var, "x", sep = ".")
      varlook <- paste(var, type, sep = '')
      to <- max(dens.ob1[[varx]])
      from <- min(dens.ob1[[varx]])
      num = length(dens.ob1[[varx]])
      by = (to - from) / num; 
      search = x[, j]; 
      search = as.numeric(as.character(search));
      if (is.na(search[[z]])) {
      #  return(0)
      #  print("IS NA");
        p[j] = 0;
        next;
      }
      #print(search[[z]])
      #print(from)
      #print(by)
      bin = floor((search[[z]] - from) / by) + 1
      if (by == 0) {
        #Suggests that this is an invariant variable in this dataset
        bin = 1
      }
      #print(bin);
     # print(num);
      if (bin > num) {
        bin = num
      }
      if(bin < 1) {
        bin = 1;
      }
      lr = (dens.ob1[[varlook]][bin]); 
      if(is.na(lr)){
        lr = 0;
      }
      p[j] = (lr*by)/weight[[j]]; 
      if(is.na(p[j])){p[j] = 0;}
     # p[j] = (lr*by)/weight[[j]];
      #if(p[j] == 0){ return(list(p[[j]], j, z))}
    }
   # set[z] = prod(stats::na.omit(p));
    set[z] = sum((log(p)));
  #  if(set[z] == 0){return(list(z, search, p))}
    
  }
  return(list(set, weight));
}

#multiv_likelihood = compiler::cmpfun(multiv_likelihood);

#' Filter a set of occurrence data based on the multivariate likelihoods.
#' 
#' Given a table including localities and climate data,and a value of alpha this 
#' function will return the same table with statistical outliers
#' removed according to the alpha-confidence interval of the multivariate likelihoods.
#' 
#' @param ext_ob A data.frame of climate values (i.e., extracted from the climate raster object). Will be passed to multiv_likelihood().
#' @param clim A raster object of climate data (matching x)
#' @param dens_ob A density object from vegdistmod::densform(), vegdistmod::and_fun(), or vegdistmod::or_fun();
#' @param min The minimum likelihood to be kept. Will override the value of alpha given. Optional.
#' @param alpha The value of alpha you would like to use for confidence interval construction. Default to 0.01 for a 99 percen confidence interval.
#' @param type Designate either ".gauss" or ".kde".
#' @param w To weight log-likelihoods by coefficient of variation or not. 
#' 
#' @export
#' @examples \dontrun{
#' 
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='raw');
#' dens.abies = densform(ext.abies, climondbioclim);
#' f <- 
#'  filter_dist(ext.abies, dens.abies, 
#'  climondbioclim, alpha = 0.01, type = '.kde')
#' }

filter_dist <- function(ext_ob, dens_ob, clim, min = 0, alpha = 0.01, type = '.kde', w=FALSE) {
    
    if (type == '') {
      type = '.kde'
      
    }
    if (class(ext_ob) == 'list') {
      ext.abies.filter = list()
      for (n in 1:length(ext_ob)) {
        r <- clim[[n]]
        ext.abies = ext_ob[[n]]
        head <- which(colnames(ext.abies) %in% 'cells') + 1
        dens.ab <- dens_ob[[n]]
        for (i in 1:length(ext.abies[, 1])) {
          lookup <- (ext.abies[i, head:length(ext.abies[1,])])
          p <- multiv_likelihood(lookup, r, dens.ab, type = type, w=w)
          ext.abies[i, 'prob'] = p[[1]];
        }
        
        ext.abies.filter[[n]] <- ext.abies
      }
      isnow = 0
      probs = 0
      #Sum likelihoods across data sets
      for (z in 1:length(ext.abies.filter)) {
        probs = probs + ext.abies.filter[[z]]$prob
      }
      if (is.na(min) | min == 0) {
        min = mean(probs) + stats::qnorm(alpha / 2) * (stats::sd(probs) / sqrt(length(probs)))
      }
      filt = list()
      for (x in 1:length(ext.abies.filter)) {
        ext.abies.filter[[x]]$prob = probs
        filt[[x]] <- subset(ext.abies.filter[[x]], probs >= min)
      }
    } else {
      dens.ab <- dens_ob
      ext.abies = ext_ob
      r = clim
      head <- which(colnames(ext.abies) %in% 'cells') + 1
      for (i in 1:length(ext.abies[, 1])) {
        lookup <- (ext.abies[i, head:length(ext.abies[1, ])])
        p <- multiv_likelihood(lookup, r, dens.ab, type = type, w=w)
        ext.abies[i, 'prob'] = p[[1]];
      }
      if (min == 0) {
        min = mean(ext.abies$prob) + stats::qnorm(alpha / 2) * (stats::sd(ext.abies$prob) / sqrt(length(ext.abies$prob)))
      }
        filt <- subset(ext.abies, ext.abies$prob >= min)
    }
    return(filt)
  }

#filter_dist <- compiler::cmpfun(filter_dist);


#Hidden function to find the distance between two points
# .distance <- function(lon1, lat1, lon2, lat2) {
#   R = 6378.137
#   #pideg = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989;
# 
#  # lon1 = as.numeric(as.character(lon1));
#  # lon2 = as.numeric(as.character(lon2));
# #  lat1 = as.numeric(as.character(lat1));
#  # lat2 = as.numeric(as.character(lat2));
#   #cat(lon1, ", ", lat1, ", ", lon2, ", ", lat2, "\n")
#   toRad = pi/180;
#   lon1 = lon1 * toRad
#   
#   lon2 = lon2 * toRad
#   
#   lat1 = lat1 * toRad
#   
#   lat2 = lat2 * toRad
#   
#   dlon = lon2 - lon1
#   
#   dlat = lat2 - lat1
#   
#   a = (sin(dlat / 2) ^2) + (cos(lat1) * cos(lat2) * (sin(dlon / 2) ^2))
#   
#   d = 2 * atan2(sqrt(a), sqrt(1 - a)) * R
#   return(d)
#   
# }

#.distance = compiler::cmpfun(.distance);


#Rcpp_code <- "
#include<iostream> 
#include<cmath> 
#include <Rcpp.h>
#using namespace std;
#// [[Rcpp::export]]
#
#float distance(double lon1, double lat1, double lon2, double lat2)
#{
    # float R = 6378.137;
#
  #float toRad = 3.14159/180;
  #lon1 = lon1 * toRad;
  #lon2 = lon2 * toRad;
  #lat1 = lat1 * toRad;
  #lat2 = lat2 * toRad;
  #float dlon = lon2 - lon1;
  #float dlat = lat2 - lat1;
  
  #double a = pow(sin(dlat / 2), 2) + (cos(lat1) * cos(lat2) * pow(sin(dlon / 2),2));
  
  #double d = 2 * atan2(sqrt(a), sqrt(1 - a)) * R;
  
  #return d ;
  #}
  #
#"
     
#  Rcpp::sourceCpp(code=Rcpp_code)
 # .distance = distance;
  
Rcpp_code2 <- "
  #include<iostream> 
  #include<cmath> 
  #include <Rcpp.h>
  using namespace std;
  //' Multiply a number by two
  //'
  //' @param lon1 Starting point longitute.
  //' @param lat1 Starting point latitute.
  //' @param lon2 Ending point longitute.
  //' @param lat2 Ending point latitute.
  //' @export
  // [[Rcpp::export]]
  
  Rcpp::NumericVector findcoord(double lon, double lat, double dist, double brng) 
  {
  float R = 6378.137;
  float pi = 3.14159;
  brng = brng * (pi / 180);
  lat = lat * pi / 180;
  lon = lon * pi / 180;
  
  float lat2 = asin(sin(lat) * cos(dist / R) + cos(lat) * sin(dist / R) * cos(brng));
  float lon2 = lon + atan2(sin(brng) * sin(dist / R) * cos(lat), cos(dist / R) - sin(lat) * sin(lat2));
  lat2 = lat2 / pi * 180;
  lon2 = lon2 / pi * 180;
//  return lat2;

  Rcpp::NumericVector ll(2);
   ll[0] = lon2;
  ll[1] = lat2;
  return ll ;
  
  }
  "
  Rcpp::sourceCpp(code=Rcpp_code2)
  
  findcoord = findcoord;

  
  
#Hidden function to get coordinates given a direction and bearing from start point.
# .findcoord <- function(lon, lat, dist, brng) {
#   R = 6737.137
#  # pi = 3.14159265359;
#   
#   #lon = as.numeric(as.character(lon));
#   #lat = as.numeric(as.character(lat));
#   #dist = as.numeric(as.character(dist));
#   #brng = as.numeric(as.character(brng));
#   # print(c(lon, lat));
#   brng = brng * (pi / 180)
#   #print(brng);
#   lat <- lat * pi / 180
#   
#   lon <- lon * pi / 180
#   
#   lat2 = asin(sin(lat) * cos(dist / R) + cos(lat) * sin(dist / R) * cos(brng))
#   lon2 = lon + atan2(sin(brng) * sin(dist / R) * cos(lat),
#                      cos(dist / R) - sin(lat) * sin(lat2))
#   lat2 = lat2 / pi * 180
#   
#   lon2 = lon2 / pi * 180
#   return(c(lon2, lat2))
#   
#   
# }

#.findcoord = compiler::cmpfun(.findcoord);



#' Find nearby localities that are likely occurrences given environmental PDFs.
#' 
#' Given an occurrence set, find nearby localities that have an equal 
#' or greater multivariate likelihood given a set of environmental PDFs. 
#' Will return up to one simulated record per given occurrence, however,
#'  it is unlikely that so many suitable points will be found.
#' 
#' @param ext_ob A data.frame of climate values (i.e., extracted from the climate raster object). Will be passed to multiv_likelihood().
#' @param clim A raster object of climate data (matching x)
#' @param dens_ob A density object from vegdistmod::densform(), vegdistmod::and_fun(), or vegdistmod::or_fun();
#' @param type Designate either ".gauss" or ".kde".
#' @param name Optional. Give taxon name.
#' @param w To weight log-likelihoods by coefficient of variation or not.
#'
#' @export
#' @examples \dontrun{
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='flat', factor = 2);
#' dens.abies = densform(ext.abies, climondbioclim);
#' n <- near2(ext.abies, climondbioclim, dens.abies, type = '.kde');
#' }

near2 <- function(ext_ob, clim, dens_ob, type, name = 'NULL', w=FALSE) {
  cells = ext_ob;
  cells$prob = as.numeric(rep(0, length(cells[,1])));
  ras = clim;
  dens = dens_ob;
  if (name == 'NULL') {
    if (class(cells) == "list") {
      name = as.character(cells[[1]][1, 2])
      
    } else {
      name = as.character(cells[1, 2])
    }
    cat("NOTE: using auto name selection ", name, "\n")
  }
  if (class(cells) == "list") {
    for (n in 1:length(cells)) {
      if (colnames(cells[[n]])[ncol(cells[[n]])] == 'prob') {
        cells[[n]] = cells[[n]][, -ncol(cells[[n]])]
      }
    }
  } else {
    if (colnames(cells)[ncol(cells)] == 'prob') {
      cells <- cells[, -ncol(cells)]
    }
  }
  j = 1
  minc = 0
  maxc = 0
  maxdist = 0
  count = 0;
  newcount = 0;
  if (class(cells) == "list") {
    minc <-
      cells[[1]][grep(min(cells[[1]]$cells), cells[[1]]$cells), ]
    maxc <-
      cells[[1]][grep(max(cells[[1]]$cells), cells[[1]]$cells), ]
    maxdist <-
      distance(minc[1, 'lon'], minc[1, 'lat'], maxc[1, 'lon'], maxc[1, 'lat'])
    maxdist = 2*maxdist; #Consider calculating the actual distance matrix and selecting the maximum...
    newrecord= list();
    for(a in 1:length(cells)){
      #newrecord[[a]] <-
       # matrix(nrow = nrow(cells[[a]]),
        #       ncol = ncol(cells[[a]]))
      count = length(cells[[a]][,1]);
    }
  } else {
    minc <- cells[grep(min(cells$cells), cells$cells), ]
    maxc <-
      cells[grep(max(cells$cells), cells$cells), ]
    maxdist <-
      distance(minc[1, 'lon'], minc[1, 'lat'], maxc[1, 'lon'], maxc[1, 'lat'])
    count = length(cells[,1]);
    #newrecord = matrix(nrow = 100 * nrow(cells), ncol = ncol(cells))
  }
  #maxdist = 500; #FIXED MAXIMUM DISTANCE
 # dist.v <- (seq(2, as.integer(0.3 * maxdist))); #print(maxdist);
  dist.v <- seq(2, maxdist);
  dist.p <- (1 / sqrt(dist.v)); #print(dist.p);
 # dist.p <- ((max(dist.v) - dist.v)/max(dist.v)); #linear scaling of probabilities.
  j = 1;
  #return(list(cbind(dist.v, dist.p), maxdist))
  #dist.p <- dist.p/max(dist.p); 
  it = 1;
  origl = count;
  while(it <= count) {
    dir <- as.integer(stats::runif(1, min = 0, max = 360)); #get random bearing
    dist <- sample(dist.v,
                   size = 1,
                   replace = F,
                   prob = dist.p); #get "random" distance.
    new = 0;
    if (class(cells) == "list") {
      new = findcoord(cells[[1]][it, 'lon'], cells[[1]][it, 'lat'], dist =
                         dist, brng = dir)
      pold = vector();
      pnew = vector();
      newextr = list();
      newer = as.list(rep(NA, length(ras)));
      #print(newer);
      for (a in 1:length(cells)) {
        extr <- raster::extract(ras[[a]], cbind(new[1], new[2]), cellnumbers = T); #print("JUST EXTRACTED");print(extr);
        if(anyNA(extr) == FALSE){
          newextr[[a]] <- extr;
        } else {
          #cat("NA RECORD RETURNED\n");
          newextr[[a]] = 0;
          break;
        }
        
       #   
        #if(is.na(newextr[[a]][,1])==TRUE | length(newextr[[a]][,1])<1){next;}
       # if(length(newextr[[a]][,1]) <1){return(newextr);}
       # print(newextr[[a]])
        newer[[a]] <-
          cbind("0000", name, new[2], new[1], newextr[[a]]);
        newer[[a]] <- as.data.frame(newer[[a]]); 
        for (i in 3:length(newer[[a]][1, ])) {
          newer[[a]][, i] <- as.numeric(as.character(newer[[a]][, i]))
        }
       
        if (anyNA(newer[[a]])) {
          break;
        }
        pold[[a]] <-
          (multiv_likelihood(cells[[a]][it, 6:length(cells[[a]][1, ])], ras[[a]],  dens[[a]], type = type, w=w))[[1]]
        pnew[[a]] <-
          (multiv_likelihood(newer[[a]][1, 6:length(newer[[a]][1, ])], ras[[a]], dens[[a]], type = type, w=w))[[1]]
      }
      pnew = sum(pnew)
      pold = sum(pold)
      #return(newer);
      if(anyNA(unlist(newer)) | length(newer) < 1){  } else {
      for (b in 1:length(cells)) {
        if (pnew >= (pold)) { 
          # newrecord[j, 1:ncol(newer)] = (newer); 
          # print(newer);
          #colnames(newrecord) = colnames(cells);
          if(b == 1){
            count = count+1; 
          }
          # newrecord[j, 3] = as.numeric(as.character(newrecord[j, 3]))
          #  newrecord[j, 4] = as.numeric(as.character(newrecord[j, 4]))
          newer[[b]][1] = '0000';
          newer[[b]][2] = name;
          #   cells[count,] = newrecord[j,];
          
          cells[[b]][count,] = newer[[b]];
          #print(cells[count,]);
          cells[[b]][count,1] = '0000'
          #print(cells[count,])
          
          #cat('number of records is:', length(cells[,1]), "vs.", origl, "\n")
          j = j + 1 #I dont think I'm using j anymore.
          
          #if(length(newer[[b]])<1){next;}
         # if(anyNA(newer[[b]]) == TRUE){next;}
         # print(newer[[b]]); cat("LENGTH OF NEWER IS: ", length(newer[[b]]), "for", j, "\n");
#          newer[[b]] = stats::na.omit(newer[[b]]); #print(newer[[b]]); 
          
          #newrecord[j, 1:ncol(cells[[b]])] = (newer[[b]])
#          count = count+1;
          #return(newrecord)
          #newrecord[j, 'lat'] = as.numeric(as.character(newrecord[[b]][j, 'lat']));
          #newrecord[j, 'lon'] = as.numeric(as.character(newrecord[[b]][j, 'lon']));
          #return(newer[[b]])
          #newer[[b]]$lon = as.numeric(as.character(newer[[b]]$lon));
          #newer[[b]]$lat = as.numeric(as.character(newer[[b]]$lat));
        #  print(newer[[b]])
#          cells[[b]][count,] = as.matrix(newer[[b]]);
#         for(nex in 3:length(newer[[b]])){
 #           cells[[b]][count, nex] = as.numeric(as.character(cells[[b]][count,nex]))
            
  #        }
   #       print(j)
          #j = j+1;
          
        } else {
          
        }
      }
      }
    } else {
      new = findcoord(cells[it, 'lon'], cells[it, 'lat'], dist =
                         dist, brng = dir)
      newextr <- raster::extract(ras, cbind(new[1], new[2]), cellnumbers = T)
      newer <- cbind("0000", name, new[2], new[1], newextr)
      newer <- as.data.frame(newer)
      for (i in 3:length(newer[1, ])) {
        newer[, i] <- as.numeric(as.character(newer[, i]))
      }
      if (is.na(newer[, ncol(newer)])) {
        next
      }
      pold <-
        (multiv_likelihood(cells[it, 6:length(cells[1, ])], ras,  dens, type = type, w=w))[[1]]
      pnew <-
        (multiv_likelihood(newer[1, 6:length(newer[1, ])], ras, dens, type = type, w=w))[[1]]
      #  cat(pnew, "::", pold, "\n")
       # print(c(cells[it,c(4:3, 6:10) ], newer[1, c(4:3, 6:10)]))
      if (pnew >= (pold)) {
       # newrecord[j, 1:ncol(newer)] = (newer); 
       # print(newer);
        #colnames(newrecord) = colnames(cells);

        count = count+1; 
       # newrecord[j, 3] = as.numeric(as.character(newrecord[j, 3]))
      #  newrecord[j, 4] = as.numeric(as.character(newrecord[j, 4]))
        newer[1] = '0000';
        newer[2] = name;
     #   cells[count,] = newrecord[j,];
        
        cells[count,] = newer;
        #print(cells[count,]);
        cells[count,1] = '0000';
        #print(cells[count,])
        
        #cat('number of records is:', length(cells[,1]), "vs.", origl, "\n")
        j = j + 1
      } else {
        
      }
      #print(count); print(origl);
      if(count > 2*origl){break;}
      
    }
    it = it + 1;
  }
  return(cells);
  
  # if (class(cells) == 'list') {
  #   returnob = list();
  #   for(c in 1:length(cells)){
  #     newrecord[[c]] <- stats::na.omit(newrecord[[c]])
  #     newrecord[[c]] <- as.data.frame(newrecord[[c]])
  #     colnames(newrecord[[c]]) <- colnames(cells[[c]])
  #     returnob[[c]] <- rbind(cells[[c]], newrecord[[c]])
  #     returnob[[c]]$lon <- as.numeric(as.character(returnob[[c]]$lon))
  #     returnob[[c]]$lat <- as.numeric(as.character(returnob[[c]]$lat))
  #     for (i in 4:length(returnob[[c]][1, ])) {
  #       returnob[[c]][, i] <- as.numeric(as.character(returnob[[c]][, i]))
  #     }
  #   }
  #   return(returnob)
  #   
  # } else{
  #   newrecord <- stats::na.omit(newrecord)
  #   newrecord <- as.data.frame(newrecord)
  #   colnames(newrecord) <- colnames(cells)
  #   returnob <- rbind(cells, newrecord)
  #   returnob$lon <- as.numeric(as.character(returnob$lon))
  #   returnob$lat <- as.numeric(as.character(returnob$lat))
  #   for (i in 4:length(returnob[1, ])) {
  #     returnob[, i] <- as.numeric(as.character(returnob[, i]))
  #   }
  #   return(returnob)
  # }
  
}

#near2 <- compiler::cmpfun(near2);

#' Search for likely occurrences
#' 
#' This is a Markov Chain like search procedure to construct a set of likely occurrences 
#' for a given sample of known occurrences. Simulated occurrences are accepted if each one 
#' is more likely than the parent locality used to find the simulated point AND if the 
#' multivariate likelihood of the entire sample is improved. The improvement of the entire
#' sample is calculated after each parent locality is visited once in the search process. 
#' The process is repeated a given number of times (think Markov Chain again), to hopefully 
#' reach a near optimal set of simulated occurrences. 
#' 
#' @param ext_ob A data.frame of climate values (i.e., extracted from the climate raster object). Will be passed to multiv_likelihood().
#' @param clim A raster object of climate data (matching x)
#' @param type Designate either ".gauss" or ".kde".
#' @param manip Specify either 'reg', 'condi', or 'bayes. Default is 'condi'.
#'  See documentation for densform() for descriptions of these.
#' @param maxiter Maximum number of search iterations.
#' @param alpha The value of alpha to be used to calculate the initial confidence interval 
#'  for removing climatic outliers in the sample(s).
#' @param searchrep How many times to search for simulated localities 
#'  (per parent occurrence) per iteration. Recommend 1, but feel free to tune this parameter.
#' @param factor To be passed to the extraction() function for post search 
#'  thinning of data to limit overfitting. Set to 1 to ignore.
#' @param n An integer value of the number of bins to use for Kernel Density Estimation
#' @param w To weight the log-likelihoods by coefficient of variation or not.

#'
#' @export
#' @examples \dontrun{
#' 
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='flat', factor =4);
#' sea <- findlocal(
#'  ext.abies, climondbioclim, 
#'  type = '.kde', 
#'  maxiter = 15, searchrep = 1, 
#'  manip = 'condi');
#'  plot_clim(sea[[1]], climondbioclim[[5]])
#'  }



findlocal <-
  function(ext_ob,
           clim,
           type = '.kde',
           maxiter = 50,
           searchrep = 3,
           manip = 'condi',
           alpha = 0.05,
           factor = 4,
           n = 1024,
           w = FALSE) {
    nrec <- 100000
    
    
    ext = ext_ob
    
    if (class(ext) == 'list') {
      best = list()
    } else {
      best = ext
      
    }
    bestp = vector()
    
    searchp = vector()
    
    bc = 1
    
    nlast <- 0
    
    samecount = 0
    
    currdist <- ext
    
    
    name = currdist[1, 2]
    
    iter = 1
    pnew = 0
    r = clim
    currdist[, 1] <- as.numeric(as.character(currdist[, 1]))
    
    dens <- densform(currdist, r, manip = manip, n = n)
    vporig <- vector()
    currdist <- stats::na.omit(currdist)
    for (i in 1:length(currdist[, 1])) {
      vporig[i] <-
        (multiv_likelihood(currdist[i, 6:length(currdist[1, ])],
                           r,
                           dens,
                           type = type, w = w))[[1]]
    }
    currdist = cbind(currdist, vporig)
    colnames(currdist) = c(colnames(currdist[,-ncol(currdist)]), 'prob');
    #return(currdist)
    
    plast <- mean(vporig)
    
    if (length(vporig) >= 100) {
      vporig <- sort(vporig)
      
      origmin <- stats::quantile(vporig, probs = alpha / 2)
      
    } else {
      origmin <-
        plast + stats::qnorm(alpha / 2) * (stats::sd(vporig) / sqrt(length(vporig)))
    }
    # print(alpha)
    
    porig = plast
    print(origmin)
    last <- currdist
   # f = filter_dist(currdist, dens, r, min = origmin, type = type)
    
   # currdist <- f
    
    
    currdist = subset(currdist, currdist$prob >= origmin);
    
    dens <- densform(currdist, r, manip = manip, n=n)
    vporig <- vector()
    currdist <- stats::na.omit(currdist)
    for (i in 1:length(currdist[, 1])) {
      vporig[i] <-
        (multiv_likelihood(currdist[i, 6:length(currdist[1, ])],
                           r,
                           dens,
                           type = type, w = w))[[1]]
    }
    currdist$prob = vporig
    #return(currdist)
    
    plast <- mean(vporig)
    porig = plast
    origmin <- min(vporig)
    
    print(origmin)
    cat("porig is: ", porig, '\n')
    
    last <- currdist
    best = last
    
    best.dens = dens
    
    
    orig_num = length(currdist[, 1])
    
    
    
    bestp[bc] = porig
    
    searchp[bc] = porig
    
    bc = bc + 1
    print(bc)
    
    while (porig < 0) {
      last <- currdist
      currdist = best
      
      
      for (i in 1:searchrep) {
        currdist <- near2(currdist, r, dens, name = name, type = type)
      }
      
      
      vporig <- vector()
      bep <- vector()
      
      currdist <- stats::na.omit(currdist)
      if (sum(colnames(currdist) %in% 'vporig') > 0) {
        currdist = currdist[, -ncol(currdist)]
        
      }
      sim <- subset(currdist, currdist[, 1] == "0000")
      
      if (length(sim[, 1]) > 10) {
        # print("spThin");
        print(length(sim[, 1]))
        ext.sim <-
          extraction(
            sim[, 1:(which(colnames(sim) == 'cells') - 1)],
            clim,
            schema =
              'flat',
            factor = factor,
            rm.outlier = FALSE
          )
        
        sim <- ext.sim
        
        sim[, 1] = rep("0000", length(sim[, 1]))
        
      } else {
        sim <- sim[!duplicated(sim[, "cells"]), ]
        
        
      }
      
      sub <- subset(currdist, currdist[, 1] != '0000')
     # sub <- sub[,1:(ncol(sub)-1)]
    #  sim <- sim[,1:(ncol(sim)-1)]
      
      currdist <- rbind(sub, sim)
      
      dens = densform(currdist, r,  manip = manip, n=n)
      for (i in 1:length(currdist[, 1])) {
        vporig[i] <-
          (multiv_likelihood(currdist[i, 6:length(sub[1, ])],
                             r,
                             dens,
                             type = type, w = w))[[1]]
      }
      currdist = cbind(currdist, vporig)
      
      currdist <- as.data.frame(currdist)
      
      
      colnames(currdist) = c(colnames(currdist[, -ncol(currdist)]), 'prob')
      #return(currdist)
      # print(cbind(currdist[,1], vporig));
      for (i in 1:length(best[, 1])) {
        bep[i] <-
          (multiv_likelihood(best[i, 6:length(sub[1, ])],
                             r,
                             dens,
                             type = type, w = w))[[1]]
      }
      # return(cbind(sub, vporig))
      
      # return( vporig)
      origmin <- min(vporig[which(currdist[, 1] != '0000')])
      print(origmin)
      # f = subset(currdist, vporig>= origmin); # This should not be necessary bc near2 only picks points better than the occ set.
      
      
      # currdist <- f
      
      p <- vector()
      
      nrec = length(currdist[1, ])
      pnew <- mean(vporig)
      
      iter = iter + 1
      cat('Like: pnew', pnew, ' vs. pprimary', mean(bep), '\n')
      if (is.na(pnew)) {
       # return(currdist)
        break;
      } ##IF this line triggers something went very wrong
      if (mean(bep) < pnew) {
        best  <- currdist
        best.dens <- dens
        samecount = 0
        
        bestp[bc] = pnew
        
        searchp[bc] = pnew
        
        bc = bc + 1
        #return(list(best, best.dens));
      } else {
        turbo = TRUE; #Allow user to set in function
        if(turbo == TRUE){
          sim <- subset(currdist, currdist[,1] == '0000');
          sub <- subset(currdist, currdist[,1] != '0000')
          sim <- subset(sim, sim$prob >= mean(bep))
          currdist <- rbind(sub,sim);
          if(mean(currdist[,ncol(currdist)])>mean(bep)){
            best  <- currdist
            best.dens <- dens
            samecount = 0
            
            bestp[bc] = mean(currdist[,ncol(currdist)])
            
            searchp[bc] = mean(currdist[,ncol(currdist)])
            
            bc = bc + 1
          } else {
            bestp[bc] = bestp[bc - 1]
            currdist = best
            searchp[bc] = pnew
            bc  = bc + 1
          }
        } else {
          bestp[bc] = bestp[bc - 1]
          currdist = best
          searchp[bc] = pnew
          bc  = bc + 1
        }
      }
      #THis is not how we want to count this anymore. 
      #No progress now is if there is no improvement. I think this works still, but is not right.
      if (bestp[[bc - 2]] == bestp[[bc - 1]]) {
        samecount = samecount + 1
        
      }
      plast = pnew
      print(iter)
      if (samecount > 10) {
        break
      }
      
      if (iter >= maxiter) {
        break
      }
    }
    return(list(best, bestp, searchp))
  }

#findlocal <- compiler::cmpfun(findlocal);


#' Search for likely simulated localities with geographic subsetting.
#' 
#' This wraps the search() function within a set of commands that randomly splits the distribution 
#' in to 4 quadrants 'n' number of times. For widespread, genetically diverse species this has 
#' the effect of applying more locally relevant probability functions to the data.
#' 
#' @param ext_ob A data.frame of climate values (i.e., extracted from the climate raster object). Will be passed to multiv_likelihood().
#' @param clim A raster object of climate data (matching x)
#' @param type Designate either ".gauss" or ".kde".
#' @param manip Specify either 'reg', 'condi', or 'bayes. Default is 'condi'.
#'  See documentation for densform() for descriptions of these.
#' @param maxiter Maximum number of search iterations.
#' @param searchrep How many times to search for simulated localities 
#'  (per parent occurrence) per iteration. Recommend 1, but feel free to tune this parameter.
#' @param alpha The value of alpha to be used to calculate the initial confidence interval 
#'  for removing climatic outliers in the sample(s).
#' @param divisions How many times should the data be split into quadrants? 
#'  Default is 5 resulting in 20 (5x4quads) geographically oriented samples to be selected.
#' @param factor To be passed to the extraction() function for post search 
#'  thinning of data to limit overfitting. Set to 1 to ignore.
#' @param parallel True or False to use parallel computing. This implements 
#'  each division as embarassingly parallel processes. However,
#'  note that this option changes the meaning of the divisions object 
#'  to the maximum number of division iterations. Fewer iterations may 
#'  be returned if too few occurrences are selected.
#' @param nclus If parallel is TRUE then set the maximum number of cores 
#'  to use in the compute cluster. Default is 2.
#' @param w To weight log-likelihoods by the coefficient of variation or not.
#' 
#' @export
#' @examples \dontrun{
#' data(abies);
#' head(abies)
#' ext.abies = extraction(abies, climondbioclim, 
#'  schema='flat', factor =4,  
#'  rm.outlier=FALSE);
#' sea <- geo_findlocal(ext.abies, climondbioclim, 
#'  type = '.kde', maxiter = 5, 
#'  searchrep = 1, manip = 'condi', 
#'  divisions = 8, parallel =TRUE, nclus = 4)
#'  plot_clim(sea, climondbioclim[[5]])
#'  }

 

geo_findlocal <-
  function(ext_ob,
           clim,
           type,
           maxiter = 10,
           searchrep = 1,
           manip = 'condi',
           alpha = 0,
           divisions = 10,
           factor = 4,
           parallel = FALSE,
           nclus = 2,
           w = FALSE) {
    ext = ext_ob
    
    if (w == TRUE) {
    }
    
    
    # search = list();
    if (parallel == TRUE) {
      cl <- parallel::makeCluster(nclus, type = "SOCK")
      # doParallel::registerDoParallel(cl)
      doSNOW::registerDoSNOW(cl)
      
      
      search <-
        foreach::foreach(i = 1:divisions,
                         .combine = 'rbind',
                         .packages = 'vegdistmod') %dopar% {
                           source('~/Desktop/cracle_testing/vegdistmod/R/search_fun.R')
                           
                           n = 0
                           
                           sub.nw = NA
                           
                           sub.ne = NA
                           
                           sub.se = NA
                           
                           sub.sw = NA
                           
                           while (n < 1) {
                             sam = ext[sample(nrow(ext), 1), ]
                            
                             sam.lat <- sam$lat
                             
                             sam.lon <- sam$lon
                             
                             sub.nw <-
                               subset(ext, ext$lon <= sam.lon &
                                        ext$lat >= sam.lat)
                             sub.ne <-
                               subset(ext, ext$lon >= sam.lon &
                                        ext$lat >= sam.lat)
                             sub.sw <-
                               subset(ext, ext$lon <= sam.lon &
                                        ext$lat <= sam.lat)
                             sub.se <-
                               subset(ext, ext$lon >= sam.lon &
                                        ext$lat <= sam.lat)
                             if (length(sub.nw[, 1]) < 5 |
                                 length(sub.ne[, 1]) < 5 |
                                 length(sub.se[, 1]) < 5 | 
                                 length(sub.sw[, 1]) < 5) {
                               
                             } else{
                               n = 1
                               
                             }
                           }
                         #  return(sub.se);
                           searchit = matrix(nrow = 0, ncol = ncol(ext) + 1)
                           
                           searchit <- as.data.frame(searchit)
                           
                           colnames(searchit) = c(colnames(ext), 'prob')
                           if (length(sub.sw[, 1]) >= 5) {
                             search.sw <-
                               findlocal(
                                 sub.sw,
                                 clim,
                                 type,
                                 maxiter = maxiter,
                                 searchrep = searchrep,
                                 manip = manip,
                                 alpha = alpha,
                                 w = w
                               )
                             searchit = rbind(searchit, search.sw[[1]])

                           }
                           
                           
                           if (length(sub.nw[, 1]) >= 5) {
                             search.nw <-
                               findlocal(
                                 sub.nw,
                                 clim,
                                 type,
                                 maxiter = maxiter,
                                 searchrep = searchrep,
                                 manip = manip,
                                 alpha = alpha,
                                 w = w
                               )
                             searchit = rbind(searchit, search.nw[[1]])
                             
                             
                           }
                           
                           
                           if (length(sub.ne[, 1]) >= 5) {
                             search.ne <-
                               findlocal(
                                 sub.ne,
                                 clim,
                                 type,
                                 maxiter = maxiter,
                                 searchrep = searchrep,
                                 manip = manip,
                                 alpha = alpha,
                                 w = w
                               )
                             searchit = rbind(searchit, search.ne[[1]])
                             
                             
                           }
                           
                           
                           if (length(sub.se[, 1]) >= 5) {
                             search.se <-
                               findlocal(
                                 sub.se,
                                 clim,
                                 type,
                                 maxiter = maxiter,
                                 searchrep = searchrep,
                                 manip = manip,
                                 alpha = alpha,
                                 w = w
                               )
                             searchit = rbind(searchit, search.se[[1]])
                             
                             
                           }
                           
                           # return(search.se);
                           
                           
                           
                  
                           
                           #  to.search = rbind(stats::na.omit(search.nw[[1]]), stats::na.omit(search.se[[1]]),
                           # stats::na.omit(search.sw[[1]]), stats::na.omit(search.ne[[1]]))
                           
                           return(stats::na.omit(searchit))
                         }
      
      
      
      parallel::stopCluster(cl)
      
      return(search)
      
    } else {
      search = list()
      
      i = 0
      
      while (i < divisions) {
        print(i)
        
        n = 0
        
        sub.nw = NA
        
        sub.ne = NA
        
        sub.se = NA
        
        sub.sw = NA
        
        while (n < 1) {
          sam = ext[sample(nrow(ext), 1), ]
          sam.lat <- sam$lat
          
          sam.lon <- sam$lon
          
          sub.nw <-
            subset(ext, ext$lon <= sam.lon & ext$lat >= sam.lat)
          sub.ne <-
            subset(ext, ext$lon >= sam.lon & ext$lat >= sam.lat)
          sub.sw <-
            subset(ext, ext$lon <= sam.lon & ext$lat <= sam.lat)
          sub.se <-
            subset(ext, ext$lon >= sam.lon & ext$lat <= sam.lat)
          if (length(sub.nw[, 1]) < 10 |
              length(sub.ne[, 1]) < 10 |
              length(sub.se[, 1]) < 10 | length(sub.sw[, 1]) < 10) {
            
          } else{
            n = 1
            
          }
        }
        search.nw = NA
        
        search.sw = NA
        
        search.ne = NA
        
        search.se = NA
        
        searchit = matrix(nrow = 0, ncol = ncol(ext) + 1)
        
        searchit <- as.data.frame(searchit)
        
        colnames(searchit) = c(colnames(ext), 'prob')
        #return(searchit)
        if (length(sub.nw[, 1]) >= 5) {
          search.nw <-
            findlocal(
              sub.nw,
              clim,
              type,
              maxiter = maxiter,
              searchrep = searchrep,
              manip = manip,
              alpha = alpha,
              w = w
            )
          # return(search.nw)
          # return(list(search.nw, searchit))
          colnames(search.nw[[1]]) = colnames(searchit)
          searchit = rbind(stats::na.omit(search.nw[[1]]), searchit)
          
          
        }
        
        if (length(sub.sw[, 1]) >= 5) {
          search.sw <-
            findlocal(
              sub.sw,
              clim,
              type,
              maxiter = maxiter,
              searchrep = searchrep,
              manip = manip,
              alpha = alpha,
              w = w
            )
          colnames(search.sw[[1]]) = colnames(searchit)
          searchit = rbind(stats::na.omit(search.sw[[1]]), searchit)
          
        }
        
        if (length(sub.ne[, 1]) >= 5) {
          search.ne <-
            findlocal(
              sub.ne,
              clim,
              type,
              maxiter = maxiter,
              searchrep = searchrep,
              manip = manip,
              alpha = alpha,
              w = w
            )
          colnames(search.ne[[1]]) = colnames(searchit)
          searchit = rbind(stats::na.omit(search.ne[[1]]), searchit)
          
        }
        
        if (length(sub.se[, 1]) >= 5) {
          search.se <-
            findlocal(
              sub.se,
              clim,
              type,
              maxiter = maxiter,
              searchrep = searchrep,
              manip = manip,
              alpha = alpha,
              w = w
            )
          colnames(search.se[[1]]) = colnames(searchit)
          searchit = rbind(stats::na.omit(search.se[[1]]), searchit)
          
        }
        
        i = i + 1
        
        search[[i]] = (searchit)
        
        #search[[i]] = rbind(stats::na.omit(search.nw[[1]]), stats::na.omit(search.se[[1]]),
        #                 stats::na.omit(search.sw[[1]], stats::na.omit(search.ne[[1]])));
        
      }
      
      
      
      
    }
    #  hold = matrix(ncol = ncol(search[[1]]));
    # return(search)
    hold = search[[1]]
    
    
    for (n in 2:length(search)) {
    #  print(n)
      
      hold = rbind(hold, search[[n]])
      
    }
    return(hold)
    
  }

#geo_findlocal = compiler::cmpfun(geo_findlocal);

#' Plot lat/long points on a raster map.
#' 
#' Plotting with fancy colors. Canned so you don't have to think too hard about it.
#' 
#' @param ext_ob A data.frame of climate values (i.e., extracted from the climate raster object). MUST include columns named 'lon' and 'lat'.
#' @param clim A raster object of climate data (matching ext_ob)
#' @param boundaries A shapefile (e.g., GADM country or state outlines)
#' @param file If a file path is set this function will try to write the plot as a png to that path
#' @param col Color of points to plot.
#' @param legend TRUE or FALSE to plot the legend on the map.
#' @param l.cex cex parameter to pass to legend function
#'
#' @export
#' @examples \dontrun{
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='raw');
#' plot_clim(ext.abies, climondbioclim[[5]]);
#' }
#' 
plot_clim <- function(ext_ob, clim, boundaries ='', file='', col = 'red', legend = TRUE, l.cex = 0.9) {
  #require(RColorBrewer);
  #require(classInt);
  poi = ext_ob;
  usa <- boundaries;
  nclr = 8;
  breaks <- round((raster::maxValue(clim) - raster::minValue(clim))/nclr, digits = 3)
  plotclr = ( grDevices::topo.colors(1000));
  plotvar <- seq(raster::minValue(clim), raster::maxValue(clim), by = breaks);
  class <- classInt::classIntervals(plotvar, nclr, style = 'fixed', fixedBreaks = plotvar, 1)
  colcode <- classInt::findColours(class, plotclr);

  if(file != ''){
    grDevices::png(
      file,
      units = 'in',
      height = 4,
      width = 5,
      res = 400
    )
  }
  graphics::par(mai = c(0.5, 0.5, 0.5, 0))
  
  raster::plot(
    clim,
    #main = names(clim),
    #col = rev(rainbow(1000, start = 0, end = 0.7))
    col = (grDevices::topo.colors(1000)), colNA = 'black',
    
    breaks = seq(raster::minValue(clim), raster::maxValue(clim), length.out = 1000),
    legend = F,
    # legend.width = 1,
    axis.args = list(
      at = seq(raster::minValue(clim), raster::maxValue(clim), by = 100),
      labels = round(seq(raster::minValue(clim), raster::maxValue(clim), by = 100), 0),
      cex.axis = 0.9
    )
  )
  points(poi$lon,
         poi$lat,
         col = col,
         pch = 15,
         cex = 0.5)
  if(boundaries!= ''){
    graphics::plot(usa, add = T)
  }
  if(legend == TRUE){
  graphics::legend(
    "topleft",
    legend = names(attr(colcode, "table")),
    title = names(clim),
    text.col = 'white',
    fill = attr(colcode, "palette"),
    cex = l.cex,
    bty = 'n'
  )
  }
  
  
  if(file!=''){
    grDevices::dev.off()
  }
  
  
}




#' Convert log-likelihoods of any vegdistmod PDF to a raster heat-map.
#' 
#' This is not an ENM. What it does do is provides a geographic representation of the likelihood functions
#'  that are used in CRACLE and the spatial defragmentation functions. 
#' 
#' @param clim A raster object of climate data (matching ext_ob)
#' @param dens A vegdistmod density object (see densform())
#' @param parallel Make use of multicore architecture
#' @param nclus If parallel is TRUE, how many cores should be allocated.
#' @param type Which PDF should be used, .gauss or .kde
#' @param w TRUE or FALSE should variable PDFs be weighted by relative niche breadth.
#' @author Robert Harbert, \email{rharbert@amnh.org}
#' @author Avery Hill
#'
#' @export
#' @examples \dontrun{
#' data(abies);
#' data(climondbioclim);
#' ext.abies = extraction(abies, climondbioclim, schema='raw', rm.outlier=TRUE, alpha = 0.005);
#' dens <- densform(ext.abies, climondbioclim, manip = 'condi', 
#'                  kern = 'gaussian', n = 128, bg.n = 1000)
#' h = heat_up(climondbioclim, dens, parallel=FALSE, type = '.kde', nclus =4)
#' hs = sum(h)
#' ex.h = raster::extract(hs, ext.abies[,4:3])
#' plot(hs>sort(ex.h)[ceiling(0.01*length(ex.h))])
#' points(ext.abies[,4:3], col ='green')
#'  
#' ##Bootstrap: With train/test subsetting and model evaluation
#'
#' binary = list();
#' bin.auc = list();
#' ev.auc = list();
#' data.ex = extraction(abies, climondbioclim, schema='flat', factor =2, rm.outlier=TRUE, alpha = 0.01)
#' for (i in 1:100){
#' pick = as.numeric(sample(data.ex[,1], 0.5*length(data.ex[,1]), replace =F));
#' train = data.ex[which(data.ex[,1] %in% pick),]
#' test = data.ex[-which(data.ex[,1] %in% pick),]
#' d.train <- densform(train, climondbioclim, 
#'      manip = 'condi', kern = 'gaussian', n = 128, bg.n = 1000)
#' h.t = heat_up(climondbioclim, d.train, parallel=TRUE, type = '.kde', nclus =4)
#' hs.t = sum(h.t)
#' ex.h = raster::extract(hs.t, test[,4:3])
#' binary[[i]] = hs.t>sort(ex.h)[ceiling(0.1*length(ex.h))];
#' plot(binary[[i]])
#' points(test[,4:3], col ='purple', pch = 20)
#' 
#' bg <- rad_bg(test[,4:3], climondbioclim[[1]], radius = 2000, n = 200)
#' bg.e <- raster::extract(hs.t, bg[,4:3]);
#' ev <- evaluate(ex.h, bg.e)
#' print(ev)
#' ev.auc[[i]] = ev@auc;
#' 
#' bg.bin <- raster::extract(binary[[i]], bg[,4:3]);
#' ex.bin <- raster::extract(binary[[i]], test[,4:3]);
#' ev.bin <- evaluate(ex.bin, bg.bin)
#' print(ev.bin)
#' bin.auc[[i]] = ev.bin@auc
#' }
#' bin.stack = stack(unlist(binary))
#' bin.sum = sum(bin.stack)
#' bin.weightave = sum(bin.stack * unlist(bin.auc))/sum(unlist(bin.auc))
#' plot(bin.sum>50); #bootstrap consensus
#' plot(bin.weightave>0.5); #weighted (bin.auc) consensus
#' points(data.ex[,4:3], pch = 20, cex =0.5)
#' 
#' 
#' }


heat_up <- function(clim, dens, parallel = FALSE, nclus =4, type = '.kde', w = FALSE){
  #whole = .get_bg(clim);
  whole.ex=raster::extract(clim,raster::extent(clim),cellnumbers=T,df=T) #climate values for climate raster
  
  #whole.ex = extraction(whole, clim, schema='raw')
  #cells = which(colnames(whole.ex)=='cells')+1;
  if(parallel ==TRUE){
    
    cl <- parallel::makeCluster(nclus, type = "SOCK")
    doSNOW::registerDoSNOW(cl);
    npart = ceiling(nlayers(clim)/nclus)
    lvec <-
      foreach::foreach(i = 1:nlayers(clim),
              .packages = 'vegdistmod', .combine='cbind'
                    ) %dopar% {
#        source('~/Desktop/cracle_testing/vegdistmod/R/cracle_build.R')
        vp = .vecprob(whole.ex[,i+2],
                      dens[paste(names(clim[[i]]), 'x', sep = ".")][[1]], 
                      dens[paste(names(clim[[i]]), type, sep ='')][[1]])
            whole.ex[,i+2] = (vp); #replace values with probabilities
          
        return(as.numeric(whole.ex[,i+2]));
      }
    parallel::stopCluster(cl)
    
    whole.ex = lvec;
    for(z in 1:nlayers(clim)){
      r=raster::raster(nrows=nrow(clim),ncol=ncol(clim), crs="+proj=longlat +datum=WGS84",ext=raster::extent(clim)) #make raster of same extent and dimensions as original climate raster
      r=raster::setValues(r,values=log(whole.ex[,z])) #set the values as the new probability values
      clim[[z]] = r; #replace raster layers from clim with prob layers.
    }
    names(clim) = colnames(whole.ex)
    return(clim)
    
   # lvec = unlist(lvec);
    #return(lvec)
  } else {
    
    for(i in 3:ncol(whole.ex)){
      m = .vecprob(whole.ex[,i], 
                            dens[paste(names(clim[[i-2]]), 'x', sep = ".")][[1]], 
                            dens[paste(names(clim[[i-2]]), type, sep ='')][[1]])
      whole.ex[,i] = m; #replace values with probabilities
    }
  
  
  ###extract all cells from raster
  
  # get multiv_likelihood for each cell
  
  # Write vector of likelihoods to raster of same extent:
  for(z in 1:nlayers(clim)){
    r=raster::raster(nrows=nrow(clim),ncol=ncol(clim), crs="+proj=longlat +datum=WGS84",ext=raster::extent(clim)) #make raster of same extent and dimensions as original climate raster
    r=raster::setValues(r,values=log(whole.ex[,z+2])) #set the values as the new probability values
    clim[[z]] = r; #replace raster layers from clim with prob layers.
  }
  names(clim) = colnames(whole.ex[,3:ncol(whole.ex)])
  return(clim)
  }
  
}

#heat_up <- compiler::cmpfun(heat_up);


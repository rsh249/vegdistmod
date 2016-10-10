#' @import classInt
#' @import grDevices
#' @import doSNOW
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
  curl_string = "curl --user rsh249:Roanoke1999 http://cloud.diversityoflife.org/cgi-div/tmp_mat_get.pl?taxon="
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
  e <- extent(clim)
  uplon <- e[1]
  uplat <- e[4]
  nrow = d[1]
  ncol = d[2]
  xres = (e[2] - e[1])/ncol;
  yres = (e[4] - e[3])/nrow;
  

  r <- res(clim)
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
    #dim <- dim(clim[[1]]);
    #n = dim[1]*dim[2];
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
    
    
    #   dim <- dim(clim);
    #  n = dim[1]*dim[2];
#    bg <- randomPoints(r, n)
#    bg <- cbind(rep('0000', length(bg[,1])), rep('bg', length(bg[,1])), bg[,2], bg[,1])
#    colnames(bg) <- c('ind_id', 'tax', 'lat', 'lon');
#    bg = as.data.frame(bg);
#   bg$lon <- as.numeric(as.character(bg$lon));
#    bg$lat <- as.numeric(as.character(bg$lat));
#    ret = bg;
  }
  #return(ret);
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
#'
#' @export
#' @examples
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='raw');
#' dens.abies = densform(ext.abies, climondbioclim);
#' for(i in 1:length(ext.abies[,1])){
#'   m = multiv_likelihood(ext.abies[i,6:length(ext.abies[1,])],
#'      climondbioclim, dens.abies, type = '.kde');
#'    ext.abies[i, 'prob'] = m;
#' }
#' print(ext.abies$prob)

multiv_likelihood <- function(x, clim, dens, type) {
  dens.ob1 <- dens;
  varlist <- names(dens.ob1);
  varlist <- (varlist[1:((length(varlist) - 1) / 5)]);
  varlist <- sub(".kde", "", varlist);
  p <- vector()
  weight = vector();
  x <- as.data.frame(x)
  sumit = vector();
  for (n in 1:length(varlist)) {
    var = varlist[[n]]
    varx <- paste(var, "x", sep = ".")
    varlook <- paste(var, type, sep = '')
    to <- max(dens.ob1[[varx]])
    from <- min(dens.ob1[[varx]])
    num = length(dens.ob1[[varx]])
    by = (to - from) / num;
    sumit[[n]] = 1;
  }
  weight = 1/(sumit/max(sumit));
  for (j in 1:length(varlist)) {
    var = varlist[[j]]
    varx <- paste(var, "x", sep = ".")
    varlook <- paste(var, type, sep = '')
    to <- max(dens.ob1[[varx]])
    from <- min(dens.ob1[[varx]])
    num = length(dens.ob1[[varx]])
    by = (to - from) / num
    search = x[, j]
    search = as.numeric(as.character(search));
    if (is.na(search)) {
      return(0)
    }
    bin = round((search - from) / by) + 1
    if (by == 0) {
      #Suggests that this is an invariant variable in this dataset
      bin = 1
    }
    if (bin > num) {
      bin = num
    }
    lr = dens.ob1[[varlook]][bin];
    p[j] = log(lr*by);
  }
  return(list(sum(p), weight));
}

#multiv_likelihood = compiler::cmpfun(multiv_likelihood);

#' Filter a set of occurrence data based on the multivariate likelihoods.
#' 
#' Given a table including localities and climate data,and a value of alpha this 
#' function will return the same table with statistical outliers
#' removed according to the alpha-confidence interval of the likelihoods.
#' 
#' @param ext_ob A data.frame of climate values (i.e., extracted from the climate raster object). Will be passed to multiv_likelihood().
#' @param clim A raster object of climate data (matching x)
#' @param dens_ob A density object from vegdistmod::densform(), vegdistmod::and_fun(), or vegdistmod::or_fun();
#' @param min The minimum likelihood to be kept. Will override the value of alpha given. Optional.
#' @param alpha The value of alpha you would like to use for confidence interval construction. Default to 0.01 for a 99 percen confidence interval.
#' @param type Designate either ".gauss" or ".kde".
#' 
#' @export
#' @examples
#' 
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='raw');
#' dens.abies = densform(ext.abies, climondbioclim);
#' f <- 
#'  filter_dist(ext.abies, dens.abies, climondbioclim, alpha = 0.01, type = '.kde')
#' 

filter_dist <- function(ext_ob, dens_ob, clim, min = 0, alpha = 0.01, type = '.kde') {
    
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
          p <- multiv_likelihood(lookup, r, dens.ab, type = type)
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
        p <- multiv_likelihood(lookup, r, dens.ab, type = type)
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


#HIdden function to find the distance between two points
.distance <- function(lon1, lat1, lon2, lat2) {
  R = 6737
  lon1 = as.numeric(as.character(lon1));
  lon2 = as.numeric(as.character(lon2));
  lat1 = as.numeric(as.character(lat1));
  lat2 = as.numeric(as.character(lat2));
  
  lon1 = lon1 * pi / 180
  
  lon2 = lon2 * pi / 180
  
  lat1 = lat1 * pi / 180
  
  lat2 = lat2 * pi / 180
  
  dlon = lon2 - lon1
  
  dlat = lat2 - lat1
  
  a = (sin(dlat / 2)) ^ 2 + cos(lat1) * cos(lat2) * (sin(dlon / 2)) ^ 2
  
  c = 2 * atan2(sqrt(a), sqrt(1 - a))
  
  d = R * c
  
  return(d)
  
}

#.distance = compiler::cmpfun(.distance);


#Hidden function to get coordinates given a direction and bearing from start point.
.findcoord <- function(lon, lat, dist, brng) {
  R = 6737
  lon = as.numeric(as.character(lon));
  lat = as.numeric(as.character(lat));
  dist = as.numeric(as.character(dist));
  brng = as.numeric(as.character(brng));
  # print(c(lon, lat));
  brng = brng * (pi / 180)
  #print(brng);
  lat <- lat * pi / 180
  
  lon <- lon * pi / 180
  
  lat2 = asin(sin(lat) * cos(dist / R) + cos(lat) * sin(dist / R) * cos(brng))
  lon2 = lon + atan2(sin(brng) * sin(dist / R) * cos(lat),
                     cos(dist / R) - sin(lat) * sin(lat2))
  lat2 = lat2 / pi * 180
  
  lon2 = lon2 / pi * 180
  return(c(lon2, lat2))
  
  
}

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
#'
#' @export
#' @examples
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='raw');
#' dens.abies = densform(ext.abies, climondbioclim);
#' n <- near2(ext.abies, climondbioclim, dens.abies, type = '.kde');

near2 <- function(ext_ob, clim, dens_ob, type, name = 'NULL') {
  cells = ext_ob;
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
      .distance(minc[1, 'lon'], minc[1, 'lat'], maxc[1, 'lon'], maxc[1, 'lat'])
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
      .distance(minc[1, 'lon'], minc[1, 'lat'], maxc[1, 'lon'], maxc[1, 'lat'])
    count = length(cells[,1]);
    #newrecord = matrix(nrow = 100 * nrow(cells), ncol = ncol(cells))
  }
  #maxdist = 500; #FIXED MAXIMUM DISTANCE
  #dist.v <- (seq(2, as.integer(0.3 * maxdist))); #print(maxdist);
  dist.v <- seq(5, maxdist);
  #dist.p <- (1 / sqrt(dist.v)); #print(dist.p);
  dist.p <- ((max(dist.v) - dist.v)/max(dist.v)); #linear scaling of probabilities.
  j = 1;
  
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
      new = .findcoord(cells[[1]][it, 'lon'], cells[[1]][it, 'lat'], dist =
                         dist, brng = dir)
      pold = vector();
      pnew = vector();
      newextr = list();
      newer = as.list(rep(NA, length(ras)));
      #print(newer);
      for (a in 1:length(cells)) {
        extr <- extract(ras[[a]], cbind(new[1], new[2]), cellnumbers = T); #print("JUST EXTRACTED");print(extr);
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
          (multiv_likelihood(cells[[a]][it, 6:length(cells[[a]][1, ])], ras[[a]],  dens[[a]], type = type))[[1]]
        pnew[[a]] <-
          (multiv_likelihood(newer[[a]][1, 6:length(newer[[a]][1, ])], ras[[a]], dens[[a]], type = type))[[1]]
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
      new = .findcoord(cells[it, 'lon'], cells[it, 'lat'], dist =
                         dist, brng = dir)
      newextr <- extract(ras, cbind(new[1], new[2]), cellnumbers = T)
      newer <- cbind("0000", name, new[2], new[1], newextr)
      newer <- as.data.frame(newer)
      for (i in 3:length(newer[1, ])) {
        newer[, i] <- as.numeric(as.character(newer[, i]))
      }
      if (is.na(newer[, ncol(newer)])) {
        next
      }
      pold <-
        (multiv_likelihood(cells[it, 6:length(cells[1, ])], ras,  dens, type = type))[[1]]
      pnew <-
        (multiv_likelihood(newer[1, 6:length(newer[1, ])], ras, dens, type = type))[[1]]
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
#' @param bg An object of background point climate data matching the 
#'  output of extraction(). Generate random backround points and then use extraction().

#'
#' @export
#' @examples
#' 
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='raw');
#' dens.abies = densform(ext.abies, climondbioclim);
#' sea <- findlocal(
#'  ext.abies, climondbioclim, 
#'  type = '.kde', 
#'  maxiter = 5, searchrep = 1, 
#'  manip = 'condi');

findlocal <-
  function(ext_ob,
           clim,
           bg = 0,
           type,
           maxiter = 50,
           searchrep = 3,
           manip = 'condi',
           alpha = 0.05,
           factor = 4) {
    nrec <- 100000
    
    
    ext = ext_ob
    
    if (length(bg) < 2) {
      if (class(ext) == 'list') {
        bg.hold = .get_bg(clim[[1]])
        
        bg = list()
        
        for (i in 1:length(clim)) {
          bg[[i]] = extraction(bg.hold, clim[[i]], schema = 'raw')
          
        }
        
        
      } else {
        bg.e = .get_bg(clim)
        #bg.e = get_bg(clim)
        
        bg = extraction(bg.e, clim, schema = 'raw')
        
      }
      
    }
    if (class(ext) == 'list') {
      best = list()
    } else {
      best = 0
    }
    bestp = vector()
    
    searchp = vector()
    
    bc = 1
    
    nlast <- 0
    
    samecount = 0
    
    currdist <- ext
    
    if (class(ext) == 'list') {
      name = currdist[[1]][1, 2]
      
    } else{
      name = currdist[1, 2]
    }
    iter = 1
    pnew = 0
    r = clim
    if (class(ext) == 'list') {
      dens = list()
      
      vporig = matrix(nrow = length(ext[[1]][, 1]), ncol = length(ext))
      for (a in 1:length(ext)) {
        currdist[[a]][,1] = as.numeric(as.character(currdist[[a]][,1]));
        currdist[[a]] <- stats::na.omit(currdist[[a]])
        dens[[a]] <-
          densform(currdist[[a]], r[[a]], bg = bg[[a]], manip = manip)
        currdist[[a]] <- stats::na.omit(currdist[[a]])
        for (i in 1:length(currdist[[a]][, 1])) {
          vporig[i, a] <-
            (multiv_likelihood(currdist[[a]][i, 6:length(currdist[[a]][1,])],
                               r[[a]],
                               dens[[a]],
                               type = type))[[1]]
        }
      }
      vporig = apply(vporig, 1, sum)
      
    } else {
      currdist[,1] <- as.numeric(as.character(currdist[,1])); 
      dens <- densform(currdist, r, bg = bg, manip = manip)
      vporig <- vector()
      currdist <- stats::na.omit(currdist)
      for (i in 1:length(currdist[, 1])) {
        vporig[i] <-
          (multiv_likelihood(currdist[i, 6:length(currdist[1,])],
                             r,
                             dens,
                             type = type))[[1]]
      }
    }
    plast <- mean(vporig)
    
    if (length(vporig) >= 100) {
      vporig <- sort(vporig)
      
      origmin <- stats::quantile(vporig, probs = alpha / 2)
      
    } else {
      origmin <-
        plast + stats::qnorm(alpha / 2) * (stats::sd(vporig) / sqrt(length(vporig)))
    }
    print(alpha)
    
    porig = plast
    print(origmin)
    last <- currdist
    f = filter_dist(currdist, dens, r, min = origmin, type = type)
    
    currdist <- f
    if (class(ext) == 'list') {
      dens = list()
      
      vporig = matrix(nrow = length(currdist[[1]][, 1]), ncol = length(ext))
      for (a in 1:length(ext)) {
        currdist[[a]] <- stats::na.omit(currdist[[a]])
        dens[[a]] <-
          densform(currdist[[a]], r[[a]], bg = bg[[a]], manip = manip)
        currdist[[a]] <- stats::na.omit(currdist[[a]])
        for (i in 1:length(currdist[[a]][, 1])) {
          vporig[i, a] <-
            (multiv_likelihood(currdist[[a]][i, 6:length(currdist[[a]][1,])],
                               r[[a]],
                               dens[[a]],
                               type = type))[[1]]
        }
      }
      vporig <- stats::na.omit(vporig)
      vporig = apply(vporig, 1, sum)
      
      for (d in 1:length(currdist)) {
        currdist[[d]]$prob = vporig
        
      }
    } else {
      dens <- densform(currdist, r, bg = bg, manip = manip)
      vporig <- vector()
      currdist <- stats::na.omit(currdist)
      for (i in 1:length(currdist[, 1])) {
        vporig[i] <-
          (multiv_likelihood(currdist[i, 6:length(currdist[1,])],
                             r,
                             dens,
                             type = type))[[1]]
      }
      currdist$prob = vporig
      
    }
    plast <- mean(vporig)
    porig = plast
    origmin <- min(vporig)
    
    print(origmin)
    cat("porig is: ", porig, '\n')
    
    last <- currdist
    best = last
    
    best.dens = dens
    
    if (class(currdist) == 'list') {
      orig_num = length(currdist[[1]][, 1])
    } else {
      orig_num = length(currdist[, 1])
      
      
    }
    
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
      
      if (class(ext) == 'list') {
        dens = list()
        
        
        
        ##DEFAULT TO BE SET?
        #Spatial filter is to be done on the first environmental data layer set.
        #User should have the choice, or this behavior needs to be documented
        
        sim <- subset(currdist[[1]], currdist[[1]][, 1] == "0000")
        if(length(sim[,1])>5){
        ext.sim <-
          extraction(sim[, 1:(which(colnames(sim) == 'cells') - 1)], clim[[1]], schema = 'flat', factor = factor)
        
        sim = ext.sim
        
        sim[, 1] = rep("0000", length(sim[, 1]))
        } else {
          sim = sim[!duplicated(sim[,"cells"]),];

        }        
        sub = subset(currdist[[1]], currdist[[1]][, 1] != '0000')
        
        currdist[[1]] = rbind(sub, sim)
        
        
        
        for (z in 2:length(clim)) {
          #problem here if there are fewer than 5 new records the extraction() function returns nothing.
          ext.this = extraction(sim[, 1:(which(colnames(sim) == 'cells') - 1)], clim[[z]], schema = 'raw')
          
          ext.this[, 1] = rep("0000", length(ext.this[, 1]))
          
          sub = subset(currdist[[z]], currdist[[z]][, 1] != '0000')
          #print(length(sub[,1]));
          currdist[[z]] = rbind(sub, ext.this)
          
        }
        sub = 0
        
        vporig = matrix(nrow = length(currdist[[1]][, 1]),
                        ncol = length(currdist))
        
        for (a in 1:length(currdist)) {
          currdist[[a]] <- stats::na.omit(currdist[[a]])
          
          sub <- subset(currdist[[a]], currdist[[a]][, 1] != "0000")
          
          dens[[a]] = densform(currdist[[a]], r[[a]], bg = bg[[a]], manip = manip)
          
          for (i in 1:length(currdist[[a]][, 1])) {
            vporig[i, a] <-
              (multiv_likelihood(sub[i, 6:length(sub[1,])],
                                 r[[a]],
                                 dens[[a]],
                                 type = type))[[1]]
            
          }
          
        }
        vporig = apply(vporig, 1, sum)
        
      } else {
        vporig <- vector()
        currdist <- stats::na.omit(currdist)
        sim <- subset(currdist, currdist[, 1] == "0000")
        
        if (length(sim[, 1]) > 5) {
          print("spThin");
          ext.sim <-
            extraction(sim[, 1:(which(colnames(sim) == 'cells') - 1)], clim, schema =
                         'flat', factor = factor)
          
          sim <- ext.sim
          
          sim[, 1] = rep("0000", length(sim[, 1]))
          
        } else {
          sim <- sim[!duplicated(sim[,"cells"]),];
          
        }
        
        sub <- subset(currdist, currdist[, 1] != '0000')
        
        currdist <- rbind(sub, sim)
        
        
        dens = densform(currdist, r, bg = bg, manip = manip)
        
        for (i in 1:length(currdist[, 1])) {
          vporig[i] <-
            (multiv_likelihood(sub[i, 6:length(sub[1,])],
                               r,
                               dens,
                               type = type))[[1]]
        }
      } 
      origmin <- min(vporig)
      print(origmin)
      f = filter_dist(currdist, dens, r, min = origmin, type = type)
      
      
      currdist <- f
      p <- vector()
      if (class(currdist) == 'list') {
        nrec = length(currdist[[1]][1,])
        pnew <- mean(currdist[[1]][, ncol(currdist[[1]])])
      } else {
        nrec = length(currdist[1,])
        pnew <- mean(currdist[, ncol(currdist)])
      }
      iter = iter + 1
      cat('Like: ', pnew, ' ', bestp[bc - 1], '\n')

      if (bestp[[bc - 1]] < pnew) {
        best  <- currdist
        best.dens <- dens
        samecount = 0;
        bestp[bc] = pnew
        
        searchp[bc] = pnew
        
        bc = bc + 1
        
      } else {
        bestp[bc] = bestp[bc - 1]
        
        currdist = best
        
        
        searchp[bc] = pnew
        
        bc  = bc + 1
        
      }
      if (bestp[[bc - 2]] == bestp[[bc-1]]) {
        samecount = samecount + 1
        
      }
      plast = pnew
      print(iter)
      if(samecount >10){break;}
      
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
#' @param bg An object of background point climate data matching the 
#'  output of extraction(). Generate random backround points and then use extraction().
#' @param parallel True or False to use parallel computing. This implements 
#'  each division as embarassingly parallel processes. However,
#'  note that this option changes the meaning of the divisions object 
#'  to the maximum number of division iterations. Fewer iterations may 
#'  be returned if too few occurrences are selected.
#' @param nclus If parallel is TRUE then set the maximum number of cores 
#'  to use in the compute cluster. Default is 2.
#' @export
#' @examples
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='flat', factor = 4);
#' dens.abies = densform(ext.abies, climondbioclim);
#' sea <- geo_findlocal(ext.abies, climondbioclim, 
#'  type = '.kde', maxiter = 5, 
#'  searchrep = 1, manip = 'condi', 
#'  divisions = 2)

 

geo_findlocal <- function(ext_ob, clim, type, maxiter = 10, bg = 0, searchrep = 1, manip = 'condi', alpha = 0.05, divisions = 10, factor = 4, parallel = FALSE, nclus = 2){
  ext = ext_ob;
  if(length(bg)<2){
    if(class(ext)=='list'){
      bg.hold = .get_bg(clim[[1]]);
      bg = list();
      for(i in 1:length(clim)){
        bg[[i]]= extraction(bg.hold, clim[[i]], schema='raw');
      }
      
      
    } else {
      bg.e = .get_bg(clim)
      #bg.e = get_bg(clim);
      bg = extraction(bg.e, clim, schema='raw');
    }
    
  }
  
 # search = list();
  if (parallel == TRUE) {
    nn = 1
    
    cl <- parallel::makeCluster(nclus, type = "SOCK")
   # doParallel::registerDoParallel(cl)
    doSNOW::registerDoSNOW(cl);
    
    search <-
      foreach(i = 1:divisions,
              .combine = 'rbind',
              .packages = 'vegdistmod') %dopar% {
                ##Try the "when()" function instead of while()
                #when(nn <= divisions)
                #while(i<divisions){
                #print(i)
                #if (nn >= divisions) {return(NA)}
                if (class(ext) == 'list') {
                  # for(zzz in 1:length(ext)){
                  sam.lat <- sample(ext[[1]]$lat, 1)
                  
                  sam.lon <- sample(ext[[1]]$lon, 1)
                  
                  search.nw.l <- list()
                  
                  search.ne.l <- list()
                  
                  search.sw.l <- list()
                  
                  search.se.l <- list()
                  
                  for (zzz in 1:length(ext)) {
                    sub.nw <-
                      subset(ext[[zzz]], ext[[zzz]]$lon <= sam.lon &
                               ext[[zzz]]$lat >= sam.lat)
                    
                    sub.ne <-
                      subset(ext[[zzz]], ext[[zzz]]$lon >= sam.lon &
                               ext[[zzz]]$lat >= sam.lat)
                    sub.sw <-
                      subset(ext[[zzz]], ext[[zzz]]$lon <= sam.lon &
                               ext[[zzz]]$lat <= sam.lat)
                    
                    sub.se <-
                      subset(ext[[zzz]], ext[[zzz]]$lon >= sam.lon &
                               ext[[zzz]]$lat <= sam.lat)
                    
                    if (length(sub.nw[, 1]) <= 5) {
                      return(NA)
                    }
                    
                    if (length(sub.sw[, 1]) <= 5) {
                      return(NA)
                    }
                    
                    if (length(sub.ne[, 1]) <= 5) {
                      return(NA)
                    }
                    
                    if (length(sub.se[, 1]) <= 5) {
                      return(NA)
                    }
                    search.ne.l[[zzz]] <- sub.nw
                    
                    search.nw.l[[zzz]] <- sub.ne
                    
                    search.se.l[[zzz]] <- sub.se
                    
                    search.sw.l[[zzz]] <- sub.sw
                    
                  }
                  nn = nn + 1
                  
                  print(length(search.ne.l))
                  
                  search.ne <-
                    findlocal(
                      search.ne.l,
                      clim,
                      type,
                      maxiter = maxiter,
                      searchrep = searchrep,
                      manip = manip,
                      factor = factor
                    )
                  search.sw <-
                    findlocal(
                      search.sw.l,
                      clim,
                      type,
                      maxiter = maxiter,
                      searchrep = searchrep,
                      manip = manip,
                      factor = factor
                    )
                  search.nw <-
                    findlocal(
                      search.nw.l,
                      clim,
                      type,
                      maxiter = maxiter,
                      searchrep = searchrep,
                      manip = manip,
                      factor = factor
                    )
                  search.se <-
                    findlocal(
                      search.se.l,
                      clim,
                      type,
                      maxiter = maxiter,
                      searchrep = searchrep,
                      manip = manip,
                      factor = factor
                    )
                  
                  to.search = rbind(search.nw, search.se, search.sw, search.ne)
                  
                  return(to.search)
                  #}
                } else {
                  sam.lat <- sample(ext$lat, 1)
                  
                  sam.lon <- sample(ext$lon, 1)
                  
                  sub.nw <-
                    subset(ext, ext$lon <= sam.lon & ext$lat >= sam.lat)
                  
                  sub.ne <-
                    subset(ext, ext$lon >= sam.lon & ext$lat >= sam.lat)
                  sub.sw <-
                    subset(ext, ext$lon <= sam.lon & ext$lat <= sam.lat)
                  
                  sub.se <-
                    subset(ext, ext$lon >= sam.lon & ext$lat <= sam.lat)
                  search.ne = NA
                  search.se = NA
                  search.sw = NA
                  search.nw = NA
                  if (length(sub.nw[, 1]) >= 5) {
                    search.nw <-
                      findlocal(
                        sub.nw,
                        clim,
                        type,
                        maxiter = maxiter,
                        bg = bg,
                        searchrep = searchrep,
                        manip = manip,
                        alpha = alpha
                      )
                  }
                  
                  if (length(sub.sw[, 1]) >= 5) {
                    search.sw <-
                      findlocal(
                        sub.sw,
                        clim,
                        type,
                        maxiter = maxiter,
                        bg = bg,
                        searchrep = searchrep,
                        manip = manip,
                        alpha = alpha
                      )
                  }
                  
                  if (length(sub.ne[, 1]) >= 5) {
                    search.ne <-
                      findlocal(
                        sub.ne,
                        clim,
                        type,
                        maxiter = maxiter,
                        bg = bg,
                        searchrep = searchrep,
                        manip = manip,
                        alpha = alpha
                      )
                  }
                  
                  if (length(sub.se[, 1]) >= 5) {
                    search.se <-
                      findlocal(
                        sub.se,
                        clim,
                        type,
                        maxiter = maxiter,
                        bg = bg,
                        searchrep = searchrep,
                        manip = manip,
                        alpha = alpha
                      )
                  }
                  
                  
                
                
                 
                
                  nn = nn + 1
                  
                  
                  to.search = rbind(search.nw[[1]], search.se[[1]], search.sw[[1]], search.ne[[1]])
                  
                  return(to.search)
                }
                
                
              }
    parallel::stopCluster(cl)
    
    return(search)
    
  } else {
    search = list();
    i=0;
    while(i<divisions){
      print(i)
      if(class(ext)=='list'){
        # for(zzz in 1:length(ext)){
        sam.lat <- sample(ext[[1]]$lat, 1);
        sam.lon <- sample(ext[[1]]$lon, 1);
        search.nw.l <- list();
        search.ne.l <- list();
        search.sw.l <- list();
        search.se.l <- list();
        for(zzz in 1:length(ext)){
          sub.nw <- subset(ext[[zzz]], ext[[zzz]]$lon <= sam.lon & ext[[zzz]]$lat >= sam.lat);
          sub.ne <- subset(ext[[zzz]], ext[[zzz]]$lon >= sam.lon & ext[[zzz]]$lat >= sam.lat)
          sub.sw <- subset(ext[[zzz]], ext[[zzz]]$lon <= sam.lon & ext[[zzz]]$lat <= sam.lat);
          sub.se <- subset(ext[[zzz]], ext[[zzz]]$lon >= sam.lon & ext[[zzz]]$lat <= sam.lat); 
          if(length(sub.nw[,1])>= 5){
            search.nw.l[[zzz]] <- sub.ne;
          };
          if(length(sub.sw[,1])>= 5){
            search.sw.l[[zzz]] <- sub.sw;
          };
          if(length(sub.ne[,1])>=5){
            search.ne.l[[zzz]] <- sub.nw;
          };
          if(length(sub.se[,1])>=5){
            search.se.l[[zzz]] <- sub.se;
          };
        }
        i=i+1;
        print(length(search.ne.l));
        search.ne <- findlocal(search.ne.l, clim, type, maxiter=maxiter, searchrep=searchrep, manip = manip, factor = factor)
        search.sw <- findlocal(search.sw.l, clim, type, maxiter=maxiter, searchrep=searchrep, manip = manip, factor = factor)
        search.nw <- findlocal(search.nw.l, clim, type, maxiter=maxiter, searchrep=searchrep, manip = manip, factor = factor)
        search.se <- findlocal(search.se.l, clim, type, maxiter=maxiter, searchrep=searchrep, manip = manip, factor = factor)
        
        search[[i]] = rbind(search.nw, search.se, search.sw, search.ne);
      } else {
        sam.lat <- sample(ext$lat, 1);
        sam.lon <- sample(ext$lon, 1);
        sub.nw <- subset(ext, ext$lon <= sam.lon & ext$lat >= sam.lat);
        sub.ne <- subset(ext, ext$lon >= sam.lon & ext$lat >= sam.lat)
        sub.sw <- subset(ext, ext$lon <= sam.lon & ext$lat <= sam.lat);
        sub.se <- subset(ext, ext$lon >= sam.lon & ext$lat <= sam.lat); 
        search.nw = NA;
        search.sw = NA;
        search.ne = NA;
        search.se = NA;
        if(length(sub.nw[,1])>= 5){
          search.nw <- findlocal(sub.nw, clim, type, maxiter=maxiter,  bg=bg, searchrep=searchrep, manip = manip, alpha = alpha)
          
        };
        if(length(sub.sw[,1])>= 5){
          search.sw <- findlocal(sub.sw, clim, type, maxiter=maxiter, bg=bg,  searchrep=searchrep, manip = manip, alpha = alpha)
          
        };
        if(length(sub.ne[,1])>=5){
          search.ne <- findlocal(sub.ne, clim, type, maxiter=maxiter, bg=bg, searchrep=searchrep, manip = manip, alpha = alpha)
        };
        if(length(sub.se[,1])>=5){
          search.se <- findlocal(sub.se, clim, type, maxiter=maxiter,  bg=bg, searchrep=searchrep, manip = manip, alpha = alpha)
          
        };
        i=i+1;
        
        search[[i]] = rbind(search.nw[[1]], search.se[[1]], search.sw[[1]], search.ne[[1]]);
        
      }
      
      
    } 
    
  }
  hold = matrix(ncol = ncol(search[[1]]));
  hold = search[[1]];

  for(n in 2:length(search)){
    print(n)

    hold = rbind(hold, search[[n]])
    
  }
  return(hold);
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
#'
#' @export
#' @examples
#' data(abies);
#' ext.abies = extraction(abies, climondbioclim, schema='raw');
#' plot_clim(ext.abies, climondbioclim[[5]]);
#' 
plot_clim <- function(ext_ob, clim, boundaries ='', file='', col = 'red') {
  #require(RColorBrewer);
  #require(classInt);
  poi = ext_ob;
  usa <- boundaries;
  nclr = 8;
  breaks <- round((maxValue(clim) - minValue(clim))/nclr, 0)
  plotclr = ( grDevices::topo.colors(1000));
  plotvar <- seq(minValue(clim), maxValue(clim), by = breaks);
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
  graphics::par(mai = c(0.5, 0.5, 0.5, 0), bg = 'darkgrey')
  
  raster::plot(
    clim,
    #main = names(clim),
    #col = rev(rainbow(1000, start = 0, end = 0.7))
    col = (grDevices::topo.colors(1000)),
    
    breaks = seq(minValue(clim), maxValue(clim), length.out = 1000),
    legend = F,
    # legend.width = 1,
    axis.args = list(
      at = seq(minValue(clim), maxValue(clim), by = 100),
      lebels = round(seq(minValue(clim), maxValue(clim), by = 100), 0),
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
  graphics::legend(
    "topleft",
    legend = names(attr(colcode, "table")),
    title = names(clim),
    text.col = 'white',
    fill = attr(colcode, "palette"),
    cex = 0.5,
    bty = 'n'
  )
  
  
  if(file!=''){
    grDevices::dev.off()
  }
  
  
}


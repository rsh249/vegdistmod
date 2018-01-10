#' @import grDevices
#' @import dismo
#' @import doParallel
#' @import foreach
#' @import iterators
#' @import doSNOW
#' @import raster
NULL

#' Extract environmental data
#' 
#' This function is a feature-added wrapper for raster::extract();
#' @param data Distribution data. A data.frame that should include (at least) a column for species/taxon name named 'tax', latitude named 'lat', longitude named 'lon', and an optional column 'sub' if it is necessary to define subgroups (i.e., if 'tax' corresponds to genera but sampling needs to know which records belong to which species. See param 'schema' below). 
#' @param clim A raster object (see raster::raster() and raster::stack() documentation for reading raster files into R).
#' @param schema A string of value "raw", "flat", or "species" to define the sampling protocol. In "raw", all records are counted (including duplicate exact localities). In "flat", all unique localities will be counted where a unique locality is defined as a raster grid cell. Under the "flat" sampling strategy two records in the same raster grid cell will be counted as one. The option "species", only applies when taxa are identified as genera and species identities are represented in the "sub" column of the data object. In "species", each unique locality is counted for each species within the group (taxon). This weighs more diverse localities higher. Default is "raw".
#' @param factor An integer value for the methods "flat" and "spec" to increase the systematic sampling grid size to courser resolutions than the given climate grid. The value of factor corresponds to the number of rows and columns to aggregate into each courser grid cell. Default is 0 which will not be processed.
#' @param rm.outlier TRUE or FALSE. Indicate whether to remove points that are climatic outliers for at least one variable given a normal 95 percent confidence interval.
#' @param alpha Confidence level (i.e., 0.05) for clipping out outlier records.
#' @param nmin Minimum number of records allowed. Taxa or groups with fewer records will not be returned.
#' @export
#' @examples
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' extr.flat = extraction(data=distr, clim= climondbioclim, schema='flat');
#' extr.spec = extraction(data=distr, clim= climondbioclim, schema='species');
extraction <- function(data, clim, schema = "raw", factor = 0, rm.outlier = FALSE,  alpha = 0.01, nmin = 5){

	if(length(data[,1]) < 5){cat('ERR: Too few records\n'); return(NA);}

	mat.larr <- data;
	phytoclim <- clim;
  #nclat <- which(colnames(mat.larr)=='lat');
	#nclon <- which(colnames(mat.larr)=='lon');
    
    #    if(parallel==FALSE){
        extr.larr <- raster::extract(phytoclim, cbind(mat.larr$lon, mat.larr$lat), cellnumbers=T);
        # } else {
        #  bloc = round(nrow(mat.larr)/nclus);
        #cl <- parallel::makeCluster(nclus, type = "SOCK")
        #doSNOW::registerDoSNOW(cl);
        
        # extr.larr <-
        #foreach::foreach(i = 1:nclus,
        #    .packages = 'vegdistmod',
        #    .combine = 'rbind') %dopar% {
        #        start = i;
        #       end = start+bloc;
        #        subs = mat.larr[i:(i+bloc-1),];
        #        e = raster::extract(phytoclim, cbind(subs$lon, subs$lat), cellnumbers=T);
        #        return(e);
        #    }
            
            #        parallel::stopCluster(cl)

            
        
        #}
    ##
    if(nrow(stats::na.omit(extr.larr))<5){
	  cat("ERR: Records out of study area\n")
	  return(NA)
	}
	extr.larr <- cbind(mat.larr, extr.larr);
	if(schema != 'raw'){
		if(factor == 0){} else {
			r2 <- raster::aggregate(phytoclim, fact = factor, fun=mean);
			tmp.ext <- raster::extract(r2, cbind(mat.larr$lon, mat.larr$lat), cellnumbers=T);
			extr.larr[,(ncol(extr.larr)+1)] = extr.larr[,'cells'];
			extr.larr[,'cells'] = tmp.ext[,'cells'];
			
		}
		
		
	}
	extr.larr <- stats::na.omit(extr.larr)
	if(schema == "raw"){
		holder <- data.frame();
		tlist <- unique(extr.larr$tax);
		for(i in 1:length(tlist)){
			set <- subset(extr.larr, extr.larr$tax == tlist[i]);
		
			#if(length(set[,1])>=5){
				holder <- rbind(holder, set);
				#print(length(holder[,1]))
				
			#}			
		}
		extr.larr = holder;
	} else {
		holder <- data.frame();
		tlist <- unique(extr.larr$tax);
		for(i in 1:length(tlist)){
			set <- subset(extr.larr, extr.larr$tax == tlist[i]);
			if(schema == "flat"){
				sub = set
				sub <- sub[!duplicated(sub[,"cells"]),];
				#if(length(sub[,1])>=5){
					holder <- rbind(holder, sub);
			#	}	
			}
			if(schema == "species"){
				glist <- unique(set$sub);
				for(n in 1:length(glist)){
					sub <- subset(set, set$sub == glist[n]);
					sub <- sub[!duplicated(sub[,"cells"]),];
					if(length(sub[,1])>=5){
						holder <- rbind(holder, sub);
					}
				}		
			}
		}
		extr.larr <- holder; 
	}
	if(schema != 'raw'){

		extr.larr[,'cells'] = extr.larr[,ncol(extr.larr)];

		extr.larr = extr.larr[,-ncol(extr.larr)];
	}
  #print("EXTRACTION MONITOR:")
	#  print(length(holder[,1]));
	
	#print(length(extr.larr[,1]));
	head = which(colnames(extr.larr)=='cells')-1;
	#print(head)
	
  extr.larr[,1] = as.numeric(as.character(extr.larr[,1]))
  if(rm.outlier== TRUE){
    
    for(nn in 1:raster::nlayers(phytoclim)){
      n.mean <- mean(as.numeric(extr.larr[,(head+nn)]));
      n.sd <- stats::sd(as.numeric(extr.larr[,(head+nn)]));
      rn <- length(extr.larr[,(head+nn)]);
      t = stats::qt((1-(alpha/2)), rn-1);
      minci = n.mean-(t*n.sd);
      maxci = n.mean+(t*n.sd);
      extr.larr <- subset(extr.larr, extr.larr[,(head+nn)] >= minci);
      extr.larr <- subset(extr.larr, extr.larr[,(head+nn)] <= maxci);
      
    }
    
  }
  t.list = unique(extr.larr$tax);
  
  if(length(t.list)>1){
  hold = data.frame();

  for(zz in 1:length(t.list)){
    sub <- subset(extr.larr, extr.larr$tax == t.list[[zz]]);
    if(nrow(sub) < nmin){
      
    } else {
      hold = rbind(hold, sub);
    }
  }
  colnames(hold) = colnames(extr.larr);
  } else {
    hold = extr.larr;
  }
	
	return(hold);
};

#extraction = compiler::cmpfun(extraction);



#' Generate standard probability density functions for each taxon/variable
#' 
#' This function takes extracted climate data (from an object generated by the vegdistmod::extraction() function) for one taxon/species and generates probability density functions for variable using both a Gaussian (normal) approximation and a Gaussian Kernel Density estimator.
#' @param ex An object derived from the extraction() function.
#' @param clim A raster object (see raster::raster() and raster::stack() documentation for reading raster files into R).
#' @param name A character string describing (preferably) the group for which PDFs are being constructed (i.e., a species binomial). If none is supplied, a value of column "tax" is selected as a default.
#' @param bw A bandwidth compatible with stats::density(). Options include "nrd", "nrd0", "ucv", "bcv", etc.. Default (and recommended) value is "nrd0".
#' @param kern Type of Kernel to smooth with. Recommend 'gaussian', 'optcosine', or 'epanechnikov'. See: stats::density for options.
#' @param n Number of equally spaced points at which the probability density is to be estimated. Defaults to 1024. A lower number increases speed but decreases resolution in the function. A higher number increases resolution at the cost of speed. Recommended values: 512, 1024, 2048, ....
#' @param manip Character string of 'reg' for intersectional likelihood, 'condi' for conditional likelihood statement.
#' @param from vector of starting points by variable. Default is the variable layer minimum.
#' @param to vector of ending points by variable. Default is the variable layer maximum.
#' @param clip A character string of value "range" or "95conf" or "99conf". Should the probability functions be clipped to either the empirical range or the 95 or 99 percent confidence interval? 
#' @param bg.n If manip = 'condi'. How many background points PER OCCURRENCE record should be sampled. Default is 1000.
#' @param bg Optionally send a matrix of an extraction object to use as the background sample.
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' extr.sub = subset(extr.raw, extr.raw$tax == extr.raw[5,'tax']);
#' dens.sub = densform(extr.sub, clim = climondbioclim, bw = 'nrd0', n = 128, bg.n=25);
#' densplot(dens.sub, names(climondbioclim[[1]]));
#' }

densform <- function(ex, clim, 
                     name = '', bw = "nrd0", kern = 'gaussian',
                     
                     manip = 'condi', n = 1024, 
                     from = 0, to = 0, clip = 0,
                     bg.n = 10, bg=NULL){
#  kern = 'gaussian'
  cut = 0;
 # adjust = (512*60)/n;
  adjust = 1;
  condi = FALSE;
  if(manip == 'condi') {
    condi = TRUE; #print("Conditional Likelihood")
  }

	data = ex;
	if(name == ''){
		name = data[2,'tax'];
	};
	pi = 3.14159265359;
	extr.larr <- data;
	head = which(colnames(ex) %in% 'cells') - 1;
	phytoclim <- clim;
		larr.den <- data.frame();
 		larr.den.x <- data.frame();
		larr.den.gauss <- data.frame();
		larr.mean <- data.frame();
		larr.sd <- data.frame();
	  larr.w <- data.frame();
		eval <- data.frame();
		bg.eval = data.frame();
	  ncoords = length(extr.larr$lon);
	if(condi==TRUE){ 
	dmatrix = matrix(ncol = ncoords,
	                 nrow = ncoords);
	#  print("Getting distance matrix");
	if(is.null(bg)){
	for(xx in 1:ncoords){
	  for(yy in xx:ncoords){
	   # dmatrix[xx,yy] <- vegdistmod:::.distance(extr.larr$lon[xx], extr.larr$lat[xx], 
	       #                         extr.larr$lon[yy], extr.larr$lat[yy]);
	    dmatrix[xx,yy] <- distance(extr.larr$lon[xx], extr.larr$lat[xx], 
	                                                                  extr.larr$lon[yy], extr.larr$lat[yy]);
	    
	    
	  }
	}
	#NOTE: The background is selected from a radius around each occurrence
	#record within 5x the mean distance between all points in the sample.
	#This threshold is arbitrary and needs to be empirically tested, but does seem to work.
	bg.rad = 3*max(stats::na.omit(dmatrix));
#	print(bg.rad)
	#  bg.rad = max(stats::na.omit(dmatrix));
		 # print("before rad_bg");
	bg.ex <- rad_bg(cbind(extr.larr$lon, extr.larr$lat), 
	                phytoclim, 
	                radius = bg.rad, 
	                n = bg.n)
	} else {
	  bg.ex = bg;
	}
	#	 print("after")
	}
	

		for(i in 1:length(names(phytoclim))){	
	  	fr <- raster::minValue(phytoclim[[i]]); ### At issue with from and to is reproducibility. If defined globally from the raster then it will always be compatible.
			t <- raster::maxValue(phytoclim[[i]]);
    	if(condi == TRUE){
     	  bg.vec = bg.ex[,names(phytoclim[[i]])];
    	  
	    	bg.den <- stats::density(c(as.numeric(bg.vec), as.numeric(extr.larr[,names(phytoclim[[i]])])), 
			                           n = n, kernel = kern, 
			                           adjust = 1, from = fr, 
			                           to = t, bw = bw, 
			                           na.rm = TRUE); 
	    	o.vec <- as.numeric(extr.larr[,names(phytoclim[[i]])]);
			  weights = 1/.vecprob(o.vec, bg.den$x, bg.den$y)
        weights = weights/sum(stats::na.omit(weights)); #So the weights sum to 1. Could consider other scaling. i.e., by rank order
       # print(weights);
        o.vec <- cbind(o.vec, weights);
        o.vec = stats::na.omit(o.vec);
       # if(sum(o.vec[,2]) != 1){print(names(phytoclim[[i]]));print(o.vec);}
        
        #For the KDE PDF use the weights option in the stats::density function to estimate the conditional probability
			  den <- stats::density(o.vec[,1], 
			                        n = n, kernel = kern, adjust = 1.2,
			                        from = fr,  to = t, weights = o.vec[,2],
			                        bw = bw, na.rm = TRUE);
		  	bg.mean <- mean(as.numeric(bg.ex[,names(phytoclim[[i]])]));
		  	bg.sd <- stats::sd(bg.ex[,names(phytoclim[[i]])]);
        x = extr.larr[,names(phytoclim[[i]])]
        
        #For the gaussian PDF estimates use the weighted mean and sd:
        mean <- stats::weighted.mean(as.numeric(o.vec[,1]), as.numeric(o.vec[,2]));
        sd <- sqrt(sum(weights * (x - mean)^2))
        rn <- length(extr.larr[,names(phytoclim[[i]])]);
        
        if(sd == 0 || is.na(sd) == "TRUE"){
          sd = 0.01;
        };
        for(num in 1:length(den$x)){
          eval[num,1] <- ((1/(sqrt((2*pi)
                                   *(sd^2)))
                           *(2.71828^(-1*((den$x[num] - mean)^2)
                                      /(2*sd^2)))));
          
        };
        
			} else {
			  den <- stats::density(as.numeric(extr.larr[,names(phytoclim[[i]])]), 
			                        n = n, kernel = kern, 
			                        from = fr,  to = t, 
			                        bw = bw, na.rm = TRUE);
			  mean <- mean(extr.larr[,names(phytoclim[[i]])]);
			  sd <- stats::sd(extr.larr[,i+head+1]);
			  rn <- length(extr.larr[,names(phytoclim[[i]])]);
			  
			  if(sd == 0 || is.na(sd) == "TRUE"){
			    sd = 0.01;
			  };
			  for(num in 1:length(den$x)){
			    eval[num,1] <- ((1/(sqrt((2*pi)*(sd^2)))*(2.71828^(-1*((den$x[num] - mean)^2)/(2*sd^2)))));
			  };
			}
	  	
	  	##Below is for clipping the PDFs to fit empirical ranges or inferred confidence intervals
			minci = 0;maxci=0;
			if(clip == "95conf"){
			  t = stats::qt(0.975, rn-1);
			  minci = mean-(t*sd);
			  maxci = mean+(t*sd);
			}
			if(clip == "99conf"){

			  t = stats::qt(0.995, rn-1);
			  minci = mean-(t*sd);
			  maxci = mean+(t*sd);
			}
			if(clip == "999conf"){
			  
			  t = stats::qt(0.9995, rn-1);
			  minci = mean-(t*sd);
			  maxci = mean+(t*sd);
			}
			if(clip == 'range'){
			  minci = min(extr.larr[,names(phytoclim[[i]])]);
			  maxci = max(extr.larr[,names(phytoclim[[i]])]);
			}
			if(minci == 0 & maxci==0){} else{
			  
			  if(den$x[[1]] < minci){
			    larr.den[1, i] = 0;
			    eval[1,1] =0;
			  }
			  if(den$x[[1]] > maxci){
			    larr.den[1, i] = 0;
			    eval[1,1] = 0;
			  }
			  
		  	for(zz in 2:length(den$x)){
			  
	  		  if(den$x[[zz]] < minci){
	  		    larr.den[(zz-1), i] = 0;
			      eval[(zz-1),1] =0;
			    }
		  	  if(den$x[[zz]] > maxci){
		  	    larr.den[(zz-1), i] = 0;
			      eval[(zz-1),1] = 0;
		  	  }
			  }
			}
			larr.den.x[1:n, i] <- den$x;
			larr.den[1:n, i] <- den$y;
			
			larr.den.gauss[1:n, i] <- eval[,1];
			larr.mean[1,i] <- mean;
			larr.sd[1,i] <- sd;
			range <- max(larr.den.x[,i]) - min(larr.den.x[,i]);
			#print(sd); print(range);
			w <- sd/range;
			w = 1/w;
			larr.w[1,i] = w;
	    #weight = as.numeric(larr.w[1,])
	   # weight =  weight - (min(weight));
	   # weight = weight/(0.5*max(weight));
	    #larr.w[1,] = weight;
		}
	  larr.w[1,] = larr.w[1,]/(0.5*max(larr.w[,1]));
		colnames(larr.den.gauss) <- c(paste(names(phytoclim), "gauss", sep = "."));
		colnames(larr.mean) <- c(paste(names(phytoclim), "mean", sep = "."));
		colnames(larr.sd) <- c(paste(names(phytoclim), "sd", sep = "."));
	  colnames(larr.w) <- c(paste(names(phytoclim), "w", sep = "."));
	
		colnames(larr.den) <- c(paste(names(phytoclim), "kde", sep = "."));
		colnames(larr.den.x) <- c(paste(names(phytoclim), "x", sep = "."));
		name = data.frame(name);
		larr.mean = data.frame(larr.mean);
		larr.sd = data.frame(larr.sd);
		colnames(name) <- "name";
		fin <- c(larr.den, larr.den.x, larr.den.gauss, larr.mean, larr.sd, larr.w, name);
		fin <- .makeaucone(fin);
		return(fin);
	
};

#densform = compiler::cmpfun(densform);

.vecprob <- function(search, x,y){
  to <- max(x); 
  from <- min(x);
  num <- length(x);
  by = (to - from)/num; 
  bin = floor((search - from) / by)+1; 
  
  if(length(bin)>1){
    for(nn in 1:length(bin)){
      if (is.na(bin[[nn]])) { bin[[nn]] = 1; }
      if (bin[[nn]] <= 1) { bin[[nn]] = 1; }
      if (bin[[nn]] > num) { bin[[nn]] = num; }
    }
  } else{
    if (bin <= 1){bin =1;}
    if (bin > num){
      bin = num
    }
  }
  ret = y[bin]*by;

 # ret[ret == NA] <- min(stats::na.omit(ret));
  ret[ret== -Inf] <- NA;
  ret[is.na(search)] = NA;
  return(ret);
}

#' Generate backround data within a radius of occurrence records
#' 
#' This function generates a sample of background coordinates that are within n kilometers of one occurrene record up to x background points per record.
#' @param coords A two-column matrix or data.frame of coordinates with the first column being longitude and the second being latitude.
#' @param clim A raster object to extract background data from.
#' @param radius A distance in km to define a radius around each occurrence to sample from.
#' @param n Number of background points to generate per occurrence point. Note that background points will be thinned so that only one occurs per grid cell, so the total number returned will be less than n * length(coords[,1]).
#' @export
#' @examples \dontrun{
#' data(distr);
#' bg.ext <- rad_bg(distr[,4:3], climondbioclim, radius=100, n = 50)
#' }

rad_bg <- function(coords, clim, radius, n){
  extr = coords;
  ##Premise: IF a species PDF - P(sp) - is approximately equal to the background PDF
  #P(bg) then P(condi) = P(sp)/P(bg) should be a uniform distribution and
  #not contribute to the overall estimation of climate.
  
  bg.mat <- matrix(ncol = 2, nrow = n * nrow(extr));
  for(i in 1:length(extr[,1])){
   # print(i)
    for(zz in 1:n){
     # print(paste('zz', zz))
      dir = sample(1:360, 1)
      dist = sample(1:radius, 1);
      bg.mat[(i*zz),1:2] = findcoord(extr[i,1], extr[i,2], dist, dir)
      
    }
  }
  bg.mat = cbind(rep(1111, length(bg.mat[,1])), rep('bg', length(bg.mat[,1])), bg.mat[,2], bg.mat[,1])
#  bg.mat = stats::na.omit(bg.mat)
  colnames(bg.mat) = c('ind_id', 'tax', 'lat', 'lon')
  bg.mat <- data.frame(bg.mat)
  bg.mat$lon = as.numeric(as.character(bg.mat$lon))
  bg.mat$lat = as.numeric(as.character(bg.mat$lat))
# bg.mat = unique(bg.mat);
  bg.ext <- raster::extract(clim, bg.mat[,4:3], cellnumbers=T)
  bg.ext <- cbind(bg.mat, bg.ext);
  bg.ext <- stats::na.omit(bg.ext);
  
  bg.ext = bg.ext[!duplicated(bg.ext[,"cells"]),]
  
  return(bg.ext);
  
}

#' A wrapper for vegdistmod::densform where a multi-taxon extraction object can be passed to densform one taxon at a time.
#' 
#' This function takes extracted climate data (from an object generated by the vegdistmod::extraction() function) and generates probability density functions for each taxon/variable pair using both a Gaussian (normal) approximation and a Gaussian Kernel Density estimator.
#' @param ex An object derived from the extraction() function.
#' @param clim A raster object (see raster::raster() and raster::stack() documentation for reading raster files into R).
#' @param bw A bandwidth compatible with stats::density(). Options include "nrd", "nrd0", "ucv", "bcv", etc.. Default (and recommended) value is "nrd0".
#' @param kern Type of Kernel to smooth with. Recommend 'gaussian', 'optcosine', or 'epanechnikov'. See: stats::density for options.
#' @param n Number of equally spaced points at which the probability density is to be estimated. Defaults to 1024. A lower number increases speed but decreases resolution in the function. A higher number increases resolution at the cost of speed. Recommended values: 512, 1024, 2048, ....
#' @param clip A character string of value "range" or "95conf" or "99conf". Should the probability functions be clipped to either the empirical range or the 95 or 99 percent confidence interval? 
#' @param manip Character string of 'reg' for straight likelihood, 'condi' for conditional likelihood statement.
#' @param parallel TRUE or FALSE. Make use of multicore architecture.
#' @param nclus Number of cores to allocate to this function
#' @param bg.n If there is not a background matrix, how many background points PER OCCURRENCE record should be sampled. Default is 1000.
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, 
#'  schema='flat', factor = 4, rm.outlier=FALSE);
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, 
#'  manip = 'condi', bg.n = 200, bw = 'nrd0', n = 1024);
#' multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' }

dens_obj <- function(ex, clim, manip = 'condi', bw = "nrd0", kern='optcosine',
                     clip = 0, n = 1024, parallel = FALSE, 
                     nclus = 4, bg.n = 200) {
	rawbioclim = clim;
	ex <- data.frame(ex);
	condi = FALSE;
	bayes = FALSE;
	head = which(colnames(ex) == 'cells');
  from = raster::minValue(clim);
  to = raster::maxValue(clim);
	####
	
	
	
	if(manip == 'condi') {
	  condi = TRUE; #print("Conditional Likelihood")
	}
	dens.list <- list();
	nlist <- vector();
	site.ex <- "NOSITE";

	site.coord = 0;
	if(ex[1,2] == "SITECOORD"){	
		site.coord <- ex[1,];
		ex <- subset(ex, ex$tax != "SITECOORD");
		site.ex <- ex[1,];
		
	};
	tax.list <- unique(ex$tax);
	tax.list <- stats::na.omit(tax.list);
	if(parallel == TRUE){
	  
	  cl <- parallel::makeCluster(nclus, type = "SOCK")
	  doSNOW::registerDoSNOW(cl);
	  
	  dens.list <-
	    foreach::foreach(i = 1:length(tax.list),
	           
	            .packages = 'vegdistmod') %dopar% {
	           #   source('~/Desktop/cracle_testing/vegdistmod/R/search_fun.R')
	           #   source('~/Desktop/cracle_testing/vegdistmod/R/cracle_build.R')
	              s.ex <- subset(ex, ex$tax == tax.list[[i]]);
	              
	              s.ex <- stats::na.omit(s.ex);
	              
	              dlist <- (densform(s.ex, rawbioclim, name = tax.list[[i]], 
	                                 manip = manip, bw = bw, kern=kern, 
	                                 n=n, from = from, to = to,
	                                 clip = clip, bg.n = bg.n));
	              
	              len <- length(dlist);
	              if(len <= 1) {
	                dlist <- NULL;
	              };
	              return(dlist);
	              
	              
	            }
	  parallel::stopCluster(cl)
	  
	} else {
	  
	
	for(i in 1:length(tax.list)){	
	 # print(i);
	  	s.ex <- subset(ex, ex$tax == tax.list[[i]]);
		
	  	s.ex <- stats::na.omit(s.ex);
		
	  	nlist[[i]] <- length(s.ex[,1])

  		dens.list[[i]] <- (densform(s.ex, rawbioclim, name = tax.list[[i]], 
  		                            manip = manip, bw = bw, kern = kern, 
  		                            n=n, from = from, to = to, 
  		                            clip = clip, bg.n = bg.n));

	 	  len <- length(dens.list[[i]]);
		  if(len <= 1) {
		  	dens.list[[i]] <- NULL;
	  	};
  	};
	};
	return(dens.list);
}

#dens_obj <- compiler::cmpfun(dens_obj);

#' P(A | B) = P(A) + P(B)
#' 
#' Using an object from the vegdistmod::dens_obj() function. Create a single density object (i.e., like that produced by vegdistmod::densform()) where the probability curves correspond to the probability density function of any one taxon/species from the original set occurring. This is not actually used in the implementation of finding the maximum joint likelihood in a CRACLE analysis, but is a good companion to the vegdistmod::and_fun() function.
#' @param dens.oblist An object derived from the vegdistmod::dens_ob() function.
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");

#' #OR:
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' or <- or_fun(dens.list.raw);
#' addplot(or, names(climondbioclim[[1]]), col ='black');
#' }

or_fun <- function(dens.oblist){
	varlist <- names(dens.oblist[[1]]);
	varlist <- (varlist[1:((length(varlist)-1)/6)]);
	varlist <- sub(".kde", "", varlist);

	field <- list();
	gfield <- list();
	xfield <- list();
	meanadjust <- list();
	variances <- list();
	name = "ADDITION";
	for (n in 1:length(varlist)){
		var = varlist[n];
		varx <- paste(var, "x", sep = ".");
		vargauss <- paste(var, "gauss", sep = ".");
		varkde <- paste(var, "kde", sep = ".");

		varmean <- paste(var, "mean", sep = ".");
		varsd <- paste(var, "sd", sep = ".");
		meanlist <- list();
		sdlist <- list();
		dens.obcurr <- dens.oblist[[1]];
		to <- max(dens.obcurr[[varx]]);
		from <- min(dens.obcurr[[varx]]);
		num = length(dens.obcurr[[varx]]);
		by = (to - from)/num;
		meanlist[[1]] <- as.numeric(dens.obcurr[[varmean]]);
		sdlist[[1]] <- as.numeric(dens.obcurr[[varsd]])^2;
		prod <- as.numeric(dens.obcurr[[varkde]]);
		prod.gauss <- as.numeric(dens.obcurr[[vargauss]]);

		for(i in 2:length(dens.oblist)){
			dens.obnow <- dens.oblist[[i]];
			prod <- prod + (as.numeric(dens.obnow[[varkde]]));
			prod.gauss <- prod.gauss + (as.numeric(dens.obnow[[vargauss]]));
			prod.area <- sum(prod)*by;
			prod <- prod / prod.area;
			prod.gauss.area <- sum(prod.gauss)*by;
			prod.gauss <- prod.gauss / prod.gauss.area;
			meanlist[[i]] <- as.numeric(dens.obnow[[varmean]]);
			sdlist[[i]] <- as.numeric(dens.obnow[[varsd]])^2;
		};
		prod.area <- sum(prod)*by;
		prod <- prod / prod.area;
		prod.gauss.area <- sum(prod.gauss)*by;
		prod.gauss <- prod.gauss / prod.gauss.area;
		field[[n]] <- prod;
		gfield[[n]] <- prod.gauss;
		xfield[[n]] <- dens.obcurr[[varx]];
		meanadjust[[n]] <- as.numeric(meanlist)/as.numeric(sdlist);
		variances[[n]] <- 1/as.numeric(sdlist);
	};
	meansum <- lapply(meanadjust, sum);
	varisum <- lapply(variances, sum);
	wmeans <- mapply("/", meansum, varisum);
	wsd <- mapply("/", 1, varisum);
	wsd <- lapply(wsd, sqrt);
	field <- data.frame(field);
	gfield <- data.frame(gfield);
	xfield <- data.frame(xfield);
	colnames(field) <- (paste(varlist, "kde", sep = "."));
	colnames(gfield) <- paste(varlist, "gauss", sep = ".");
	names(wmeans) <- paste(varlist, "mean", sep = ".");
	names(wsd) <- paste(varlist, "sd", sep = ".");
	colnames(xfield) <- (paste(varlist, "x", sep = "."));
	name = data.frame(name);
	colnames(name) <- "name";
	fin <- c(field, xfield, gfield, wmeans, wsd, name);
	fin <- .makeaucone(fin);
	return(fin);
};

#and_fun = compiler::cmpfun(and_fun);

#' P(A | B) = P(A) * P(B)
#' 
#' Using an object from the vegdistmod::dens_obj() function. Create a single density object (i.e., like that produced by vegdistmod::densform()) where the probability curves correspond to the probability density function of ALL taxa/species from the original set occurring. 
#' @param dens.oblist An object derived from the vegdistmod::dens_ob() function.
#' @param w Weight importance of probability functions
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' or <- or_fun(dens.list.raw);
#' addplot(or, names(climondbioclim[[1]]), col ='black');
#' and <- and_fun(dens.list.raw);
#' addplot(and, names(climondbioclim[[1]]), col ='black');
#' }

and_fun <- function(dens.oblist, w = FALSE){
	dens.oblist <- .scramble(dens.oblist);
	varlist <- names(dens.oblist[[1]]); #print(varlist)
	varlist <- (varlist[1:((length(varlist)-1)/6)]) ;
	varlist <- sub(".kde", "", varlist);

	field <- list();
	gfield <- list();
	xfield <- list();
	meanadjust <- list();
	variances <- list();
	name = "PRODUCT";
	for (n in 1:length(varlist)){ #print(varlist[n])
		var = varlist[n]; 
		varx <- paste(var, "x", sep = ".");
		varkde <- paste(var, "kde", sep = ".");

		vargauss <- paste(var, "gauss", sep = ".");
		varmean <- paste(var, "mean", sep = ".");
		varsd <- paste(var, "sd", sep = ".");
		varw <- paste(var, "w", sep = '.');
		meanlist <- list();
		sdlist <- list();

		dens.obcurr <- dens.oblist[[1]];
		
		if(w == TRUE){ 
		  we = dens.obcurr[[varw]]; #print(we);
		} else {
		  we <- 1;
		}
		
		to <- max(stats::na.omit(dens.obcurr[[varx]]));
		from <- min(stats::na.omit(dens.obcurr[[varx]]));
		num = length((dens.obcurr[[varx]]));
		by = (to - from)/num;
		meanlist[[1]] <- as.numeric(dens.obcurr[[varmean]]);
		sdlist[[1]] <- as.numeric(dens.obcurr[[varsd]])^2;
		prod <- as.numeric(dens.obcurr[[varkde]]) ^ we;
		prod <- prod*by;
		prod.gauss <- as.numeric(dens.obcurr[[vargauss]])*by;
		for(i in 2:length(dens.oblist)){
		  dens.obnow <- dens.oblist[[i]];
		  if(w == TRUE){
  		  we = dens.obnow[[varw]]; #print(we);
		  } else {
		    we = 1;
		  }
		  if(sum(stats::na.omit(dens.obnow[[varkde]]*by)) == 0) {next;}
  		prod <- prod * ((as.numeric(dens.obnow[[varkde]])*by) ^ we);
			prod.area <- sum(stats::na.omit(prod))*by;
			prod <- prod / prod.area;
			prod.gauss <- prod.gauss * ((as.numeric(dens.obnow[[vargauss]])*by)^we);
			prod.gauss.area <- sum(stats::na.omit(prod.gauss))*by;
			prod.gauss <- prod.gauss / prod.gauss.area;
			meanlist[[i]] <- as.numeric(dens.obnow[[varmean]]);
			sdlist[[i]] <- as.numeric(dens.obnow[[varsd]])^2;
		};
		prod.area <- sum(stats::na.omit(prod))*by;
		prod <- prod / prod.area;
		prod.gauss.area <- sum(stats::na.omit(prod.gauss))*by;
		prod.gauss <- prod.gauss / prod.gauss.area;
		field[[n]] <- prod;
		gfield[[n]] <- prod.gauss;
		xfield[[n]] <- dens.obcurr[[varx]];
		meanadjust[[n]] <- as.numeric(unlist(meanlist))/as.numeric(unlist(sdlist));
		variances[[n]] <- 1/as.numeric(unlist(sdlist));
	};
  meansum <- lapply(meanadjust, sum);
	varisum <- lapply(variances, sum);
	wmeans <- mapply("/", meansum, varisum);
	wsd <- mapply("/", 1, varisum);
	wsd <- lapply(wsd, sqrt);
	field <- data.frame(field);
	gfield <- data.frame(gfield);
	xfield <- data.frame(xfield);
	colnames(field) <- paste(varlist, "kde", sep = ".");
	colnames(gfield) <- paste(varlist, "gauss", sep = ".");
	names(wmeans) <- paste(varlist, "mean", sep = ".");
	names(wsd) <- paste(varlist, "sd", sep = ".");
	colnames(xfield) <- (paste(varlist, "x", sep = "."));
	name = data.frame(name);
	colnames(name) <- "name";
	fin <- c(field, xfield, gfield, wmeans, wsd, name);
	fin <- .makeaucone(fin);
	return(fin);
};

#or_fun = compiler::cmpfun(or_fun);

#get_optim() takes an object output from the densform function or and_fun or or_fun and finds optimal values for each PDF
#' Find PDF optim(a)um
#' 
#' Using an object from the vegdistmod::dens_obj() function. Create a single density object (i.e., like that produced by vegdistmod::densform()) where the probability curves correspond to the probability density function of ALL taxa/species from the original set occurring. 
#' @param dens.ob An object derived from the vegdistmod::dens_ob() function.
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' and <- and_fun(dens.list.raw);
#' addplot(and, names(climondbioclim[[1]]), col ='black');
#' optim.and <- get_optim(and);
#' abline(v=optim.and$means[paste(names(climondbioclim[[1]]), 'mean', sep = ".") ])
#' }

get_optim <- function(dens.ob){
	dens.ob1 <- dens.ob;
	varlist <- names(dens.ob1);
	varlist <- (varlist[1:((length(varlist)-1)/5)]);
	varlist <- sub(".kde", "", varlist);
	conintkde <- list();
	conintgauss <- list();
	dirconint <- list();
	origk <- list();
	origg <- list();
	means <- list();
	sds <- list();
	for (j in 1:length(varlist)){
	  
		var = varlist[[j]];
		varx <- paste(var, "x", sep = ".");
		#print(varx)
		vargauss <- paste(var, "gauss", sep = ".");
		varkde <- paste(var, "kde", sep = ".");

		varmean <- paste(var, "mean", sep = ".");
		varsd <- paste(var, "sd", sep = ".");
		cumulkde <- list();
		cumulgauss <- list();
		cikde <- list(0,0);
		cigauss <- list(0,0);
		runkde <- 0;
		rungauss <- 0;
		to <- max(stats::na.omit(dens.ob1[[varx]]));
		from <- min(stats::na.omit(dens.ob1[[varx]]));
		num = length(stats::na.omit(dens.ob1[[varx]]));
		by = (to - from)/num;
		for (i in 1:length(dens.ob1[[varkde]])){
			runkde = runkde + (dens.ob1[[varkde]][i]*by);
			cumulkde[[i]] <- runkde;
			if(is.na(cumulkde[[i]])){cumulkde[[i]] = 0;}
			if(i==1){ 
				if(cumulkde[[i]] >= 0.025){
					cikde[[1]] <- dens.ob1[[varx]][i];
				};
				if(cumulkde[[i]] >= 0.975){
					cikde[[1]] <- dens.ob1[[varx]][i];
					cikde[[2]] <- dens.ob1[[varx]][i];
				};
			} else {
				if(cumulkde[[i-1]] < 0.025 && cumulkde[[i]] >= 0.025){
					cikde[[1]] <- dens.ob1[[varx]][i];
				};
			if(cumulkde[[i-1]] < 0.975 && cumulkde[[i]] >= 0.975){
				cikde[[2]] <- dens.ob1[[varx]][i];
			};
		};
		rungauss = rungauss + (dens.ob1[[vargauss]][i]*by);
		cumulgauss[[i]] <- rungauss;
		if(is.na(cumulgauss[[i]])){cumulgauss[[i]] = 0;}
		
		if(i==1){ 
			if(cumulgauss[[i]] >= 0.025){
				cigauss[[1]] <- dens.ob1[[varx]][i];
			};
			if(cumulgauss[[i]] >= 0.975){
				cigauss[[1]] <- dens.ob1[[varx]][i];
				cigauss[[2]] <- dens.ob1[[varx]][i];
			};
		} else {
			if(cumulgauss[[i-1]] < 0.025 && cumulgauss[[i]] >= 0.025){
				cigauss[[1]] <- dens.ob1[[varx]][i];
			};
			if(cumulgauss[[i-1]] < 0.975 && cumulgauss[[i]] >= 0.975){
				cigauss[[2]] <- dens.ob1[[varx]][i];
			};
		};
	};
	logkde <- ifelse(dens.ob1[[varkde]]>0, log(dens.ob1[[varkde]]*by), -Inf);
#	print(varkde);
	origkde <- subset(dens.ob1[[varx]], logkde >= max(stats::na.omit(logkde))*1.01);
	origk[[j]] <- c(min(stats::na.omit(origkde)), max(stats::na.omit(origkde)));
	loggauss <- ifelse(dens.ob1[[vargauss]]>0, log(dens.ob1[[vargauss]]*by), -Inf);
#	loggauss <- log(dens.ob1[[vargauss]]*by)
	origgauss <- subset(dens.ob1[[varx]], loggauss >= max(stats::na.omit(loggauss))*1.01); 
	origg[[j]] <- c(min(stats::na.omit(origgauss)), max(stats::na.omit(origgauss)));
	conintkde[[j]] <- c(cikde[[1]], cikde[[2]]);
	conintgauss[[j]] <- c(cigauss[[1]], cigauss[[2]]);
	dirconint[[j]] <- c((dens.ob1[[varmean]] - 1.96*dens.ob1[[varsd]]), (dens.ob1[[varmean]]+1.96*dens.ob1[[varsd]]));
	means[[j]] <- dens.ob1[[varmean]];
	sds[[j]] <- dens.ob1[[varsd]];
	};
	conintkde <- data.frame(conintkde);
	conintgauss <- data.frame(conintgauss);
	origk <- data.frame(origk);
	origg <- data.frame(origg);
	dirconint <- data.frame(dirconint);
	means <- data.frame(means);
	sds <- data.frame(sds);
	colnames(conintkde) <- paste(varlist, "cikde", sep = ".");
	colnames(conintgauss) <- paste(varlist, "cigauss", sep = ".");
	colnames(origk) <- paste(varlist, "origkde", sep = ".");
	colnames(origg) <- paste(varlist, "origgauss", sep = ".");
	colnames(dirconint) <- paste(varlist, "cidir", sep = ".");
	colnames(means) <- paste(varlist, "mean", sep = ".");
	colnames(sds) <- paste(varlist, "sd", sep = ".");
	ret <- list(conintkde, conintgauss, origk, origg, dirconint, means, sds);
	names(ret) <- c("conintkde", "conintgauss", "origk", "origg", "dirconint", "means", "sds");
	return(ret);
};

#get_optim = compiler::cmpfun(get_optim);


#makes area under any PDF curve equal 1 (good for standardizing curves to be compared). HIDDEN!
.makeaucone <- function(dens.ob1, var){ var <- names(dens.ob1);
	var <- (var[1:((length(var)-1)/6 )]);
	var <- sub('.kde', '', var)
	for(i in 1:length(var)){
		varnow <- var[[i]];
		varx <- paste(var[[i]], "x", sep = ".");
		gauss <- paste(var[[i]], "gauss", sep = ".");
		kde <- paste(var[[i]], "kde", sep = ".");

		to <- max(subset(dens.ob1[[varx]], !is.na(dens.ob1[[kde]])));
		from <- min(subset(dens.ob1[[varx]], !is.na(dens.ob1[[kde]])));
		num <- length(subset(dens.ob1[[varx]], !is.na(dens.ob1[[kde]])));
		by = (to - from)/num;
		do <- sum(stats::na.omit(dens.ob1[[kde]]))*by;
		if(do == '0'){
		  return(0)
		}
		do.gauss <- sum(stats::na.omit(dens.ob1[[gauss]]))*by;
		dens.ob1[[kde]] <- dens.ob1[[kde]]/do;
		dens.ob1[[gauss]] <- dens.ob1[[gauss]]/do.gauss;
	};
	return(dens.ob1);
};

#scramble() reorders pdfs. No real reason to do this as order does not matter for the operations being done here.
.scramble <- function(x, k=3) {
	x.s <- seq_along(x);
	y.s <- sample(x.s);
	idx <- unlist(split(y.s, (match(y.s, x.s)-1) %/% k), use.names = FALSE);
	x[idx];
};


#' Plot PDF curve of given type and variable for single density object.
#' 
#' Using an object from the densform() function. Plot a single PDF curve in a new plot window.
#' @param dens.ob An object derived from the vegdistmod::dens_ob() function.
#' @param var A character string that matches one of the layer names in the source raster object.
#' @param col A color declaration. Default is a random color.
#' @param type A character string of value either ".kde" for a Kernel Density Estimator curve, or ".gauss" for a Gaussian (normal) curve. All other values will result in errors.
#' @param w TRUE or FALSE to show weighted probability functions
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr); data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' extr.sub = subset(extr.raw, extr.raw$tax == extr.raw[5,'tax']);
#' dens.sub = densform(extr.sub, clim = climondbioclim, bw = 'nrd0', n = 512);
#' densplot(dens.sub, names(climondbioclim[[1]]));
#' }

densplot <- function(dens.ob, var, col = sample(grDevices::colours()), type = ".kde", w=FALSE) {
	
	varx <- paste(var, "x", sep = ".");
	varw <- paste(var, "w", sep = ".");
	
	graphics::par(mar= c(5,4,4,4) + 0.3);
	tempvarlist <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "MAT", "MaximumT", "MinimumT");
	if(var %in% tempvarlist){by = 10}else{by = 1};
	var <- paste(var, type, sep = "");
	if(w ==TRUE){
	  to <- max(stats::na.omit(dens.ob[[varx]]));
	  from <- min(stats::na.omit(dens.ob[[varx]]));
	  num = length(dens.ob[[varx]]);
	  lby = (to - from)/num;
	  dens.ob[[var]] = dens.ob[[var]]^dens.ob[[varw]];
	  den.area <- sum(stats::na.omit(dens.ob[[var]]))*lby;
	  dens.ob[[var]] = dens.ob[[var]]/den.area
	}
	graphics::plot(dens.ob[[varx]]/by, dens.ob[[var]], xlab = "", ylab = "", ylim = c(0, 3.5*max(stats::na.omit(dens.ob[[var]]))), type = "l", lwd = 3, col = col, frame.plot=F, axes = F);
	graphics::axis(side = 2, at = pretty(c(0, 2.5*max(stats::na.omit(dens.ob[[var]])))));
	graphics::axis(side = 1, at = pretty(range(stats::na.omit(dens.ob[[varx]])/by)));
	graphics::mtext(var, side = 1, line =3);
	graphics::mtext("Probability Density Estimation", side = 2, line = 3);
};



#' Plot PDF curves of given type and variable for density list object
#' 
#' Using an object from the vegdistmod::dens_obj() function. Plot a series PDF curves in a new plot window.
#' @param dens.oblist An object derived from the vegdistmod::dens_obj() function.
#' @param var A character string that matches one of the layer names in the source raster object.
#' @param col A color vector of the same length as the dens.oblist object. Default is 'grDevices::heat.colors(length(dens.oblist))'.
#' @param type A character string of value either ".kde" for a Kernel Density Estimator curve, or ".gauss" for a Gaussian (normal) curve. All other values will result in errors.
#' @param l.pos  Legend position. Recommend 'topleft' or 'topright'. Default is 'topleft'.
#' @param l.cex  cex setting for legend. Default is 0.8.
#' @param w TRUE or FALSE to show weighted probability functions
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr); data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' }

multiplot <- function(dens.oblist, var, col = grDevices::heat.colors(length(dens.oblist)), type = ".kde", l.pos = 'topleft', l.cex = 0.8, w = FALSE){ 
	arr.dens.ob = dens.oblist;
	varx <- paste(var, "x", sep = ".");
	vart = paste(var, type, sep = '');

	current <- arr.dens.ob[[1]];
	densplot(current, var, col[1], type = type, w=w);
	max.x.hold = list(max(current[[varx]]));
	max.y.hold = list(max(current[[vart]]));
	names.hold = as.character(current[["name"]]);
	for(i in 2:length(arr.dens.ob)){
		current <- arr.dens.ob[[i]];
		addplot(current, var, col[i], type = type, w=w);
#		max.x.hold = c(max.x.hold, max(current[[varx]]));
#		max.y.hold = c(max.y.hold, max(current[[vart]]));
		names.hold = c(names.hold, as.character(current[["name"]]));
	};
#	max.x <- mean(as.numeric(as.character(max.x.hold)));
#	max.y <- mean(as.numeric(as.character(max.y.hold)));
	graphics::legend(l.pos, legend = as.character(names.hold), lty=1, lwd=2, cex=l.cex, col = col, box.col=NA);
};

#' Adds a single PDF plot to already open plot
#' 
#' Using an object from the vegdistmod::densform() or vegdistmod::and_fun() or vegdistmod::or_fun() functions, add a PDF over the existing plot. Useful for visualizing joint likelihood curves.
#' @param dens.ob An object derived from the vegdistmod::densform() or vegdistmod::and_fun() or vegdistmod::or_fun() functions.
#' @param var A character string that matches one of the layer names in the source raster object.
#' @param col A color declaration. Default is a random color.
#' @param type A character string of value either ".kde" for a Kernel Density Estimator curve, or ".gauss" for a Gaussian (normal) curve. All other values will result in errors.
#' @param w TRUE or FALSE to show weighted probability functions
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr); data(climondbioclim);
#' #extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' #dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' #multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' #or <- or_fun(dens.list.raw);
#' #addplot(or, names(climondbioclim[[1]]), col ='black');
#' #and <- and_fun(dens.list.raw);
#' #addplot(and, names(climondbioclim[[1]]), col ='black');
#' }

addplot <- function(dens.ob, var, col = sample(grDevices::colours()), type = ".kde", w=FALSE) {
	varx <- paste(var, "x", sep = ".");
	varw <- paste(var, "w", sep = ".");
	
	tempvarlist <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "MAT", "MaximumT", "MinimumT");
	if(var %in% tempvarlist){by = 10}else{by = 1};
	var <- paste(var, type, sep = "");
	if(w ==TRUE){
	  to <- max(dens.ob[[varx]]);
	  from <- min(dens.ob[[varx]]);
	  num = length(dens.ob[[varx]]);
	  lby = (to - from)/num;
	  dens.ob[[var]] = dens.ob[[var]]^dens.ob[[varw]];
	  den.area <- sum(dens.ob[[var]])*lby;
	  dens.ob[[var]] = dens.ob[[var]]/den.area
	}
	graphics::points(dens.ob[[varx]]/by, dens.ob[[var]], type = "l", lwd = 3, col = col);
};

#' Write results to file
#' 
#' Write results from a results object from either vegdistmod::get_optim() to a tab delimited table in a text file. This is done one 'method' at a time. Supported methods are "conintkde", "conintgauss", "origg", "origk", "dirconint" for the "get_optim()" type objects.
#' @param optima A results object from vegdistmod::get_optim()
#' @param siteval A vector of site values for envirnomental parameters (See example for generating).
#' @param clim The original raster object used in the vegdistmod::extraction() step.
#' @param method The result summary method. Supported methods are "conintkde", "conintgauss", "origg", "origk", "dirconint" for the "get_optim()" type objects.
#' @param filename Path to desired file.
#' @param append Should this be appended to an existing file. Default is FALSE.
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr); data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' and <- and_fun(dens.list.raw);
#' optima <- get_optim(and);
#' #Compare to optima.
#' write_results(optima, clim = climondbioclim, method = 'conintkde', file = 'raw')
#' }

write_results <- function(optima, siteval = '', clim, method, filename, append = F){
	# type and file name writes output either
	# header including: site-coord, ntax, taxlist (if type = "hdr")
	#or a table of results if type = 'res'
	# site values, and results and comparison to site value.
	
			
		nvars <- raster::nlayers(clim);
		comp = matrix(nrow = nvars, ncol = 6);
		for(i in 1:nvars){		
			comp[i,2] = as.matrix(siteval[i+6]);
			comp[i,3] = optima[[method]][[i]][1];
			comp[i,4] = optima[[method]][[i]][2];
			comp[i,5] = as.numeric(comp[i,3])-as.numeric(comp[i,2]);
			comp[i,6] = as.numeric(comp[i,4])-as.numeric(comp[i,2]);
		};
		comp[,1] = names(clim);

		colnames(comp) = c("climate_variable", "site_value", paste(method, "min", sep = "_"), paste(method, "max", sep = "_"), "min_resid", "max_resid");
		file = paste(filename, method, 'tab', sep =".")
		utils::write.table(comp, file=file, append=append, sep="\t", quote=F, row.names=F, col.names=T);

}



#' Write a companion results file providing site details and final taxon list.
#' 
#' Companion file will include the environmental data extracted for the site locality, site locality information (lat/lon), and the list of taxa used in this analysis as derived from an extraction object. This list may differ from the original list requested if some taxa were excluded due to too few records left once filtered.
#' @param site.coord A vector of the site coordinates and any other desired site metadata.
#' @param tax.list A list of unique taxon names. Usually obtained by running 'tax.list <- unique(extraction_object$tax)' on your 'extraction_object' (See example).
#' @param filename The path to desired file.
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr); data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' site.coord <- extr.raw[1,1:5]; 
#' #Assuming that the site details are given in the first row of the original datafile. 
#' tax.list <- unique(subset(extr.raw, extr.raw$tax != 'SITECOORD')$tax); 
#' 
#' #Assumes that the site details are given in the original 
#' #data file and that the entry in the 'tax' column for the site is "SITECOORD". 
#' #Substitute whatever you called your site if you did it this way. OR:
#' #tax.list <- unique(extr.raw$tax);
#' write_headfile(site.coord=site.coord, tax.list=tax.list, filename = "headfile.try.hdr")
#' }

write_headfile <- function(site.coord = '', tax.list = '', filename){
			
			file = paste(filename, 'hdr', 'txt', sep = ".")
			utils::write.table("sitecoord", file=file, append=F, sep="\t", quote=F, col.names = F,  row.names=F);
			utils::write.table(site.coord, file=file, append=T, sep="\t", quote=F, col.names = F, row.names=F);
			utils::write.table("TaxNum", file=file, append =T, sep = "\t", quote=F, col.names = F, row.names=F);
			utils::write.table(length(tax.list), file=file, append =T, sep = "\t", quote=F, col.names = F,  row.names=F);
			utils::write.table("TaxList", file=file, append =T, sep = "\t", quote=F, col.names = F, row.names=F);
			utils::write.table(tax.list, file=file, append =T, sep = "\t", quote=F, col.names = F, row.names=F);
	
}



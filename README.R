###vegdistmod README File####
##### Robert S. Harbert #####
########## 2016 #############
##Installation note: 
#If you do not have vegdistmod installed uncomment 
# the commands below and run them ==>

#install.packages('devtools');
#library(devtools);
#install_git('git://github.com/rsh249/vegdistmod.git')

###vegdistmod should now be installed...

#What is written below is code that makes use of vegdistmod
#functions to perform some of the analyses available in this
#library. This code has been tested and works on R3.3.1 and R3.2.3
# on Ubuntu and CentOS Linux distros as well as MacOSX El Capitan 10.11.6
#and likely works with other versions that have up-to-date dependencies.

### The vegdistmod package includes a handful of test datasets.
library(vegdistmod)


#Distribution for the genus Abies, a northern hemisphere conifer tree.
data(abies)
data(climondbioclim) #CliMond data for Eastern North America
head(abies) #This is the standard occurrence data format for vegdistmod
            #ALL occurrence records have an individual identifier, 
            #in matrices with simulated data the ind_id = 0000.
#View raw GBIF data
plot_clim(abies, climondbioclim[[1]])
     # sample to include one record per 8x8 grid cell area (to reduce over-sampling bias)

head(ext.abies)
#View filtered distribution:
plot_clim(ext.abies, climondbioclim[[1]])

#generate background points
bg <- vegdistmod:::.get_bg(climondbioclim[[1]])
bg.ext <- extraction(bg, climondbioclim, schema='raw');

#search for probable simulated occurrences using the standard algorithm
find.abies <- findlocal(ext.abies, climondbioclim, bg = bg.ext, type = '.kde', 
                        maxiter = 20, searchrep = 1, manip = 'reg', 
                        alpha = 0.05, factor = 4)
      #This runs up to 20 iterations of the search algorithm but may end early if 
      # 'better' model/sample pairs are not being found (10 iterations without change).
      #This also uses Kernel Density Estimation likelihood functions. Try the ".gauss" option for the type object to change this.

#Look at the simulated distribution
plot_clim(find.abies[[1]], climondbioclim[[1]])
    #We have a problem here. This is too broad for Abies. Specifically there are too many records in the Midwestern US.

#Maybe the problem is local populations (species in this case potentially)
#have different affinities for climate. Try the same search with 
# 10 iterations of geographic subsetting before searching.
#This will be run in parallel using 4 cores by default. 
# Set nclus = some other number to make use of more or less than 4 cores.

geofind.abies <- geo_findlocal(ext.abies, climondbioclim, bg = bg.ext, type = '.kde', 
                           maxiter = 10, searchrep = 1, manip = 'reg', 
                           alpha = 0.05, factor = 4, divisions = 4, parallel = TRUE, nclus = 4)
plot_clim(geofind.abies, climondbioclim[[1]])
#compare back to where we started:
plot_clim(abies, climondbioclim[[1]])

##IF you have time and the hardware to do it, try a more intense search 
#using geo_findlocal:
#NOTE if you have access to more than 4 cores now would be a good time to use them. 
#  set nclus equal to the number of cores you can use.
#
#geofind.abies <- geo_findlocal(ext.abies, climondbioclim, bg = bg.ext, type = #	'.kde', maxiter = 50, searchrep = 1,
#	 manip = 'condi', alpha = 0.05, factor = 4, divisions = 80, parallel = TRUE, #	nclus = 80)

#plot_clim(geofind.abies, climondbioclim[[1]])






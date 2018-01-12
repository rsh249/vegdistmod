#include<iostream>
#include<cmath>
#include <Rcpp.h>
using namespace std;
// [[Rcpp::export]]


float distance(double lon1, double lat1, double lon2, double lat2)
{
    //Haversine distance method
    float R = 6378.137;
    
    float toRad = 3.14159/180;
    lon1 = lon1 * toRad;
    lon2 = lon2 * toRad;
    lat1 = lat1 * toRad;
    lat2 = lat2 * toRad;
    float dlon = lon2 - lon1;
    float dlat = lat2 - lat1;
    
    double a = pow(sin(dlat / 2), 2) + (cos(lat1) * cos(lat2) * pow(sin(dlon / 2),2));
    
    double d = 2 * atan2(sqrt(a), sqrt(1 - a)) * R;
    
    return d ;
}

function X = infodistSimplex(a,b)
X = round(4 * asin( sqrt( sum( (sqrt(a) - sqrt(b)).^2 ) ) / 2 ) ,10);
return
    


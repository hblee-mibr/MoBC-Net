#' Distance Function
#' 
#' Calculate the distance through the dist matrix
#' 
#' @param distm input
#' @param x input (gene cluster 1)
#' @param  y input (gene cluster 2)
#'
#' get.shortest.dist
#' mean of shortest path lengths of all pariwised community genes 
#' @export 
#' 
get.shortest.dist = function(distm, x, y, avg.dist=NULL){
    sum(distm[x,y])/(length(x)*length(y))
}

#' get.kernel.dist
#' average of kernelized distance between communities
#' @export 
#' 
get.kernel.dist = function(distm, x, y, avg.dist=NULL){
    dist = -(sum(log(rowSums(exp(-(distm[x,y] +1))/length(y)))) + sum(log(colSums(exp(-(distm[x,y] +1))/length(x)))))/(length(x) + length(y))
    return(dist)
}

#' get.centre.dist
#' shortest path length between centers of communities
#' @export 
#' 
get.centre.dist = function(distm, x, y) {
    c1 = names(which.min(rowSums(distm[x,x])))
    c2 = names(which.min(rowSums(distm[y,y])))
    dist = distm[c1,c2]
    return(dist)
}

#'get.separation.dist
#' difference between shortest distance of communities and the average shortest distance of each community
#' @export 
#' 
get.separation.dist = function(distm, x, y, avg.dist=NULL){
    xy = distm[x,y] %>% mean
    x1 = distm[x,x] %>% mean
    y1 = distm[y,y] %>% mean
    xy-(x1+y1)/2 
}

#'get.closest.dist
#' average of sum of min of every community node based shortest path length
#' @export 
#' 
get.closest.dist = function(distm, xx1, yy1, avg.dist=NULL){
		x1 = distm[yy1,xx1] %>%apply(1,min) %>% sum 
		x2 = distm[yy1,xx1] %>%apply(2,min) %>% sum 
		(x1+x2)/(length(xx1)+length(yy1))   
}



#'get.s2b.dist
#' average of sum of min of every community node based shortest path length (< avg)
#' @export 
#' 
get.s2b.dist = function(distm, xx1, yy1, avg.dist){
        mm = distm[xx1,yy1]
        mm1 = mm[mm<avg.dist]
        sum(mm1)/length(mm1)
}



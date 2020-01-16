
# dim(x) = N*n, N has to be a power of 2
# i,j are the sorting dimensions of x

recursive_sort_3 <- function(l,x,i,j,n){

    l2 <- (l%%2)*j + ((l+1)%%2)*i 
    ind <- dim(x)[2]
    
    x <- x[order(x[,l2],decreasing=FALSE),]
    b <- floor(dim(x)[1]/2)
    x1 <- matrix(x[seq(1,b),],b,ind)
    x2 <- matrix(x[-seq(1,b),],dim(x)[1]-b,ind)
    
    if(l<(n-1)){
        x1 <- recursive_sort_3(l+1,x1,i,j,n)
        x2 <- recursive_sort_3(l+1,x2,i,j,n)
    }
    
    x <- rbind2(x1,x2)
}
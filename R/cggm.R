CGGM.mean <- function(data, h, g=NULL, silent=FALSE){
  
  result <- matrix(double(length(data)),nrow=nrow(data))

  if (is.null(g))
    g <- CGGM.autoscale(data,h,silent)
  
  .C("c_cggm",
     as.double(data),
     nrow(data),
     ncol(data),
     as.double(g),
     as.double(h),
     result,
     DUP=FALSE,
     PACKAGE="epsi")

  result
}

CGGM.lts <- function(data, h, g=NULL, trim=0, silent=FALSE){
  
  result <- matrix(double(length(data)),nrow=nrow(data))

  if (is.null(g))
    g <- CGGM.autoscale(data,h,silent)
  
  .C("c_cggm_lts",
     as.double(data),
     nrow(data),
     ncol(data),
     as.double(g),
     as.double(h),
     as.double(trim),
     result,
     DUP=FALSE,
     PACKAGE="epsi")

  result
}

CGGM.autoscale <- function(data, h, silent=FALSE){

  if (!silent)
    cat("calculating scale parameter")

  w <- ceiling(max(nrow(data),ncol(data))*h)
  
  q <- matrix(nrow=ncol(data)-2*w,ncol=ncol(data)-2*w)
  
  for (i in w:(ncol(data)-w)){
    if (!silent)
      cat(".")
    for (j in w:(nrow(data)-w))
      q[i-w,j-w] <- IQR(data[(i-w):(i+w),(j-w):(j+w)])
  }
  result <- median(q)
  if (!silent) {
    cat("\n")
    cat("scale parameter: ",result,"\n")
  }
  result
}

# PURPOSE: Convert original data to symbols
#---------------------------------------------------
# USAGE: symbz <- symbolize(x,y,z,m)
# where: x = data, dependent at t
#        y = data, dependent at t-1
#        z = data, independent
#        m = lag scope
#---------------------------------------------------
# RETURNS:
# symbz = a structure with sorted indices of
# ..$mhx, ..$mhy, ..$mhz, ..$mhxy, ..$mhxz and ...$mhxyz
#---------------------------------------------------
# COMMENTS: For joint use with TEntropy
#
# written by:
# Manuel Ruiz-Marin and Andres Romeu
# Universidad Politecnica de Cartagena and UMU

symbolize <- function(x,y,z,m) {
     
     T <- length(x)
     n <- T -m +1
     sx <- NULL
     sy <- NULL
     sz <- NULL
     sxy <- NULL
     sxz <- NULL
     sxyz <- NULL
     
     for (i in 1:n) {
          dx <- sort(x[i:(i+m-1)],decreasing=FALSE,index.return=TRUE)$ix       
          dy <- sort(y[i:(i+m-1)],decreasing=FALSE,index.return=TRUE)$ix
          dz <- sort(z[i:(i+m-1)],decreasing=FALSE,index.return=TRUE)$ix
          sx <- c(sx,dx)
          sy <- c(sy,dy)
          sz <- c(sz,dz)
          sxy <- c(sxy,dx,dy)    
          sxz <- c(sxz,dx,dz)
          sxyz <- c(sxyz,dx,dy,dz)
     }
     
     sx <- matrix(sx, nrow=m ,ncol =n )
     sy <- matrix(sy, nrow=m ,ncol =n )    
     sz <- matrix(sz, nrow=m ,ncol =n )
     sxy <- matrix(sxy, nrow=2*m, ncol=n )
     sxz <- matrix(sxz, nrow=2*m, ncol=n )
     sxyz <- matrix(sxyz, nrow=3*m, ncol=n )
     
     return( list(mhx = t(sx), mhy = t(sy),mhz = t(sz),mhxy = t(sxy), mhxz = t(sxz), mhxyz = t(sxyz)) )          
}

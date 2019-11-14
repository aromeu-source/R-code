# PURPOSE: Compute frequencies of rows in a matrix
#---------------------------------------------------
# USAGE: patf <- PattFreq(xmat)
# where  xmat = data matrix
#---------------------------------------------------
# RETURNS:
#   patf = a structure with entropy of
#              x, y, z, xy, xz and xyz
#---------------------------------------------------
# COMMENTS: For joint use with symbols(x,y,z)
#
# written by:
# Manuel Ruiz-Marin and Andres Romeu
# Universidad Politecnica de Cartagena

PattFreq <- function(xmat) {
     
     xmat <- as.data.frame(xmat)
     xacu <- as.data.frame(ftable(xmat))
     xacu <- xacu[xacu$Freq>0,]
     n <- nrow(xacu)
     k <- ncol(xacu)-1  
     
     return(list(patrones= data.matrix(xacu[,1:k]), frx=xacu$Freq/sum(xacu$Freq)))
}

# PURPOSE: Compute total entropy
#---------------------------------------------------
# USAGE: entr <- TEntropy(mh)
# where inputs are output from symbols(x,y,z) and
# x = data, dependent at t
# y = data, dependent at t-1
# z = data, independent       
#---------------------------------------------------
# RETURNS:
#    entr = a structure with entropy of
#              x, y, z, xy, xz and xyz
#---------------------------------------------------
# COMMENTS: For joint use with symbols(x,y,z)
#
# written by:
# Manuel Ruiz-Marin and Andres Romeu
# Universidad Politecnica de Cartagena

TEntropy <- function(mh) {
     
     nsx=PattFreq(mh$mhx)$frx
     nsy=PattFreq(mh$mhy)$frx
     nsz=PattFreq(mh$mhz)$frx
     nsxy=PattFreq(mh$mhxy)$frx
     nsxz=PattFreq(mh$mhxz)$frx
     nsxyz=PattFreq(mh$mhxyz)$frx
     
     entx=-sum(nsx*log(nsx))
     enty=-sum(nsy*log(nsy))
     entz=-sum(nsz*log(nsz))
     entxy=-sum(nsxy*log(nsxy))
     entxz=-sum(nsxz*log(nsxz))
     entxyz=-sum(nsxyz*log(nsxyz))
     
     entycondx=entxy-entx
     entycondxz=entxyz-entxz
     
     return( list(entx = entx, enty=enty, entz=entz, entxy=entxy, entxz=entxz, entxyz=entxyz, TE=entycondx-entycondxz) )
}

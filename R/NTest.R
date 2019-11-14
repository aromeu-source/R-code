# PURPOSE: NTE test
#---------------------------------------------------
# USAGE: nte <- NTest(data,noms,m,r,B)
# where pdataframe, first column Y second column X
# names = vector of names of variables to test
# m = dimension
# r = lags       
# B = number of bootstraps
#---------------------------------------------------
# RETURNS:
#    nte = a structure with 
#         testXcY, pvalXcY, testYcX, pvalYcX
#         ntest, p1ntest, p2ntest
#---------------------------------------------------
# written by:
# Andres Romeu
# Universidad de Murcia

NTest <- function(data,noms,m,r,B) {
     
     indice <- attributes(data)$index
     secciones <- unique(indice[[1]])
     npanels= length(secciones)
     XcYtot <- list(mhx= NULL,mhy= NULL,mhz= NULL,mhxy= NULL,mhxz= NULL,mhxyz= NULL)
     YcXtot <- list(mhx= NULL,mhy= NULL,mhz= NULL,mhxy= NULL,mhxz= NULL,mhxyz= NULL)
     
     for (j in 1:npanels) {
          
          # Selecciona panel por panel los datos -----------
          paneli <- subset(data,indice[[1]]==secciones[j])
          paneli <- subset(paneli,complete.cases(paneli))
          T=nrow(paneli);
          Y <- paneli[,noms[1]]
          X <- paneli[,noms[2]]
          
          # X causa Y --------------------------------------
          XcY <- symbolize(Y[1:T-r],Y[r+1:T],X[1:T-r],m)
          XcYtot$mhx= rbind(XcYtot$mhx,XcY$mhx)
          XcYtot$mhy= rbind(XcYtot$mhy,XcY$mhy)
          XcYtot$mhz=rbind(XcYtot$mhz,XcY$mhz)
          XcYtot$mhxy=rbind(XcYtot$mhxy,XcY$mhxy)
          XcYtot$mhxz=rbind(XcYtot$mhxz,XcY$mhxz)
          XcYtot$mhxyz=rbind(XcYtot$mhxyz,XcY$mhxyz)
          
          # X causa Y --------------------------------------
          YcX <- symbolize(X[1:T-r],X[r+1:T],Y[1:T-r],m) 
          YcXtot$mhx= rbind(YcXtot$mhx,YcX$mhx)
          YcXtot$mhy= rbind(YcXtot$mhy,YcX$mhy)
          YcXtot$mhz= rbind(YcXtot$mhz,YcX$mhz)
          YcXtot$mhxy= rbind(YcXtot$mhxy,YcX$mhxy)
          YcXtot$mhxz= rbind(YcXtot$mhxz,YcX$mhxz)
          YcXtot$mhxyz= rbind(YcXtot$mhxyz,YcX$mhxyz)
          
     }
     EstXcY <- TEntropy(XcYtot)
     EstYcX <- TEntropy(YcXtot)
     
     # Net Transfer ------------------------------------
     EstNTE=EstXcY$TE-EstYcX$TE
     
     # Inicia el bootstrap -----------------------------
     # OJO: el procedimiento tsbootstrap permite sacar unamatrix de bootstraps
     # y así evitarnos un bulce -----------------------------  
     pb <- txtProgressBar(min=0, max=B, style=3) # Barra de progreso
     BTS <- data.frame(XcYtest=double(), YcXtest=double(),NTEtest=double())
     
     for (Bit in 1:B){
          
          XcYtot <- list(mhx= NULL,mhy= NULL,mhz= NULL,mhxy= NULL,mhxz= NULL,mhxyz= NULL)
          YcXtot <- list(mhx= NULL,mhy= NULL,mhz= NULL,mhxy= NULL,mhxz= NULL,mhxyz= NULL)
          
          for (j in 1:npanels) {
               # Sacamos la muestra bootstrap ----------------------------------------
               paneli <- subset(data,indice[[1]]==secciones[j])
               paneli <- subset(paneli,complete.cases(paneli))
               T=nrow(paneli);
               Y <- paneli[,noms[1]]
               X <- paneli[,noms[2]]
               OptBlockY <- b.star(Y)[1,"BstarSB"]
               OptBlockX <- b.star(X)[1,"BstarSB"]
               if (OptBlockX==0) {OptBlockX=1}
               if (OptBlockY==0) {OptBlockY=1}
               YBoot <- tsbootstrap(Y, nb=1, b = 5, type="stationary")
               XBoot <- tsbootstrap(X, nb=1, b = 5, type="stationary")
               
               # X causa Y --------------------------------------
               XcY <- symbolize(YBoot[1:T-r],YBoot[r+1:T],XBoot[1:T-r],m)
               XcYtot$mhx= rbind(XcYtot$mhx,XcY$mhx)
               XcYtot$mhy= rbind(XcYtot$mhy,XcY$mhy)
               XcYtot$mhz=rbind(XcYtot$mhz,XcY$mhz)
               XcYtot$mhxy=rbind(XcYtot$mhxy,XcY$mhxy)
               XcYtot$mhxz=rbind(XcYtot$mhxz,XcY$mhxz)
               XcYtot$mhxyz=rbind(XcYtot$mhxyz,XcY$mhxyz)
               
               # X causa Y --------------------------------------
               YcX <- symbolize(XBoot[1:T-r],XBoot[r+1:T],YBoot[1:T-r],m) 
               YcXtot$mhx= rbind(YcXtot$mhx,YcX$mhx)
               YcXtot$mhy= rbind(YcXtot$mhy,YcX$mhy)
               YcXtot$mhz= rbind(YcXtot$mhz,YcX$mhz)
               YcXtot$mhxy= rbind(YcXtot$mhxy,YcX$mhxy)
               YcXtot$mhxz= rbind(YcXtot$mhxz,YcX$mhxz)
               YcXtot$mhxyz= rbind(YcXtot$mhxyz,YcX$mhxyz)
          }
          
          # POner BTS como un data frame que recoge todos los bootstraps
          BtsXcY <- TEntropy(XcYtot)
          BtsYcX <- TEntropy(YcXtot)
          BtsNTE = BtsXcY$TE-BtsYcX$TE
          BTS[ nrow(BTS)+1, ] <- c(BtsXcY$TE, BtsYcX$TE, BtsNTE)
          setTxtProgressBar(pb,Bit)
     }
     
     close(pb) # End progress bar
     pvalXcY = sum(BTS$XcYtest>EstXcY$TE)/B
     pvalYcX = sum(BTS$YcXtest>EstYcX$TE)/B
     if (EstNTE>0) {
          pvalNTE1t = sum(BTS$NTEtest > EstNTE)/B
     } else {
          pvalNTE1t = sum(BTS$NTEtest < EstNTE)/B
     }
     pvalNTE2t = sum(abs(BTS$NTEtest)>abs(EstNTE))/B
     
     XcYtest = c(EstXcY$TE,pvalXcY) 
     names(XcYtest) <- c("Test","P-Value")
     YcXtest = c(EstYcX$TE,pvalYcX) 
     names(YcXtest) <- c("Test","P-Value")
     NETransfer = c(EstNTE,pvalNTE1t,pvalNTE2t)
     names(NETransfer) <- c("Test","PVal1s", "PVal2s")
     
     return(list(XcYtest = XcYtest, 
                 YcXtest = YcXtest, 
                 NETransfer = NETransfer, 
                 BTS = BTS))
}

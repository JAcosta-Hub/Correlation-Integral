################################
### Librerias

 library(mvtnorm)

 dmvnorm2 <- function(x1,x2,...){
   x=cbind(x1,x2)
   y=dmvnorm(x,...)
   return(y)
 }

##############################################
########### Funciones ########################
 
# Funcion de autocorrelación del ARMA(1,1)
 rho_arma11 = function(j=1, phi=0.5, theta=0.25){
   rhoj=(1+theta*phi)*(phi+theta)/(1+2*theta*phi+theta^2)*phi^(j-1)
   return(rhoj)
 }

## Calculo de la probabilidad (Volumen)
 vol <- function(rho=0.5, r=1, m=10000){
   r=r/sqrt(2*(1-rho))
   index=1:m
   x1 <- x2 <- runif(m, -r*1.1, r*1.1)
   x=cbind(x1,x2)[index[x1^2+x2^2<r],]
   y=dmvnorm(x)
   return(mean(y)*pi*r^2)
 }  

## Calculo del volumen exacto por Riemmann
 vol2 <- function(rho=0.5, r=1){
   rj=r/sqrt(2*(1-rho))
   g <- function(x,...){ # funcion auxiliar
       den=1-rho*sin(2*x)
       num=1-exp(-rj^2/(2*(1-rho^2))*den)
       res=num/den
       return(res)
     }
   theta=seq(0, 2*pi,l=201)
   del=theta[2]-theta[1]
   gt=g(theta)
   res=(mean(gt[-1])+mean(gt[-201]))*sqrt(1-rho^2)/2 # promedio suma por exceso y defecto
  return(res)
 }

## proyeccion
 kop=function(vol, r=1, rho=0.5){
   rj=r/sqrt(2*(1-rho))
   kop=sqrt(-2*log(2*vol/rj^2)-log(1-rho^2))
   return(kop)
 }
   
 proy=function(k,rho){
   x=seq(-k*0.999,k*0.999, by=0.01)
   R=matrix(c(1,rho,rho,1), nc=2)
   R1=solve(R)
   discr=R1[1,2]^2*x^2-R1[2,2]*(R1[1,1]*x^2-k^2)
   y1=-(R1[1,2]*x+sqrt(discr))/R1[2,2]
   y2=-(R1[1,2]*x-sqrt(discr))/R1[2,2]
   return(data.frame(x=x,y1=y1,y2=y2))
 }
 
 
 ######################################################
 ## Funcion que grafica y calcula el volumen para p=2
 
 plot_areaP2 <- function(j=1, r=1, phi=-0.5, theta=0.25, opcion=1,...){
   rho=rho_arma11(j=j, phi=phi, theta=theta)
   R=matrix(c(1,rho,rho,1), nc=2)
   rj <- r/sqrt(2*(1-rho))
   V1=vol(rho=rho, r=r) # probabilidad aproximada
   V2=vol2(rho=rho, r=r) # probabilidad exacta
   k=kop(vol=V2, r=r, rho=rho)
   
   xyp=proy(rho=rho, k=k)
   alt=dmvnorm2(xyp$x[1], xyp$y1[1], sigma=R)
   V3=pi*rj^2*alt # volumen estimado por el cilindro
   print(c(V1,V2,V3))
   
   #############################
   ## Grafica 3D 
   
   x=seq(-2.5, 2.5, l=201)
   y=expand.grid(x,x)
   
   phi <- seq(0, 2*pi, len = 251)
   xr <- rj * cos(phi)
   yr <- rj * sin(phi)
   
   aux=matrix(y[,1]^2+y[,2]^2<rj^2, nc=length(x), nr=length(x))
   z=outer(x, x, dmvnorm2, sigma=R)
   
   # Grafica del volumen calculado directo
   if(opcion==1){
     
     op <- par(bg = "white", mai=c(0.1,0.3,0.1,0.1) )
     persp(x, x, z,  theta = 20, phi = 20, expand = 0.5, col = "gray99", 
           ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "X",
           ylab = "Y", zlab = "", zlim=c(0,0.20) ) -> res
     lines(trans3d(xr, yr, dmvnorm2(xr, yr, sigma=R), res), col = "red", lwd = 2)
     for(i in seq(0, 1, by=0.01)) 
       lines(trans3d(xr, yr, dmvnorm2(xr, yr, sigma=R)*i, res),
             col = heat.colors(10, alpha=0.25), lwd = 2)
     # title(main=expression(P(abs(Z)<r[h])), xlab=paste('h =',j))
   }
   
   # grafica del volumen del cilindro de por el valor medio
   if(opcion==2){
     op <- par(bg = "white", mai=c(0.1,0.3,0.1,0.1) )
     persp(x, x, z,  theta = 20, phi = 20, expand = 0.5, col = "gray99", 
           ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "X",
           ylab = "Y", zlab = "", zlim=c(0,0.20) ) -> res
     for(i in seq(0, alt, by=0.005)) 
       lines(trans3d(xr, yr, i, res), col = cm.colors(10, alpha=0.4), lwd = 2)
     #  title(main=expression(P(abs(Z)<r[h])), xlab=paste('h =',j))
     lines(trans3d(c(xyp$x,rev(xyp$x)), c(xyp$y1, rev(xyp$y2)), alt, res2), col = "blue", lwd = 2)
     lines(trans3d(c(xyp$x,rev(xyp$x)), c(xyp$y1, rev(xyp$y2)), 0, res2), col = "blue", lwd = 2)
   }
 }

##########################################
############ Ejemplos ####################
 
## Valor de la densidad de la normal multivariada 
 
 rho=rho_arma11() # por defecto
 R=matrix(c(1,rho,rho,1), nc=2)
 R
 dmvnorm(cbind(1:2,1:2), sigma=R)
 dmvnorm2(1:2,1:2, sigma=R)
 
## Comparación de la probabilidad en términos de rho
 
 rho2=0.9
 plot(seq(-rho2, rho2,by=0.01),sapply(seq(-rho2, rho2,by=0.01), vol2),
      type="l", xlab=expression(rho), ylab="vol", col="red", ylim=c(0,2))
 points(seq(-rho2, rho2,by=0.01),sapply(seq(-rho2, rho2,by=0.01), vol),
        type="l", col="blue")
 abline(h=1, lty=2)
 legend("topleft", legend=c("Vol by Riemmann", "Vol by Monte-Carlo"), 
        col=c("red","blue"), lty=1, bg="gray90", inset=0.01)
  # el metodo de monte-carlo tiene problemas para correlaciones altas

 r=2
 vol(rho=0, r=r) # Monte-Carlo (mientras menor r, mejor la aproximación)
 vol2(rho=0, r=r) # Integral de Riemmann
 1-exp(-r^2/4) # valor exacto 
 
#############################################
## Graficos

png("arma1_p2.png", width = 1024, height=480, res=120)
par(mfrow=c(1,3))
plot_areaP2(j=1, r=2, phi=-0.5, theta=0.25, opcion=1)
plot_areaP2(j=2, r=2, phi=-0.5, theta=0.25, opcion=1)
plot_areaP2(j=3, r=2, phi=-0.5, theta=0.25, opcion=1)
dev.off()

png("arma2_p2.png", width = 1024, height=480, res=120)
par(mfrow=c(1,3))
plot_areaP2(j=1, r=2, phi=-0.5, theta=0.25, opcion=2)
plot_areaP2(j=2, r=2, phi=-0.5, theta=0.25, opcion=2)
plot_areaP2(j=3, r=2, phi=-0.5, theta=0.25, opcion=2)
dev.off()

plot_areaP2(j=200, r=2, phi=-0.5, theta=0.25, opcion=2)

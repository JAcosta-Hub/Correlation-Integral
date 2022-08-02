########################################
############ Funciones #################

## Calculo el área por Monte-carlo con Distribucion Normal
Area1 <- function(r, rho=0.5, p=1, m=1000){
   r=r/sqrt(2*(1-rho))
   y=c()
   n=0
   while(length(y)<m){
     n=n+1
     x=rnorm(1)
     if(abs(x)<r) y=c(y,1) 
   }
  return(m/n)
}

## Calculo del área por Monte-carlo con Dist Uniforme
Area2 <- function(r, rho=0.5, p=1, m=1000){
  r=r/sqrt(2*(1-rho))
  x=runif(m,-r,r)
  y=dnorm(x) 
  return(mean(y)*2*r)
}

## Calculo del area exacta
Area3=function(r=1, rho=0.5){return(2*pnorm(r/sqrt(2*(1-rho)))-1)}

## Calculo del punto del teorema del valor medio
PM=function(r=1, rho=0.5, ep=10^(-4)){
  Aj=Area3(r=r, rho=rho)
  r=r/sqrt(2*(1-rho))
  c0=r/2
  cj=r
  j=0
  while(abs(cj-c0)>ep|j<100){
    c0=cj
    dg=-2*r*c0*dnorm(c0)
    g=2*r*dnorm(c0)-Aj
    cj=c0-g/dg
    j=j+1
  }
  return(abs(cj))
}

# Funcion de autocorrelación del ARMA(1,1)
rho_arma11 = function(j=1, phi=0.5, theta=0.25){
  rhoj=(1+theta*phi)*(phi+theta)/(1+2*theta*phi+theta^2)*phi^(j-1)
  return(rhoj)
}

# Calculo y gráfica de E[delta_{i,j+1}]=P[|Z|<rj]
plot_area <- function(j=1, r=1, rho=0.5,...){
  x=seq(-4,4, by=0.01)
  y=dnorm(x)
  rho=rho_arma11(j=j,...)
  rj=r/sqrt(2*(1-rho))
  r0=PM(r=r, rho=rho)
  print(data.frame(rj=rj, cj=r0))

  A1=Area3(r=r, rho=rho) 
  A2=2*rj*dnorm(r0) # calcula area utilizando valor medio
  print(data.frame(Area_exacta=A1, Area_valor_medio=A2))
  
#################################
## Se crea el gráfico
xx <- c(x[abs(x)<rj],-x[abs(x)<rj])
yy <- c(y[abs(x)<rj], rep(0, sum(abs(x)<rj)))

xx0 <- c(-rj,rj,rj,-rj)
yy0 <- c(dnorm(-r0),dnorm(-r0),0,0)

plot(x, y, type = "l", xlab = "", ylab = "Density", xlim=c(-3,3), ylim=c(0,0.5))
abline(h=0)
polygon(xx, yy, col = "gray", border = "red", lwd=2)
polygon(xx0, yy0, col=cm.colors(10, alpha=0.2), border = "blue", lwd=1.5)
lines(c(r0,r0), c(0, dnorm(r0)), lty=2, lwd=1.5)
lines(-c(r0,r0), c(0, dnorm(r0)), lty=2, lwd=1.5)
title(main=expression(P(abs(Z)<r[h])), xlab=paste('h =',j))
legend("topleft", legend=c(paste("Exact Area =", round(A1,4)),
                            paste("MV Theorem =", round(A2,4))), 
       col=c("red","blue"), lty=1, bg="gray90", inset=0.01)
}


########################################
######### EjemploS #####################

## Calculo del Area

 Area1(r=1, rho=rho_arma11(j=1, phi=0.5, theta=0.9)) # Aprox 1
 Area2(r=1, rho=rho_arma11(j=1, phi=0.5, theta=0.9)) # Aprox 2
 Area3(r=1, rho=rho_arma11(j=1, phi=0.5, theta=0.9)) # Exacto

#r=1, phi=0.5 y theta=0.25 por defecto
rho=rho_arma11(j=1, phi=0.5, theta=0.9)
rho
par(mfrow=c(1,1))
plot_area(j=1, phi=-0.75, theta=0.5) 
plot_area(j=2, phi=-0.75, theta=0.5)
plot_area(j=3, phi=-0.75, theta=0.5)
plot_area(j=4, phi=-0.75, theta=0.5)
plot_area(j=5, phi=-0.75, theta=0.5)
plot_area(j=6, phi=-0.75, theta=0.5)
plot_area(j=100, phi=-0.75, theta=0.5)
plot_area(j=500, phi=-0.75, theta=0.5)

plot_area(j=1, phi=-0.5) 
plot_area(j=2, phi=-0.5)
plot_area(j=3, phi=-0.5)
plot_area(j=4, phi=-0.5)
plot_area(j=5, phi=-0.5)
plot_area(j=6, phi=-0.5)
plot_area(j=100, phi=-0.5)
plot_area(j=500, phi=-0.5)

plot_area(j=1, theta=-0.25) 
plot_area(j=2, theta=-0.25)
plot_area(j=3, theta=-0.25)
plot_area(j=4, theta=-0.25)
plot_area(j=5, theta=-0.25)
plot_area(j=6, theta=-0.25)
plot_area(j=100, theta=-0.25)
plot_area(j=500, theta=-0.25)



##################################
## grafico para el ejemplo
x=1:10
par(mfrow=c(1,1))
plot(x, rho_arma11(j=x, phi=-0.5, theta=0.25), type="h",
     main="phi=-0.5, theta=0.25", ylab=expression(rho(h)),
     xlab=expression(h), ylim=c(-0.65,0.65)); abline(h=0)


png("arma_p1.png", width = 1250, height=480, res=120)
par(mfrow=c(1,4))
plot_area(j=1, phi=-0.75, theta=0.5) 
plot_area(j=2, phi=-0.75, theta=0.5)
plot_area(j=10, phi=-0.75, theta=0.5)
plot_area(j=100, phi=-0.75, theta=0.5)
dev.off()


png("arma2_p1.png", width = 1250, height=480, res=120)
par(mfrow=c(1,4))
plot_area(j=1, phi=-0.5, theta=0.25, r=2) 
plot_area(j=2, phi=-0.5, theta=0.25, r=2)
plot_area(j=3, phi=-0.5, theta=0.25, r=2)
plot_area(j=4, phi=-0.5, theta=0.25, r=2)
dev.off()

plot_area(j=10, phi=-0.75, theta=0.5, r=2)
Area1(r=2, rho=0) # Aprox 1
Area2(r=2, rho=0) # Aprox 2
Area3(r=2, rho=0) # Exacto

Area1(r=1, rho=rho_arma11(j=1, phi=0.5, theta=0.9)) # Aprox 1
Area2(r=1, rho=rho_arma11(j=1, phi=0.5, theta=0.9)) # Aprox 2
Area3(r=1, rho=rho_arma11(j=1, phi=0.5, theta=0.9)) # Exacto



#############################################
## Comportamiento correlacion phi > theta

x=1:10
plot(x, rho_arma11(j=x, phi=0.5, theta=0.25), type="h", 
     main="phi=0.5, theta=0.25", ylab=expression(rho(h)), 
     xlab=expression(h), ylim=c(-0.65,0.65)); abline(h=0)

plot(x, rho_arma11(j=x, phi=-0.75, theta=0.5), type="h",
     main="phi=-0.5, theta=0.25", ylab=expression(rho(h)),
     xlab=expression(h), ylim=c(-0.65,0.65)); abline(h=0)

plot(x, rho_arma11(j=x, phi=0.5, theta=-0.25), type="h",
     main="phi=0.5, theta=-0.25", ylab=expression(rho(h)),
     xlab=expression(h), ylim=c(-0.65,0.65)); abline(h=0)

plot(x, rho_arma11(j=x, phi=-0.5, theta=-0.25), type="h",
     main="phi=-0.5, theta=-0.25", ylab=expression(rho(h)),
     xlab=expression(h), ylim=c(-0.65,0.65)); abline(h=0)


#############################################
## Comportamiento correlacion phi < theta

x=1:10
plot(x, rho_arma11(j=x, phi=0.25, theta=0.5), type="h", 
     main="phi=0.25, theta=0.5", ylab=expression(rho(h)), 
     xlab=expression(h), ylim=c(-0.6,0.6)); abline(h=0)

plot(x, rho_arma11(j=x, phi=-0.25, theta=0.5), type="h",
     main="phi=-0.25, theta=0.5", ylab=expression(rho(h)),
     xlab=expression(h), ylim=c(-0.6,0.6)); abline(h=0)

plot(x, rho_arma11(j=x, phi=0.25, theta=-0.5), type="h",
     main="phi=0.25, theta=-0.5", ylab=expression(rho(h)),
     xlab=expression(h), ylim=c(-0.6,0.6)); abline(h=0)

plot(x, rho_arma11(j=x, phi=-0.25, theta=-0.5), type="h",
     main="phi=-0.25, theta=-0.5", ylab=expression(rho(h)),
     xlab=expression(h), ylim=c(-0.6,0.6)); abline(h=0)



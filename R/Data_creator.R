

timefunction<- function (x)(  (x^2)/3 + x/2)
#'export
Data_sim=function(tfuntion=timefunction, depths=1:30,supp=15,thick=1,Fi=100,rho=0,errors=TRUE,crea.file=T,folder=T){
  dptS=depths[length(depths)]
  lambda=0.03114
  if (typeof(supp)=="character"){
    supported=supp(depths)
  }else{
    supported=rep(supp,length(depths))}
  ########Functions #########

  t1 =tfuntion#<- function (x)(  (x^2)/3 + x/2) #  8*x-20*sin(x/(1.*pi)) )#
  curve(t1,0,(dptS*.8))
  abline(h=250)
  abline(v=30)

  dt <- function(x){}
  body(dt) <- D(body(t1), 'x')
  if(typeof(rho)=="double"){
    rho <- function(x){
      1.5-.5*cos(x/30*pi)
    }
  }else{rho=rho}
  curve(rho,0,dptS)

  r.sed <- function(x){ rho(x)/dt(x)}
  curve(r.sed,0,dptS)

  C0<- function(x){  Fi/r.sed(x)}
  curve(C0,0,dptS)


  f <- function(x){

    return ((C0(x)*exp(-lambda*t1(x))))
  }
  curve(f,0,dptS)

  f2 <- function(x){
    return  ( (C0(x)*exp(-lambda*t1(x)))*rho(x))
  }
  curve(f2,0,dptS)

  ##########Sampling #######

  sample1=c()
  density=c()
  k=1
  for (i in depths){
    sample1=c(sample1,(integrate(f,(i-thick),i)$value+supported[k]))
    density=c(density,integrate(rho,(i-thick),i)$value)

    k=k+1
  }

  if(errors==TRUE){
    uncer=seq(10,1,length.out = length(sample1))
    sample1un=c()
    unc=c()
    suppo=abs(rnorm(length(supported),supported,rep(3,length(supported))) )
    unc.cons=3
    k=1
    for (i in sample1){
      sample1un=c(sample1un,i+rnorm(1,0,uncer[k]))
      k=k+1
    }
  }else{
    if(any(typeof(errors)=="double",typeof(errors)=="integer") ){
      if (length(errors)==1){
        sample1un=c()
        suppo=abs(rnorm(length(supported),supported,rep(3,length(supported))) )
        unc.cons=3
        k=1
        uncer=rep(errors,length(depths))
        for (i in sample1){
          sample1un=c(sample1un,i+rnorm(1,0,errors))
          k=k+1
        }
      } else{if(length(errors)==length(depths)){
        sample1un=c()
        suppo=abs(rnorm(length(supported),supported,rep(3,length(supported))) )
        unc.cons=3
        uncer=errors
        k=1
        for (i in sample1){
          sample1un=c(sample1un,i+rnorm(1,0,errors[k]))
          k=k+1
        }
    }else{
      print("Errors should either have the same lenth as the depths or be a constant")
        }
      }
    }else{
      sample1un=sample1
      uncer=rep(0,length(depths))
      suppo=abs(rnorm(length(supported),supported,rep(3,length(supported))) )
    }
  }


  plot(depths,sample1un,pch=16,cex=.6,main="Simulated data",xlab="Depth (cm)",ylab="Bq/kg",ylim=c(0,(max(sample1un)+max(uncer))) )
  segments(depths, sample1un+uncer, x1 = depths, y1 = sample1un-uncer)
  points(depths,suppo,pch=16,col="red",cex=.6)
  segments(depths, suppo+3, x1 = depths, y1 = suppo-3,col="red")


  sim_data=matrix(c(depths,(density/10),sample1un,uncer,rep(thick,length(depths)),suppo,rep(3,length(suppo)) ),nrow=length(depths),byrow = FALSE)
  
  
  if(crea.file==TRUE){
    if(folder==TRUE){
      if(!dir.exists("~/Plum/")){
        dir.create("~/Plum/")
      }
      ranind=floor(runif(1,10,99))
      newfolder=paste("~/Plum/","Simulation-",ranind,sep="" )
      if(!dir.exists(newfolder)){
        dir.create(newfolder)
        folder=paste(newfolder,sep="")
      }else{folder=paste(newfolder,sep="")}

    }else {ranind=floor(runif(1,10,99)) }
    
    Col.names=c("Depth (cm)","Density g/cm^3","210Pb (Bq/kg)","sd(210Pb)","Thickness (cm)","226Ra (Bq/kg)","sd(226Ra)")
    datname=paste("Simulation-",ranind,".csv",sep="")
    write.table(sim_data,file = paste(folder,"/",datname,sep=""),sep=",",col.names=Col.names, row.names=F)
  }
  print("Simulated data is located at" )
  print(folder)
  return(c(folder,datname))
  
}


#
#
# t <- function (x)(  (x^2)/3 + x/2) #  8*x-20*sin(x/(1.*pi)) )#
# supp<- function(x)(10+abs(20*sin(x/pi))+.9*x/10)
# curve(supp,0,30)
# ro=function(x)(.9+cos(x/5.80*pi)+.085*x)
# curve(ro,0,30)
# Data_sim(t,seq(2,30,1),"supp",thick=2,Fi=50,rho="ro")
#


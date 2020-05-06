

timefunction<- function (x)(  (x^2)/3 + x/2)
#'export
Data_sim=function(tfuntion=timefunction, Core.name=TRUE,depths=1:30,supp=15,thick=1,Fi.t=100,rho=0,errors=TRUE, error_supp=TRUE,crea.file=T,folder=T,epsi=.03,y.scat=1.5,sigmax=1.){
  dptS=depths[length(depths)]
  lambda=0.03114
  if (typeof(supp)=="character"){
    supported=supp(depths)
  }else{
    supported=rep(supp,length(depths))}
  ########Functions #########

  if (typeof(Fi.t) =="double"){
    Fi <- function(x)(Fi.t*(x/x))
  }else{Fi=Fi.t}
  curve(Fi,0,30)
  
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

  C0<- function(x){  Fi(x)/r.sed(x)}
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
  if(error_supp==TRUE){
    suppsd=supp*epsi*1.5 
  }else{suppsd=error_supp
  }

  if(errors==TRUE){
    sample1un=c()
    unc=c()
    suppo=abs(rnorm(length(supported),supported,rep(suppsd,length(supported))) )
    k=1
    for (i in sample1){
      sample1un=c(sample1un,i+rnorm(1,0,max(sigmax,epsi*y.scat*i)))
      k=k+1
    }
    uncer=epsi*y.scat*sample1un
    uncmin=which(uncer<sigmax)
    uncer[uncmin]=1
  }else{
    if(any(typeof(errors)=="double",typeof(errors)=="integer") ){
      if (length(errors)==1){
        sample1un=c()
        suppo=abs(rnorm(length(supported),supported,rep(suppsd,length(supported))) )
        k=1
        uncer=rep(errors,length(depths))
        for (i in sample1){
          sample1un=c(sample1un,i+rnorm(1,0,errors))
          k=k+1
        }
      } else{if(length(errors)==length(depths)){
        sample1un=c()
        suppo=abs(rnorm(length(supported),supported,rep(suppsd,length(supported))) )
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

      suppo=abs(rnorm(length(supported),supported,rep(suppsd,length(supported))) )
    }
  }

  
  plot(depths,sample1un,pch=16,cex=.6,main="Simulated data",xlab="Depth (cm)",ylab="Bq/kg",ylim=c(min(suppo)-suppsd,(max(sample1un)+max(uncer))) )
  segments(depths, sample1un+uncer, x1 = depths, y1 = sample1un-uncer)
  points(depths,suppo,pch=16,col="red",cex=.6)
  segments(depths, suppo+suppsd, x1 = depths, y1 = suppo-suppsd,col="red")
  lab <- c()
  for (i in depths){lab <- c(lab,paste0('Sim-',i))}

  sim_data=matrix(c(lab,depths,(density/10),sample1un,uncer,rep(thick,length(depths)),suppo,rep(suppsd,length(suppo)) ),nrow=length(depths),byrow = FALSE)

  
  if(crea.file==TRUE){
    if(folder==TRUE){
      if(!dir.exists("~/Bacon_runs/")){
        dir.create("~/Bacon_runs/")
      }
        if(Core.name==TRUE){
          ranind=floor(runif(1,10,99))
          newfolder=paste("~/Bacon_runs/","Simulation-",ranind,sep="" )
        }else{
          newfolder=paste0("~/Bacon_runs/",Core.name)
        }
      
      if(!dir.exists(newfolder)){
        dir.create(newfolder)
        folder=paste(newfolder,sep="")
      }else{folder=paste(newfolder,sep="")}

    }else {ranind=floor(runif(1,10,99)) }
    
    Col.names=c('Lab',"Depth (cm)","Density g/cm^3","210Pb (Bq/kg)","sd(210Pb)","Thickness (cm)","226Ra (Bq/kg)","sd(226Ra)")

    if(Core.name==TRUE){
    datname=paste("Simulation-",ranind,".csv",sep="")
    }else{
      datname=paste(Core.name,".csv",sep="")
    }
    
    
    write.table(sim_data,file = paste(folder,"/",datname,sep=""),sep=",",col.names=Col.names, row.names=F)
  }
  print("Simulated data is located at" )
  #print(folder)
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


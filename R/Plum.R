#' MCMC to obtain the posterior distribution of the Chronology.
#'
#' It takes the Data of eat each depth and creates samples from the
#' posterior distribution using the Twalk methodology.
#'
#' @param itera is the number of iterations that Twalk will do (1e+6 is the recomended)
#' @param Datos is the activity at each depth and organize starting from the surface.
#' @param sdDatos is the standar deviation of each activity.
#' @param Conti is a vector showing which data was sample discontinuosly (1 would mean the data was sample discontinusly and 0 means the data was sampled continously).
#' @param Bqkg is True when the data is in the form Bq/kg, if it False the data would be consider as ...

#' @export  
runPlum=function(Core.name=TRUE,folder=TRUE,iterations=2e+3,by=TRUE,
                 number_supported=FALSE,detection_limit=.1,Bqkg=TRUE,
                 Cs=TRUE,Sampledate=2017,Cs_date=1968,
                 memory_shape=4., memory_mean=.4,
                 acc_shape=1.5,acc_mean=15,fi_mean=50,fi_acc=2,
                 As_mean=10,As_acc=2,resolution=200,seeds=12345678,thin=30,burnin=7000){
  
  
  ##checks if data needs to be simulated
  if (Core.name==TRUE & folder==TRUE){
    folder=Data_sim()[1]
    temp1=unlist(strsplit(folder,split = "/"))
    Core.name=temp1[length(temp1)]
  }else if(Core.name==TRUE & typeof(folder)=="character"){
    temp1=unlist(strsplit(folder,split = "/"))
    Core.name=temp1[length(temp1)]
    #Core.name= "Simulation"
  }else if (folder==TRUE & typeof(Core.name)=="character"){
    folder=paste("~/Plum/",Core.name,sep="")
  }
  #remove(temp1)
  # if(Core.name==TRUE){
  #   temp1=unlist(strsplit(folder,split = "/"))
  #   Core.name=temp1[length(temp1)]
  # }  
  folder=paste(normalizePath(folder),"/",sep="")
  dir.create(paste(folder,"Results ",Core.name,sep = ""))

  Lead=read.table(paste(folder,Core.name,".csv",sep=""),sep=",",header = T)
  write.table(Lead,paste(folder,"Results ",Core.name,"/",Core.name,".csv",sep=""),sep=",",col.names = F,row.names = F)
  Lead=read.table(paste(folder,"Results ",Core.name,"/",Core.name,".csv",sep=""),sep=",")
  
  if (as.character(by)==TRUE){
   by=(Lead[length(Lead[,1]),1])/20#25
  }
  




print("working folder is")
print(folder)
if(number_supported==FALSE){
  if(length(Lead[1,])==5){
    cat("You do not have 226Ra cocentrations.")
    usemod  =1
    n.check=check.equi(Lead)
    print(n.check[1])
    usrresp=1
    while(!(usrresp=="Yes"||usrresp=="No"||usrresp=="no"||usrresp=="yes")){
      cat("Are you sure these data represent the supported 210Pb for this site?")
      usrresp=readline( "Yes or No \n ")
      if(usrresp=="Yes"||usrresp=="No"||usrresp=="no"||usrresp=="yes"){
        if(usrresp=="Yes"||usrresp=="yes"){
          number_supported=n.check[1]}
        if(usrresp=="No"||usrresp=="no"){
          checksupp="notoy"
          while(typeof(checksupp)!="double"){
            cat("Please indicate how many data points whould be use for estimating the supported 210Pb")
            number_supported=scan("",n=1)

            }
        }
      }
    }
    usemod=1
  }else if(length(Lead[1,])==7){
    cat("You have 226Ra data. \n")
    plot(Lead$V1,Lead$V6,pch=16,ylim=c(min(Lead$V6-Lead$V7),max(Lead$V6+Lead$V7)), ylab="Concentration of 226Ra", xlab="Depth (cm)")
    segments(Lead$V1, Lead$V6-Lead$V7, x1 = Lead$V1, y1 = Lead$V6+Lead$V7)
    cat("Plum can assume to have a constant supported 210Pb and use the 226Ra data to infer this one value\n")
    cat("Plum can also assume individual supported 210Pb per data point.\n
        It is important to consider that this will greatly increases the computing time and it should only be use when clear patters are observed in the 226Ra data.\n")
    cat("\n If you want to use the constant supported 210Pb press 1, if you want to use the individual 210Pb press 2\n")
    usemod=0
    while(!(usemod==1||usemod==2)){
      usemod=scan("",n=1)
    }
    if (usemod==1){number_supported=0}
    print(number_supported)
    }
}else{usemod=1}



settings.file=c(Core.name,iterations, by, number_supported, detection_limit, memory_shape,memory_mean,
                acc_shape,acc_mean,fi_mean,fi_acc,As_mean,As_acc,resolution,seeds,thin,burnin,usemod)
settings.file=matrix(settings.file,ncol=1)
write.table(x = settings.file,file = paste(folder,"settings.txt",sep=""),sep = ",",col.names = F,row.names = F)  







modirec=path.package("Plum", quiet = T)


#twalk=paste(modirec,"/","pytwalk.py",sep="")
if (usemod==1){
  MCMC=paste(modirec,"/","ModelMCMC.py",sep="")
}else if(usemod==2){
  MCMC=paste(modirec,"/","ModelMCMC_Ra.py",sep="")
}


python.load(MCMC)
#python.load(twalk)


python.call("plumMCMC",folder,Core.name,FALSE,    number_supported   ,    detection_limit   ,  iterations,
   by ,memory_shape     ,memory_mean    ,acc_shape       ,acc_mean,fi_mean,fi_acc,As_mean,As_acc
   ,resolution,seeds,thin,burnin,Bqkg,Cs,Sampledate,Cs_date)

Lead=read.table(paste(folder,Core.name,".csv",sep=""),sep=",")

if(length(Lead[1,])==5){
  Col.names=c("Depth","Density","210Pb","sd(210Pb)","Thickness")
}else if(length(Lead[1,])==7){
  Col.names=c("Depth","Density","210Pb","sd(210Pb)","Thickness","226Ra","sd(226Ra)")
}

#write.table(Lead,paste(folder,Core.name,".csv",sep=""),sep=",",col.names = Col.names,row.names = F)

##############
Data=paste(Core.name,".csv",sep="")

Output=read.table(paste(folder,"Results ",Core.name,"/Results_output.csv",sep=""),sep=",",header = T)

num_var=length(Output[0,])

dev.new()
plot(as.numeric(Output[-1,num_var]),type="l",main="Energy",xlab="",ylab="")


pdf(paste(folder,paste('Chronologylines ',Core.name,'.pdf',sep=""),sep=""))
chronologylines(folder = folder)
dev.off()


pdf(paste(folder,paste('Energy ',Core.name,'.pdf',sep=""),sep=""))
plot(as.numeric(Output[-1,num_var]),type="l",main="Energy",xlab="",ylab="")
dev.off()

pdf(paste(folder,paste('Chronology ',Core.name,'.pdf',sep=""),sep=""))
fullchronology(folder = folder)
dev.off()


par(mfrow=c(1,1))
dev.new()
fullchronology(folder = folder)
intervalfile=by_cm(folder,Core.name)
write.csv(x = intervalfile,file = paste0(folder,"ages.csv"))


}

#################Check equilibrium #########################
#' @export
check.equi = function (rawdat){
  rawdata=rawdat[,3]
  rawsd=rawdat[,4]
  deps=rawdat[,1]
  cat("Because linear regression needs at least three data points, this function will start checking for equilibrium starting from the last three data points\n")
  lendat=length(rawdata)
  numdat=as.integer(.5*length(rawdata))
  usedat=rawdata[(lendat-3):lendat]
  usesd=rawsd[(lendat-3):lendat]
  usex=1:length((lendat-3):lendat)
  usereg= lm(usedat ~ usex, weights=1/(usesd^2))
  reg=coef(summary(usereg))[2,4]
  est=coef(summary(usereg))[1,1]
  coe=3
  for (i in 1:numdat){
    usedat=rawdata[(lendat-3-i):lendat]
    usesd=rawsd[(lendat-3-i):lendat]
    usex=1:length((lendat-3-i):lendat)
    usereg= lm(usedat ~ usex, weights=1/(usesd^2))
    reg1=coef(summary(usereg))[2,4]
    est1=coef(summary(usereg))[1,1]
    if(reg1>reg){ reg=reg1;coe=(3+i);est=est1 }
  }

  cat(paste( "The regression process proposes the use of the last", as.integer(coe), " data points as estimates of the supported activity, with a p-value of" ,round((reg),3)) )
  cat(".\n")
    plot(deps,rawdata,pch=16,xlab="Depth",ylab="210Pb")
  points(deps[(lendat-coe+1):lendat],rawdata[(lendat-coe+1):lendat],col="red",pch=16,xlab="Depth (cm)",ylab="210Pb concentration")
  abline(h=est,lty=2) #mean(rawdata[(lendat-coe+1):lendat])
  return (c(coe,reg1))
}

#' @export
check.qui.Ra = function(rawdat){
  `226Ra`=rawdat[,6]
  rawsd=rawdat[,7]
  deps=rawdat[,1]
  tests=shapiro.test(`226Ra`)
  print(tests)
  if(tests$p.value<=.05){
    cat("Because the p-value of the normality test is smaller than .05 it cannot be assumed that the 226Ra data comes from a single distribution.\n")
    return(2)
  }else{
    cat("Because the p-value of the normality test is bigger than .05 it can be assumed that the 226Ra data comes from a single distribution.\n")
    return(1)
  }
}

#' @export
by_cm=function(folder,Core.name){
  #paste(folder,"Results ",Core.name,"/",Core.name,".csv",sep="")
  intv=read.table(paste(folder,"Results ",Core.name,"/intervals.csv",sep=""),sep=",")
  agedmodl=approx(c(0,intv[,1]),c(0,intv[,2]),c(1:max(intv[,1])))$y
  agedmodm=approx(c(0,intv[,1]),c(0,intv[,3]),c(1:max(intv[,1])))$y
  agedmodu=approx(c(0,intv[,1]),c(0,intv[,4]),c(1:max(intv[,1])))$y
  Age_max=matrix(c(c(1:max(intv[,1])),agedmodl,agedmodm,agedmodu),ncol = 4,byrow=F)
  colnames(Age_max) <- c("depth","min","mean","max")
  return(Age_max)
}

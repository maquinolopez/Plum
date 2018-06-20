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
runPlum=function(Data=TRUE,folder=TRUE,iterations=1.5e+3,by=1.5,number_supported=FALSE,
                 detection_limit=.05,memory_shape=4., memory_mean=.7,
                 acc_shape=1.5,acc_mean=20,fi_mean=50,fi_acc=2,
                 As_mean=20,As_acc=2,resolution=200,seeds=12345678,thin=30,burnin=10000){
  ##checks if data needs to be simulated
  if (Data==TRUE & folder==TRUE){
    folder=Data_sim() 
    Data=folder[2]
    folder=folder[1]
  }else if(Data==TRUE & typeof(folder)=="character"){
    folder=Data_sim(folder =folder) 
    Data=folder[2]
    folder=folder[1]
  }
    
  Lead=read.table(paste(folder,Data,sep=""),sep=",",header = T)
  write.table(Lead,paste(folder,Data,sep=""),sep=",",col.names = F,row.names = F)
  Lead=read.table(paste(folder,Data,sep=""),sep=",")

  




folder=paste(normalizePath(folder),"/",sep="")


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
    cat("Plum can assum to have a constant supported 210Pb and use the 226Ra data to infer this one value\n")
    cat("Plum can also assum individual supporeted 210Pb per data point.\n
        It is important to consider that this will greatly increses the computing time and it should only be use when clear parters are observed in the 226Ra data.\n")
    cat("\n If you want to use the constant supported 210Pb press 1, if you want to use the individual 210Pb press 2\n")
    usemod=0
    while(!(usemod==1||usemod==2)){
      usemod=scan("",n=1)
    }
    if (usemod==1){number_supported=0}
    print(number_supported)
    }
}else{usemod=1}

modirec=path.package("Plum", quiet = T)


#twalk=paste(modirec,"/","pytwalk.py",sep="")
if (usemod==1){
  MCMC=paste(modirec,"/","ModelMCMC.py",sep="")
}else if(usemod==2){
  MCMC=paste(modirec,"/","ModelMCMC_Ra.py",sep="")
}


python.load(MCMC)
#python.load(twalk)
dir.create(paste(folder,"Results",sep = ""))

python.call("plumMCMC",folder,Data,FALSE,    number_supported   ,    detection_limit   ,  iterations,
   by ,memory_shape     ,memory_mean    ,acc_shape       ,acc_mean,fi_mean,fi_acc,As_mean,As_acc
   ,resolution,seeds,thin,burnin)


##############


Ages=read.table(paste(folder,"Results/dates.csv",sep=""),sep=" ")
intervals=read.table(paste(folder,"Results/intervals.csv",sep=""),sep=",")
Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
Output=read.table(paste(folder,"Results/Results_output.csv",sep=""),sep=",")
Plotval=read.table(paste(folder,"Results/Graphs.csv",sep=""),sep=",")
Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
num_var=length(Output[0,])

maxA=max(Ages[,length(Ages[1,])])+.10
ageSeq=seq(from=0,to=maxA,maxA/resolution)
deptsSeq=seq(from=0,to=Depths[length(Depths)],Depths[length(Depths)]/resolution)
deptsSeq=deptsSeq
diffSep=(deptsSeq[2]-deptsSeq[1])/2
TotSeq=length(Ages[,1])
X11()
plot(as.numeric(Output[-1,num_var]),type="l",main="Energy",xlab="",ylab="")


pdf(paste(folder,'Chronologylines.pdf',sep=""))
plot(-1,-1, xlim=c(0,Depths[length(Depths)]),ylim = c(0, maxA),xlab = "Depth (cm)",ylab="Age (years)" )

for (i2 in 1:(resolution-1)){
    for (i in 1:(resolution-1) ){
      rect(deptsSeq[i2], ageSeq[i], deptsSeq[i2+1], ageSeq[i+1], density = NA, border = NA,
      col = gray((1-Plotval[i2,i])))
    }
}
lines(Depths,c(0,intervals[,2]),type="l", lty=2, lwd=1)
lines(Depths,(c(0,intervals[,4])),type="l", lty=2, lwd=1)
lines(Depths,(c(0,intervals[,3])),type="l", lty=2, lwd=1)
dev.off()


# pdf(paste(folder,'Slopes.pdf',sep=""))
# maxS=max(Slopes)+.10
# plot(-10,-10, xlim=c(0,Depths[length(Depths)]),ylim = c(0, maxS),xlab = "Depth (cm)",ylab="Slopes (Accumulations)" )
# for (i in 1:(iterations-1)){
#   lines(Depths,as.numeric(c(0,Slopes[i,])),type="l",col=rgb(0,0,0,.01), lwd=2)
# }
# dev.off()
# pdf(paste(folder,'Chronology.pdf',sep=""))
# plot(Depths,c(0,Ages[2,]),type="l",col=rgb(0,0,0,.01), lwd=2,ylim = c(0,max(Ages[,length(Ages[1,])])),xlab = "Depth (cm)",ylab="Age (years)")
# for (i in 1:(iterations-1)){
#   lines(Depths,c(0,Ages[i,]),type="l",col=rgb(0,0,0,.01), lwd=2)
# }
# lines(Depths,c(0,intervals[,2]),type="l", lty=2, lwd=1)
# lines(Depths,(c(0,intervals[,4])),type="l", lty=2, lwd=1)
# lines(Depths,(c(0,intervals[,3])),type="l", lty=2, lwd=1)
# dev.off()


# 
# 
# pdf(paste(folder,'Fi.pdf',sep=""))
# plot(as.numeric(Output[-1,1]),type="l",main="fi",xlab="",ylab="")
# dev.off()
# pdf(paste(folder,'Supported.pdf',sep=""))
# plot(as.numeric(Output[-1,2]),type="l",main="Supported Act",xlab="",ylab="")
# dev.off()
pdf(paste(folder,'Energy.pdf',sep=""))
plot(as.numeric(Output[-1,num_var]),type="l",main="Energy",xlab="",ylab="")
dev.off()
# pdf(paste(folder,'Memory.pdf',sep=""))
# plot(as.numeric(Output[-1,3]),type="l",main="Memory",xlab="",ylab="")
# dev.off()
# 
# pdf(paste(folder,'Fi hist.pdf',sep=""))
# hist(as.numeric(Output[-1,1]),breaks=50,probability=T,main="fi",xlab="")
# dev.off()
# pdf(paste(folder,'Supported hist.pdf',sep=""))
# hist(as.numeric(Output[-1,2]),breaks=50,probability=T,main="Supported Act",xlab="")
# dev.off()
# pdf(paste(folder,'Energy hist.pdf',sep=""))
# hist(as.numeric(Output[-1,num_var]),breaks=50,main="Energy",xlab="")
# dev.off()
# pdf(paste(folder,'Memory hist.pdf',sep=""))
# hist(as.numeric(Output[-1,3]),breaks=50,main="Memory",probability=T,xlim=c(0,1),xlab="")
# lines(seq(0,100,.01),dbeta(seq(0,100,.01),4,1.5),col="red")
# dev.off()
# 
# 
# pdf(paste(folder,'acc hist.pdf',sep=""))
# hist(as.numeric(unlist(Output[-1,-c(1,2,3,num_var)])),breaks=50,main="Acc",probability=T,xlab="")
# lines(seq(0,100,.5),dgamma(seq(0,100,.5),1.5,scale=20/1.5),col="red")
# dev.off()




pdf(paste(folder,'Chronology.pdf',sep=""))
fullchronology(folder = folder,Data = Data,supp_type = usemod,resolution = resolution,
               memory_shape = memory_shape,memory_mean =  memory_mean,
               acc_shape = acc_shape,acc_mean = acc_mean,
               fi_mean = fi_mean,fi_acc = fi_acc,
               As_mean = As_mean,As_acc = As_acc)
dev.off()


par(mfrow=c(1,1))
X11()
fullchronology(folder = folder,Data = Data,supp_type = usemod,resolution = resolution,
               memory_shape = memory_shape,memory_mean =  memory_mean,
               acc_shape = acc_shape,acc_mean = acc_mean,
               fi_mean = fi_mean,fi_acc = fi_acc,
               As_mean = As_mean,As_acc = As_acc)

Lead=read.table(paste(folder,Data,sep=""),sep=",",header = T)
Col.names=c("Depth (cm)","Density g/cm^3","210Pb (Bq/kg)","sd(210Pb)","Thickness (cm)","226Ra (Bq/kg)","sd(226Ra)")
write.table(Lead,paste(folder,Data,sep=""),sep=",",col.names = Col.names,row.names = F)

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


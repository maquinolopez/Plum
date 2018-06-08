#' @export

fullchronology= function(folder,Data,resolution=200,supp_type=1,
                         memory_shape=4., memory_mean=.7,
                         acc_shape=1.5,acc_mean=20,
                         fi_mean=50,fi_acc=2,
                         As_mean=20,As_acc=2){


  par(mfrow=c(1,1))
  Ages=read.table(paste(folder,"Results/dates.csv",sep=""),sep=" ")
  intervals=read.table(paste(folder,"Results/intervals.csv",sep=""),sep=",")
  Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
  Output=read.table(paste(folder,"Results/Results_output.csv",sep=""),sep=",")
  Lead=read.table(paste(folder,Data,sep=""),sep=",")
  Plotval=read.table(paste(folder,"Results/Graphs.csv",sep=""),sep=",")
  Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
  num_var=length(Output[0,])
  if(length(Lead[1,])==5){
    Ran=1  } else if(length(Lead[1,])==7){
      if(supp_type==1){Ran=1}else{
      Ran=length(Lead[,1])}}

  maxA=max(Ages[,length(Ages[1,])])+.10
  ageSeq=seq(from=0,to=maxA,maxA/resolution)
  deptsSeq=seq(from=0,to=Depths[length(Depths)],Depths[length(Depths)]/resolution)
  deptsSeq=deptsSeq
  diffSep=(deptsSeq[2]-deptsSeq[1])/2
  TotSeq=length(Ages[,1])
  memory_shape2=(memory_shape*(1-memory_mean) )/memory_mean


  layout(matrix(c(1,1,2,2,3,3,4,4,
                  1,1,2,2,3,3,4,4,
                  5,5,5,5,5,5,5,5,
                  5,5,5,5,5,5,5,5,
                  5,5,5,5,5,5,5,5,
                  5,5,5,5,5,5,5,5,
                  5,5,5,5,5,5,5,5), 7, 8, byrow = TRUE))

  d <- density(as.numeric(Output[-1,(Ran+2)]))
  plot(d, xlab="",main="Memory",ylab = "",xlim=c(0,1),xaxs="i",yaxs="i")
  polygon(d, col=gray(.6))
  lines(seq(0,1,.01),dbeta(x = seq(0,1,.01),shape1 = memory_shape,shape2 = memory_shape2),col="green")


  d <- density(as.numeric(unlist(Output[-1,-c(c(1:(Ran+2)),num_var)])))
  plot(d,xlab="",main="Acc",ylab = "",xlim=c(0,80),xaxs="i",yaxs="i")
  polygon(d, col=gray(.6))
  lines(seq(0,100,.5),dgamma(seq(0,100,.5),shape=acc_shape,scale=acc_mean/acc_shape),col="green")

  if(Ran==1){
  d <- density(as.numeric(Output[-1,2]))
  plot(d,xlab="",main="Supported Act",ylab = "",xaxs="i",yaxs="i")
  polygon(d, col=gray(.6))
  lines(seq(0,100,.05),dgamma(seq(0,100,.05),shape=As_acc,scale=As_mean/As_acc),col="green")
  }else {
    min1=min(as.numeric(unlist(Output[-1,2:(Ran+1)])))+.1
    max1=max(as.numeric(unlist(Output[-1,2:(Ran+1)])))+.1
    plot(0,0,xlim=c(Lead[1,1],Lead[Ran,1]),ylim=c(min1,max1),xlab="Depth (cm)",ylab="",main="Supported 210Pb",xaxs="i",yaxs="i")
    for (k in 1:Ran) {
      points(rep(Lead[k,1],length(Output[-1,1+k])), Output[-1,1+k],pch=19,col=rgb(0,0,0,.03) )
      points(Lead[k,1],Lead[k,6],col="red",pch=18)

    }
  }

  d <- density(as.numeric(Output[-1,1]))
  plot(d,xlab="",main="Supply of 210Pb",ylab = "",xaxs="i",yaxs="i")
  polygon(d, col=gray(.6))
  lines(seq(0,350,.05),dgamma(seq(0,350,.05),shape=fi_acc,scale=fi_mean/fi_acc),col="green")

  plot(-1,-1, xlim=c(0,Depths[length(Depths)]),ylim = c(0, maxA),xlab = "Depth (cm)",ylab="Age (years)" ,
       xaxs="i",yaxs="i",main= gsub("\\..*","",Data))
  for (i2 in 1:(resolution-1)){
    for (i in 1:(resolution-1) ){
      rect(deptsSeq[i2], ageSeq[i], deptsSeq[i2+1], ageSeq[i+1], density = NA, border = NA,
           col = gray((1-Plotval[i2,i])))
    }
  }
  lines(Depths,c(0,intervals[,2]),type="l", lty=2, lwd=1,col="red")
  lines(Depths,(c(0,intervals[,4])),type="l", lty=2, lwd=1,col="red")
  lines(Depths,(c(0,intervals[,3])),type="l", lty=2, lwd=1,col="red")
  for (i in 1:length(Lead[,1])){
    rug(Lead[i,1], lwd = Lead[i,5],col = "blue")
  }

  par(mfrow=c(1,1))
}


#' @export
chronologylinesP= function(folder,Data,...){
  Ages=read.table(paste(folder,"Results/dates.csv",sep=""),sep=" ")
  intervals=read.table(paste(folder,"Results/intervals.csv",sep=""),sep=",")
  Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
  Output=read.table(paste(folder,"Results/Results_output.csv",sep=""),sep=",")
  Lead=read.table(paste(folder,Data,sep=""),sep=",")
  Plotval=read.table(paste(folder,"Results/Graphs.csv",sep=""),sep=",")
  Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
  num_var=length(Output[0,])
  iterations=length(Ages[,1])


  plot(Depths,c(0,Ages[2,]),type="l",col=rgb(0,0,0,.01), lwd=2,ylim = c(0,max(Ages[,length(Ages[1,])])),
       xlab = "Depth (cm)",ylab="Age (years)",...)
  for (i in 1:(iterations-1)){
    lines(Depths,c(0,Ages[i,]),type="l",col=rgb(0,0,0,.01), lwd=2)
  }
  lines(Depths,c(0,intervals[,2]),type="l", lty=2, lwd=1)
  lines(Depths,(c(0,intervals[,4])),type="l", lty=2, lwd=1)
  lines(Depths,(c(0,intervals[,3])),type="l", lty=2, lwd=1)
  for (i in 1:length(Lead[,1])){
    rug(Lead[i,1], lwd = Lead[i,5],col = "blue")
  }


}


#' @export
chronologylines= function(folder,Data,...){
  Ages=read.table(paste(folder,"Results/dates.csv",sep=""),sep=" ")
  intervals=read.table(paste(folder,"Results/intervals.csv",sep=""),sep=",")
  Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
  Output=read.table(paste(folder,"Results/Results_output.csv",sep=""),sep=",")
  Lead=read.table(paste(folder,Data,sep=""),sep=",")
  Plotval=read.table(paste(folder,"Results/Graphs.csv",sep=""),sep=",")
  Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
  num_var=length(Output[0,])
  iterations=length(Ages[,1])


plot(Depths,c(0,Ages[2,]),type="l",col=rgb(0,0,0,.01), lwd=2,ylim = c(0,max(Ages[,length(Ages[1,])])),
     xlab = "Depth (cm)",ylab="Age (years)",main= gsub("\\..*","",Data),...)
for (i in 1:(iterations-1)){
  lines(Depths,c(0,Ages[i,]),type="l",col=rgb(0,0,0,.01), lwd=2)
}
lines(Depths,c(0,intervals[,2]),type="l", lty=2, lwd=1)
lines(Depths,(c(0,intervals[,4])),type="l", lty=2, lwd=1)
lines(Depths,(c(0,intervals[,3])),type="l", lty=2, lwd=1)
for (i in 1:length(Lead[,1])){
  rug(Lead[i,1], lwd = Lead[i,5],col = "blue")
}


}



#' @export
slopes= function(folder,Data){
    Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
    Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
    iterations=length(Slopes[,1])
  maxS=max(Slopes)+.10
  plot(-10,-10, xlim=c(Depths[length(Depths)]),ylim = c(0, maxS),xlab = "Depth (cm)",ylab="Slopes (Accumulations)" )
  for (i in 1:(iterations-1)){
    lines(Depths,as.numeric(c(Slopes[i,])),type="l",col=rgb(0,0,0,.01), lwd=2)
  }
  for (i in 1:length(Lead[,1])){
    rug(Lead[i,1], lwd = Lead[i,5],col = "blue")
  }
}

#' @export
ageof=function(x,interval=.95,folder){
  if(x!=0){
    Ages=read.table(paste(folder,"Results/dates.csv",sep=""),sep=" ")
    Depths=as.numeric(read.table(paste(folder,"Results/depths.csv",sep=""),sep=",") )
    Slopes=read.table(paste(folder,"Results/Slopes.csv",sep=""),sep=",")
    depfix=which(Depths<x)
    depfix=depfix[length(depfix)]
    m2=Slopes[,depfix]
    sumages=c()
    if(depfix!=1){
      for (i in 1:length(Ages[,1])){
        sumages=c(sumages, Ages[i,(depfix-1)] + Slopes[i,depfix]* (x-Depths[depfix]) )

      }
    }else{sumages= Slopes[,depfix]* (x)}
    n=length(Ages[,1])
    mean1=mean(sumages)
    inter=(1-interval)/2
    lim1=sort(sumages)[as.integer(n*inter)]
    lim2=sort(sumages)[as.integer(n*(inter+interval))]

    d<- density(sumages)
    plot(d, main=paste("Age at ",x, "cm"),ylab = "",xlab = "Age",xaxs="i",yaxs="i")
    polygon(d, col=gray(.6))
    points(x=lim1,y=(max(d$y)*.015),col="red",pch="(")
    points(x=lim2,y=(max(d$y)*.015),col="red",pch=")")
    points(x=mean1,y=(max(d$y)*.015),col="red",pch=16)
    print(paste("the age of depth",x,"is between") )
    print(c(lim1,lim2))
    print(paste("with a ", interval, "%"," condifence interval and a mean of:",sep = ""))
    print(mean1)


  }else{
    print("For depth 0, the age is equal to the collection date.")
  }


}





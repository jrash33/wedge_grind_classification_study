rm(list=ls())

#functions
LoadLibraries=function(){
  library(svDialogs)
  library(fda)
  library(Funclustering)
  library(R.matlab)
  library(ggplot2)
  library(date)
  library(data.table)
  library(readr)
  library(plyr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(zoo)
  library(outliers)
  library(rlist)
  library(tuple)
  print("The Libraries have been loaded!")
}

LoadLibraries()

#---------------------------------- FUNCTIONS -------------------------------------------------

#function to read in file (skip first 6 lines/skip NaNs/fill blanks/no headers)
freadFUN=function(fileX){
  x=fread(fileX,skip=6,na.strings=c("","NaN"),sep="\t",
        fill=TRUE, header=FALSE)
  return(x)
}

#function returns logical list showing which swings don't have ball hat captures
findbadhatFUN=function(listX){
  ballbackloc=unlist(listX[1,135])
  badhatlogical=(typeof(ballbackloc)=="double")
  return(badhatlogical)
}

#function returns logical list showing which swings have messed up impact location
#bad swings have ball hat coordinates as the last marker
#if I subtract the first 2 rows of the last marker and get zero, I know it's a bad swing
badswingFUN=function(listX){
  col1=unlist(listX[1:2,138])
  diffcol1=col1[1]-col1[2]
  difflogical=(diffcol1==0 || is.na(diffcol1==TRUE))
  return(difflogical)
}


#function to delete all points after specified point (input is a list of impact points for all swings)
DeletePostImpactFUN=function(MoCapALL,cutpointsX){
  startpoint=cutpointsX+1
  endpoint=nrow(MoCapALL)
  MoCapALL=MoCapALL[-startpoint:-endpoint,]
  return(MoCapALL)
}

#function to return the impact location of each swing (last frame seen prior to impact)
#currently set to last frame before impact
impactFUN=function(listX){
  ballhat=unlist(listX[1,135])
  V_ImpactLocation=listX[,138]
  #val=V_ImpactLocation-ballhat
  
  maxloc=which.max(matrix(unlist(V_ImpactLocation),nrow = nrow(V_ImpactLocation), ncol = length(V_ImpactLocation)))
  
  #new addition 10/29
  if (length(which(is.na(V_ImpactLocation[(maxloc-30):maxloc])))>25){
    listX=DeletePostImpactFUN(listX,maxloc-1)
    V_ImpactLocation=listX[,138]
    maxloc=which.max(matrix(unlist(V_ImpactLocation),nrow = nrow(V_ImpactLocation), ncol = length(V_ImpactLocation)))
  }
  
  if (length(which(is.na(V_ImpactLocation[maxloc:(maxloc+50)])))<20){
    listX=DeletePostImpactFUN(listX,maxloc-1)
    V_ImpactLocation=listX[,138]
    maxloc=which.max(matrix(unlist(V_ImpactLocation),nrow = nrow(V_ImpactLocation), ncol = length(V_ImpactLocation)))
  }
  #end of new addition
  
  diffvals=V_ImpactLocation-ballhat
  
  j=maxloc
  while (diffvals[j]>=0 || is.na(diffvals[j])==TRUE){
    j=j-1
  }
  
  #0 frames past impact
  j=j+0
  
  return(j)
}

#function that combines lists in a list into one matrix of ncol=number of lists/nrow=MatrixSize
CombineListMatFUN=function(List,MatrixSize){
  add.col=function(vector,newsize){
    length(vector)=newsize
    return(vector)
  }
  delete.col=function(matrix,newsize){
    matrix=matrix[-(newsize+1):-nrow(matrix),]
    return(matrix)
  }

  X=lapply(List,function(List){mapply(add.col,vector=List,newsize=MatrixSize[which.max(MatrixSize)])})
  AllSwingsMat=mapply(delete.col,matrix=X,newsize=MatrixSize)
  AllSwingsMat[[which.max(MatrixSize)]]=X[[which.max(MatrixSize)]]
  return(AllSwingsMat)
}

#function takes a list of matrices and registers each list to max point in curve OR impact
#if val=1->register to max point/if val=0->register to last point (LATER)
registerFUN=function(AllSwingsMat,SizeVec,MatrixSize,val){
  
  if (val==1){
  maxloc=lapply(AllSwingsMat,function(x){apply(x,2,function(x){which.max(x)})})
  }else{
    maxloc=SizeVec
  }
  #pinpoint longest swing by finding the swing with the furthest max location
  OffsetVal=lapply(maxloc,function(x){x[[which.max(x)]]})
  #now find the offset value for every swing
  Offset=mapply(function(maxloc,OffsetVal){OffsetVal-unlist(maxloc)},maxloc=maxloc,OffsetVal=OffsetVal,SIMPLIFY=FALSE)
  
  RegisteredSwings=rep(list(NA),length(AllSwingsMat))
  indswingnum=lapply(AllSwingsMat,function(x){ncol(x)})
  for (k in 1:length(AllSwingsMat)){
    RegisteredSwings[[k]]=matrix(0,nrow=MatrixSize[[k]]+5000,ncol=indswingnum[[k]])
    findmax=which.max(maxloc[[k]])
    for (i in 1:indswingnum[[k]]){
      if (i != findmax){
        RegisteredSwings[[k]][Offset[[k]][i]:(Offset[[k]][i]+SizeVec[[k]][[i]]-1),i]=AllSwingsMat[[k]][1:SizeVec[[k]][[i]],i]
      } else {
        RegisteredSwings[[k]][1:SizeVec[[k]][[i]],i]=AllSwingsMat[[k]][1:SizeVec[[k]][[i]],i]
      }}}
  
  #cut out all the unneeded values after registration
  sumZeroFUN=function(Matrix,OffsetVal,cutval){
    Matrix=Matrix[-OffsetVal:-(nrow(Matrix)),]
    #take out the unneeded zeros at the start
    j=nrow(Matrix)
    Matrix=Matrix[-1:-(j-cutval),]
    return(Matrix)
  }
  RegisteredSwings=mapply(sumZeroFUN,Matrix=RegisteredSwings,OffsetVal=OffsetVal,cutval=1300,SIMPLIFY=FALSE)
  
  return(RegisteredSwings)
}


#Function to fit an average FDA curve to all registered swings
fdameanFUN=function(registeredswings,rsize,rrange,rtime,nbasis){
  registeredswings=lapply(registeredswings,function(x){x=tail(x,n=rsize)})
  #Need to redo FDA curves in order to find average
  splines=create.bspline.basis(rangeval=rrange,nbasis=nbasis,norder=6)
  #functional data object
  OBJfda=lapply(registeredswings,function(x){Data2fd(y=x,argvals=rtime,basisobj = splines)})
 
  return(OBJfda)
}

#Function to convert a list of all swings for each player to a list of FDA objects for a certain club for each player
#input list with all swings of certain variable/vector of club names/club number of choice in that vector/List of all Sizes
SortClubsfdaFUN=function(ListMAT,Size,clubs,Rsize,clubnum,val,nbasis){
  Rrange = c(1,Rsize)
  Rtime = Rrange[1]:Rrange[2]
  
  ListMAT=MoCapFUN(ListMAT,players,pswings)
  Size=MoCapFUN(Size,players,pswings)
  
  CLUB=lapply(ListMAT,function(x){list.match(x,clubs[clubnum])})
  CLUB=list.clean(CLUB,function(x) length(x)==0L,TRUE)
  Size1=lapply(Size,function(x){list.match(x,clubs[clubnum])})
  Size1=list.clean(Size1,function(Size1) length(Size1)==0L,TRUE)
  MatrixSize1=lapply(Size1,function(x){x[[which.max(x)]]})
  CLUB=CombineListMatFUN(CLUB,MatrixSize1)
  RegisteredCLUB=registerFUN(CLUB,Size1,MatrixSize1,val)
  CLUBfda=fdameanFUN(RegisteredCLUB,Rsize,Rrange,Rtime,nbasis)
  names(CLUBfda)=paste(names(CLUB)) 
  
  return(CLUBfda)
}

#set of functions to delete players that did not test both clubs of interest if applicable
ClubBigCropFUN=function(club1,club2){
  club=list(club1,club2)
  lengthvec=c(length(club1),length(club2))
  if (lengthvec[1]==lengthvec[2]){
    return(club1)
  }else{
    max=which.max(lengthvec)
    min=which.min(lengthvec)
    clubbig=club[[max]]
    clubsmall=club[[min]]
    
    NonMatch=which((names(clubbig) %!in% names(clubsmall))==TRUE)
    clubbig=list.remove(clubbig,NonMatch)
    return(clubbig)}
}

ClubSmallFUN=function(club1,club2){
  club=list(club1,club2)
  lengthvec=c(length(club1),length(club2))
  if (lengthvec[1]==lengthvec[2]){
    return(club2)
  }else{
    min=which.min(lengthvec)
    max=which.max(lengthvec)
    clubbig=club[[max]]
    clubsmall=club[[min]]
    return(clubsmall)}
}


#changed!
FDAfun=function(swing,Time,nbasis,vel){
  norder=6
  #Time lists/vectors prep for FDA --> DEFLECTION TIME VECTORS
  TimeMS=Time*1000
  Size=length(TimeMS)
  TimeRange=c(0,tail(TimeMS,n=1))
  
  #create basis function (bspline) and fit to data
  func_basis = create.bspline.basis(TimeRange,nbasis=nbasis,norder=6) 
  func_smooth = smooth.basis(TimeMS, swing, func_basis)$fd
  if (vel==0){
  Smooth_Curve=eval.fd(TimeMS,func_smooth)
  }else if (vel==1){
  Smooth_Curve=eval.fd(TimeMS,func_smooth,1)
  }else if (vel==2){
  Smooth_Curve=eval.fd(TimeMS,func_smooth,2)
  }
  
  return(Smooth_Curve)
}


#FDA on individual swings function (returns smoothed matrix velocity curve)
FDAfunVelNEW=function(swing,Time,nbasis,pos){
  timems=Time*1000
  size=length(timems)
  timerange=c(0,tail(timems,n=1))
  
  #create basis function (bspline) and fit to data
  func_basis = create.bspline.basis(timerange,nbasis=nbasis,norder=6) 
  func_smoothX = smooth.basis(timems, swing[,1], func_basis)$fd
  func_smoothY = smooth.basis(timems, swing[,2], func_basis)$fd
  func_smoothZ = smooth.basis(timems, swing[,3], func_basis)$fd
  #velocities of these curves
  Smooth_VelX=eval.fd(timems,func_smoothX,1)
  Smooth_VelY=eval.fd(timems,func_smoothY,1)
  Smooth_VelZ=eval.fd(timems,func_smoothZ,1)
  #calculate magnitude
  Smooth_Vel=sqrt(Smooth_VelX^2+Smooth_VelY^2+Smooth_VelZ^2)
  
  if (pos==1){
  return(Smooth_VelX)
  }else if (pos==2){
  return(Smooth_VelY)  
  }else if (pos==3){
  return(Smooth_VelZ)  
  }else if(pos==4){
  return(Smooth_Vel)  
  }  
}

#FDA on individual swings function (returns smoothed matrix velocity curve)
FDAfunVel=function(swing,timerange,timems,nbasis){
  norder=6
  #create basis function (bspline) and fit to data
  func_basis = create.bspline.basis(timerange,nbasis,norder=6) 
  func_smoothX = smooth.basis(timems, swing[,1], func_basis)$fd
  func_smoothY = smooth.basis(timems, swing[,2], func_basis)$fd
  func_smoothZ = smooth.basis(timems, swing[,3], func_basis)$fd
  #velocities of these curves
  Smooth_VelX=eval.fd(timems,func_smoothX,1)
  Smooth_VelY=eval.fd(timems,func_smoothY,1)
  Smooth_VelZ=eval.fd(timems,func_smoothZ,1)
  #calculate magnitude
  Smooth_Vel=sqrt(Smooth_VelX^2+Smooth_VelY^2+Smooth_VelZ^2)
  
  return(Smooth_Vel)
}


FDAfunAcc=function(swing,nbasis=120,norder=6,timerange,timems){
  
  #create basis function (bspline) and fit to data
  func_basis = create.bspline.basis(timerange,nbasis=120,norder=6) 
  func_smoothX = smooth.basis(timems, swing[,1], func_basis)$fd
  func_smoothY = smooth.basis(timems, swing[,2], func_basis)$fd
  func_smoothZ = smooth.basis(timems, swing[,3], func_basis)$fd
  #velocities of these curves
  Smooth_AccX=eval.fd(timems,func_smoothX,2)
  Smooth_AccY=eval.fd(timems,func_smoothY,2)
  Smooth_AccZ=eval.fd(timems,func_smoothZ,2)
  #calculate magnitude
  Smooth_Acc=sqrt(Smooth_AccX^2+Smooth_AccY^2+Smooth_AccZ^2)
  
  return(Smooth_Acc)
}


#Function of the entire FDA process
FDAtransFUN=function(prepmarker,Size,players,pswings,Rsize,nbasis){
  #create individual FDA curve for parameter for each swing in list
  #FDAind=mapply(FDAfun,swing=prepmarker,Time=Time,nbasis=120,0)
  
  #sort lists by player and club (parameter of interest/time/size)
  SORTEDLIST=MoCapFUN(prepmarker,players,pswings)
  #Time=MoCapFUN(Time,players,pswings)
  Size=MoCapFUN(Size,players,pswings)
  
  #combine lists of parameters under each player into 1 matrix for each player
  MatrixSize=lapply(Size,function(x){x[[which.max(x)]]})
  CombinedMAT=CombineListMatFUN(SORTEDLIST,MatrixSize)
  
  #register curves
  RegisteredSwings=registerFUN(CombinedMAT,Size,MatrixSize,0)
  
  #fit FDA curve to registered swings
  Rsize=Rsize
  Rrange=c(1,Rsize)
  Rtime=Rrange[1]:Rrange[2]
  
  #now calculate FDA for each column (needed in order to find average)
  ALLfdacurves=fdameanFUN(RegisteredSwings,Rsize,Rrange,Rtime,nbasis)
  
  return(ALLfdacurves)
}


#function to translate xyz distances to magnitude
magFUN=function(x){
  x=apply(x,1,function(matrix){sqrt(matrix[1]^2+matrix[2]^2+matrix[3]^2)})
  return(x)
}

#function to add column
add.col=function(vector,newsize){
  length(vector)=newsize
  return(vector)
}

#function to delete column
delete.col=function(matrix,newsize){
  matrix=matrix[-(newsize+1):-nrow(matrix),]
  return(matrix)
}


#Function to create Impact Position Point (Do not take derivative of artificial point!)
#finds time from last frame to impact using V_ImpactLocation/BallBackCover
#can interpolate values using this time for velocity or position

polyfitFUN=function(time,matrix,frames,order,return){
  #polyfit designated frames at input order to get parameter at impact frame
  testvec=data.frame(cbind(time,matrix))
  names(testvec)=c("x","y")
  x=testvec$x
  y=testvec$y
  PolyImpact=lm(y~ poly(x,order),data=testvec)
  PolyfitVel=predict(PolyImpact)
  if (return==0){
    return(PolyfitVel)
  } else if (return==1){
    return(PolyImpact)
  }
}

time2impact=function(swing,marker,time,frames,order,vel){
  # swing=MoCapALL[[117]]
  # time=TimeN[[117]]
  # marker=Marker1[[117]]
  # frames=15
  # order=2
  # vel=0
  # 
  time=time*1000
  rowsize=nrow(swing)
  pointsnum=rowsize-frames
  swing=swing[pointsnum:rowsize,]
  time=time[pointsnum:rowsize]
  marker=marker[pointsnum:rowsize,]
  
  #Find exact XYZ dist to impact using BallBackCover and ImpactLocation Marker
  ImpactDistX=swing[frames+1,135]-swing[frames+1,138]
  ImpactDistY=swing[frames+1,136]-swing[frames+1,139]
  ImpactDistZ=swing[frames+1,137]-swing[frames+1,140]
  ImpactDist=cbind(ImpactDistX,ImpactDistY,ImpactDistZ)
  ImpactDist=magFUN(ImpactDist)
  
  #Derive velocity (Backward Difference)
  VTot=cbind(diff(swing[,138])/diff(time),diff(swing[,139])/diff(time),diff(swing[,140])/diff(time))
  VTot=magFUN(VTot)
  
  ##polyfit designated frames at input order to get parameter at impact frame
  PolyfitVel=polyfitFUN(time[1:frames],VTot,frames,order,0)
  
  #new velocity from curve fit
  LastSeenVel=tail(PolyfitVel,n=1)
  #calculate time from last frame to impact assuming constant velocity
  #units for dt in sec*1000
  dt=ImpactDist/LastSeenVel
  
  #Now we have time from last frame to impact--> use this to find impact value
  ImpactTime=data.frame('x'=c(tail(time[frames],n=1)+dt))
  
  if (LastSeenVel>30){
  #interpolate impact velocity of marker of choice
  if (vel==1){
    #Derive velocity (Backward Difference)
    VTot2=cbind(diff(marker[,1])/diff(time),diff(marker[,2])/diff(time),diff(marker[,3])/diff(time))
    marker=magFUN(VTot2)
    markerVELeq=polyfitFUN(time[1:frames],marker,frames,order,1)
    
    #Now we have time from last frame to impact--> use this to find impact speed
    ImpactSpeed=predict(markerVELeq,ImpactTime)
    IMPACT=c(unlist(ImpactTime),unlist(ImpactSpeed))
    return(IMPACT)
  }else if (vel==0){
    markerXeq=polyfitFUN(time[1:frames],marker[,1],10,2,1)
    markerYeq=polyfitFUN(time[1:frames],marker[,2],10,2,1)
    markerZeq=polyfitFUN(time[1:frames],marker[,3],10,2,1)
    markerXPOS=predict(markerXeq,ImpactTime)
    markerYPOS=predict(markerYeq,ImpactTime)
    markerZPOS=predict(markerZeq,ImpactTime)
    intmarker=cbind(ImpactTime,markerXPOS,markerYPOS,markerZPOS)
    return(intmarker)
  }
  }
}

#confidence interval
cinterval=function(fdaswings,num){
  stdev=std.fd(fdaswings)
  sqrtN=sqrt(length(fdaswings))
  z=1.96
  
  interval=z/sqrtN*stdev
  upper=mean(fdaswings)+interval
  lower=mean(fdaswings)-interval
  if (num==0){
    return(lower)
  }else if (num==1){
    return(upper)
  }
}


#function that analyzes statistical differences through the swing for a given parameter for all players
#input--> List: list of all swings with the variable of interest/Size: list of size of matrix of all swings
#clubs: string vec of club names/clubofinterest1: integer 1,2, or 3 which represents the club number (same for clubofinterest2)
clubanalysisFUNavg=function(List,Size,clubs,clubofinterest1,clubofinterest2){
  #List=VTot
  #Size=SizeS
  # club1tag=clubs[clubofinterest1]
  # club2tag=clubs[clubofinterest2]
  nbasis=101
  
  club1=SortClubsfdaFUN(List,Size,clubs,1000,clubofinterest1,0,nbasis)
  club2=SortClubsfdaFUN(List,Size,clubs,1000,clubofinterest2,0,nbasis)
  
  #ANALYSIS of clubs of interest
  #club1 vs club2
  if (length(club1)>length(club2)){
    CLUB1=ClubBigCropFUN(club1,club2) #find club w bigger player group
    CLUB2=ClubSmallFUN(club1,club2) #find club w smaller player group
  }else{
    CLUB2=ClubBigCropFUN(club1,club2) 
    CLUB1=ClubSmallFUN(club1,club2)
  }
  AVGCLUB1=lapply(CLUB1,function(x){mean(x)$coef})
  AVGCLUB2=lapply(CLUB2,function(x){mean(x)$coef})
  
  
  statdiff=mapply(function(x,y){x-y},x=AVGCLUB1,y=AVGCLUB2,SIMPLIFY = FALSE)
  
  #calculate 95% Confidence Intervals for each club
  club1lower=mapply(cinterval,fdaswings=CLUB1,num=0,SIMPLIFY = FALSE)
  club1upper=mapply(cinterval,fdaswings=CLUB1,num=1,SIMPLIFY = FALSE)
  
  club2lower=mapply(cinterval,fdaswings=CLUB2,num=0,SIMPLIFY = FALSE)
  club2upper=mapply(cinterval,fdaswings=CLUB2,num=1,SIMPLIFY = FALSE)
  
  ttest=rep(list(NA),length(CLUB1))
  for (i in 1:length(CLUB1)){
    ttest[[i]]=tperm.fd(CLUB1[[i]],CLUB2[[i]],nperm=1000,q=.05,plotres=TRUE)
  }
  
  
  #plot ALL
  par(mfrow=c(2,2)) 
  for (i in 1:length(CLUB1)){ 
    #layout(matrix(c(1,1,2,2),2,2,byrow=TRUE))
    plot(CLUB1[[i]],col="red",main=names(CLUB1)[i])
    lines(club1lower[[i]],col="blue",lwd=4,lty=3)
    lines(club1upper[[i]],col="blue",lwd=4,lty=3)
    lines(CLUB2[[i]],col="green")
    lines(club2lower[[i]],col="darkgreen",lwd=4,lty=3)
    lines(club2upper[[i]],col="darkgreen",lwd=4,lty=3)
    plot(ttest[[i]]$pvals.pts,main="P Value")
    lines(ttest[[i]]$pvals.pts)
    abline(h=.05,col="red",lty="dashed",lwd=2)
    for (j in 1:length(statdiff[[i]])){
      if(statdiff[[i]][j]>0){color="green"}
      if(statdiff[[i]][j]<0){color="red"}
      points(j,ttest[[i]]$pvals.pts[j],col=color,lwd=4)
    }
    plot(statdiff[[i]]*2.23694,main="CLUB1 - CLUB2 DIFFERENCE")
    lines(statdiff[[i]]*2.23694)
    abline(h=0,col="black",lty="dashed")
    for (j in 1:length(statdiff[[i]])){
      if(statdiff[[i]][j]>0){color="green"}
      if(statdiff[[i]][j]<0){color="red"}
      points(j,statdiff[[i]][j]*2.23694,col=color,lwd=4)
    }
    plot(CLUB1[[i]]*2.23694,xlim=c(970,1000),col="red")
    legend('bottomright',ncol=3L,title='Club Speed',text.width = 2.5,
           legend=c('Club1 AVG',signif(AVGCLUB1[[i]][nbasis]*2.23694,5),
                    'Club2 AVG',signif(AVGCLUB2[[i]][nbasis]*2.23694,5),
                    'Club1-Club2',signif(AVGCLUB1[[i]][nbasis]*2.23694-AVGCLUB2[[i]][nbasis]*2.23694,5)))
    lines(mean(CLUB1[[i]])*2.23694,col="red",lwd=4)
    lines(CLUB2[[i]]*2.23694,col="blue")
    lines(mean(CLUB2[[i]])*2.23694,col="blue",lwd=4)}
  
  
  
  return(ttest)
}


#takes the manually clicked scrubbed values "ScrubValList" and deletes them from "List"
ManualReplaceFUN=function(List,ScrubValList){
  if (is.na(ScrubValList[1])==FALSE){
    List=List[-ScrubValList]
  }else{
    List=List
  }
  return(List)
}


SimpleScrubFUN=function(V,TimeN,num){
  
  #easily delete the really bad outliers first
  if (V[which.max(V)]*2.23694<50){
    chopnum=2
    chopnum2=10
  }else{
    chopnum=4
    chopnum2=4
  }
  
  last=V[length(V)]
  last=last+chopnum2
  bad=which(V>last)
  if (length(bad)>0){
    V=V[-bad]
    TimeN=TimeN[-bad]}
  
  if (num==1){return(V)}
  if (num==2){return(TimeN)}
}


AccelScrubFUN=function(ATot,V,TimeN,num){
  
  quantiles=quantile(ATot,probs=c(.015,.985))
  
  x=ATot>quantiles[1]
  x2=which(x==FALSE)
  y=ATot<quantiles[2]
  y2=which(y==FALSE)
  bad=c(x2,y2)
  
  if (num==1){
    V=V[-bad]
    return(V)
  }
  if (num==2){
    TimeN=TimeN[-bad]
    return(TimeN)
  }
}


ManualScrubFUN=function(ATot,VTot){
  #ATot=ATot[[122]]
  #VTot=VTot[[122]]
  
  maxaccel=ATot[which.max(ATot)]
  
  if (maxaccel>1.5){
    plot(VTot)
    bad=identify(VTot)
    return(bad)
  }else{
    return(NA)
  }
}



#--------------------- Load Data and configure for modeling -----------------------------------
#load all MoCap swings into a list
temp=list.files(pattern="*.trc")
MoCapALL=lapply(temp,freadFUN)


#find bad ball hat capture swings
# badhatlogical=lapply(MoCapALL,findbadhatFUN)
# badhats=which(badhatlogical==FALSE)
# #delete these swings from MoCapALL
# MoCapALL=MoCapALL[-badhats]
# temp=temp[-badhats]


#function returns logical list showing which swings have messed up impact location
#bad swings have ball hat coordinates as the last marker
#if I subtract the first 2 rows of the last marker and get zero, I know it's a bad swing
badswingFUN=function(listX){
  col1=unlist(listX[1:2,138])
  if (typeof(col1[1])=="character" & typeof(col1[2])=="character"){
    difflogical=TRUE
  } else{
  diffcol1=col1[1]-col1[2]
  difflogical=(diffcol1==0 || is.na(diffcol1==TRUE))
  }
  
  return(difflogical)
}

#find swings that have messed up markers (No impact marker availability)
badswinglogical=lapply(MoCapALL,badswingFUN)

badswings=which(badswinglogical==TRUE)
#delete these swings from MoCapALL
#*************************************find better solution later! aka try to salvage the swings
MoCapALL=MoCapALL[-badswings]
temp=temp[-badswings]

#delete other bad swings (fix filtering but later!)
#badswings=c(117,161,163,169,177,178,179,200,219)
#MoCapALL[badswings]=NULL
#temp=temp[-badswings]

#delete huge gaps after impact occurs so that it can be accurately spline fit
cutpoints=lapply(MoCapALL,impactFUN)


# cutpoints=rep(list(NA),length(MoCapALL))
# for (i in 1:length(MoCapALL)){
# cutpoints[[i]]=impactFUN(MoCapALL[[i]])
# }


# #function to return the impact location of each swing (last frame seen prior to impact)
# #currently set to last frame before impact
# impactFUN=function(listX){
#   listX=MoCapALL[[36]]
#   ballhat=unlist(listX[1,135])
#   V_ImpactLocation=listX[,138]
#   #val=V_ImpactLocation-ballhat
#   
#   maxloc=which.max(matrix(unlist(V_ImpactLocation),nrow = nrow(V_ImpactLocation), ncol = length(V_ImpactLocation)))
#   
#   #new addition 10/29
#   if (length(which(is.na(V_ImpactLocation[(maxloc-30):maxloc])))>25){
#     listX=DeletePostImpactFUN(listX,maxloc-1)
#     V_ImpactLocation=listX[,138]
#     maxloc=which.max(matrix(unlist(V_ImpactLocation),nrow = nrow(V_ImpactLocation), ncol = length(V_ImpactLocation)))
#   }
#   
#   if (length(which(is.na(V_ImpactLocation[maxloc:(maxloc+50)])))<20){
#     listX=DeletePostImpactFUN(listX,maxloc-1)
#     V_ImpactLocation=listX[,138]
#     maxloc=which.max(matrix(unlist(V_ImpactLocation),nrow = nrow(V_ImpactLocation), ncol = length(V_ImpactLocation)))
#   }
#   #end of new addition
#   
#   diffvals=V_ImpactLocation-ballhat
#   
#   j=maxloc
#   while (diffvals[j]>=0 || is.na(diffvals[j])==TRUE){
#     j=j-1
#   }
#   
#   #0 frames past impact
#   j=j+0
#   
#   return(j)
# }



MoCapALL=mapply(DeletePostImpactFUN,MoCapALL=MoCapALL,cutpointsX=cutpoints,SIMPLIFY=FALSE)
#Also delete NA column #141
MoCapALL=lapply(MoCapALL,function(list){if (ncol(list)>140){list[,-141]}else{list=list}})

# #delete extra points at the beginning for long swings (2000 pt threshold)
# MoCapALL=lapply(MoCapALL,function(list){
#   if (nrow(list)>1998){
#   offset=nrow(list)-1998
#   list=list[offset:nrow(list)]
#   }else{list=list}
#   })

#spline interpolate missing values
MoCapALL=lapply(MoCapALL,function(list){list=apply(list,2,na.spline)})

#------------------------ Function to Sort all imported swings by player ------------------------------------

#extract player strings
temp2=str_sub(temp,8,13)
#extract club strings
temp3=str_sub(temp,18,26)
temp3[1:length(temp3)]="WedgeFull"

#create string of players that tested
players=unique(temp2)
#create string of clubs that were tested
clubs=unique(temp3)

#find number of swings of each player
pswings=match(players,temp2)

#Assign club names to swings first before organizing by player
names(MoCapALL)=paste(temp3)

#create list of lists for each player and their swings
size=length(MoCapALL)

#Use this function in each analysis to sort by player AT THE END of calcs
MoCapFUN=function(MoCapALL,players,pswings){
  #List of lists (each list is a player/lists in each list represent swings)
  MoCapALL2=rep(list(NA),length(players))
  indswingnum=rep(NA,length(players))
  for (i in 1:length(players)){
    if (i!=length(players)){
      MoCapALL2[[i]]=MoCapALL[c(pswings[i]:(pswings[i+1]-1))]
      indswingnum[i]=length(MoCapALL2[[i]])
    }else{
      MoCapALL2[[i]]=MoCapALL[c(pswings[i]:size)]
      indswingnum[i]=length(MoCapALL2[[i]])
    }
  }
  names(MoCapALL2)=paste(players)
  return(MoCapALL2)
}


#---------------------- Sort by markers -----------------------------------------------

#Adjust Time Vec bc it is rounded incorrectly in imported file
Time=lapply(MoCapALL,function(x){x[,2]})
Time=lapply(Time,function(x){seq(0,(length(x)-1)*.00133333333333333,.00133333333333333)})

#Import All Markers as Lists
Hosel1=lapply(MoCapALL,function(x){x[,3:5]})
Hosel2=lapply(MoCapALL,function(x){x[,6:8]})
UpperToe=lapply(MoCapALL,function(x){x[,9:11]})
GripTop=lapply(MoCapALL,function(x){x[,12:14]})
S1B=lapply(MoCapALL,function(x){x[,15:17]})
S1F=lapply(MoCapALL,function(x){x[,18:20]})
V_LFrameOrigGR=lapply(MoCapALL,function(x){x[,21:23]})
V_LFrameFaceGR=lapply(MoCapALL,function(x){x[,24:26]})
V_LFrameLoftGR=lapply(MoCapALL,function(x){x[,27:29]})
V_LFrameOrig=lapply(MoCapALL,function(x){x[,30:32]})
V_LFrameFace=lapply(MoCapALL,function(x){x[,33:35]})
V_LFrameLoft=lapply(MoCapALL,function(x){x[,36:38]})
V_FaceNormal=lapply(MoCapALL,function(x){x[,39:41]})
V_FaceNormalGR=lapply(MoCapALL,function(x){x[,42:44]})
V_CenterFM20HD=lapply(MoCapALL,function(x){x[,45:47]})
V_LFrameFM20GR=lapply(MoCapALL,function(x){x[,48:50]})
V_CenterFace=lapply(MoCapALL,function(x){x[,51:53]})
V_ImpactLocLE=lapply(MoCapALL,function(x){x[,54:56]})
V_ImpactLocToe=lapply(MoCapALL,function(x){x[,57:59]})
V_CenterFaceGR=lapply(MoCapALL,function(x){x[,60:62]})
V_ShaftCenter=lapply(MoCapALL,function(x){x[,63:65]})
V_HoselCenterGR=lapply(MoCapALL,function(x){x[,66:68]})
V_UpShaftGR=lapply(MoCapALL,function(x){x[,69:71]})
V_FaceParallelGR=lapply(MoCapALL,function(x){x[,72:74]})
V_GripYAxis=lapply(MoCapALL,function(x){x[,75:77]})
V_GripXAxis=lapply(MoCapALL,function(x){x[,78:80]})
V_GripZAxis=lapply(MoCapALL,function(x){x[,81:83]})
V_HoselCenter=lapply(MoCapALL,function(x){x[,84:86]})
V_UpShaft=lapply(MoCapALL,function(x){x[,87:89]})
V_FaceParallel=lapply(MoCapALL,function(x){x[,90:92]})
V_HeadYAxis=lapply(MoCapALL,function(x){x[,93:95]})
V_HeadXAxis=lapply(MoCapALL,function(x){x[,96:98]})
V_HeadZAxis=lapply(MoCapALL,function(x){x[,99:101]})
V_HoselCenterALT=lapply(MoCapALL,function(x){x[,102:104]})
V_HeadXAxisALT=lapply(MoCapALL,function(x){x[,105:107]})
V_HeadYAxisALT=lapply(MoCapALL,function(x){x[,108:110]})
V_HeadZAxisALT=lapply(MoCapALL,function(x){x[,111:113]})
V_LFrameLoftALT=lapply(MoCapALL,function(x){x[,114:116]})
V_LFrameFaceALT=lapply(MoCapALL,function(x){x[,117:119]})
V_LFrameFaceALTALT=lapply(MoCapALL,function(x){x[,120:122]})
V_FaceNormalALT=lapply(MoCapALL,function(x){x[,123:125]})
V_HCequal=lapply(MoCapALL,function(x){x[,126:128]})
V_FCequal=lapply(MoCapALL,function(x){x[,129:131]})
V_Center=lapply(MoCapALL,function(x){x[,132:134]})
V_BallBackCover=lapply(MoCapALL,function(x){x[,135:137]})
V_ImpactLocation=lapply(MoCapALL,function(x){x[,138:140]})


#----------------------------- ANALYSIS ----------------------------------




#***********************ANGULAR VELOCITIES*****************************
ShaftCS_origin=lapply(V_HoselCenterGR,function(x){x/1000})
HeadCS_origin=lapply(V_HoselCenter,function(x){x/1000})
Greg_GripX=lapply(V_GripXAxis,function(x){x/1000})
Greg_GripY=lapply(V_GripYAxis,function(x){x/1000})
Greg_GripZ=lapply(V_GripZAxis,function(x){x/1000})
Greg_HeadX=lapply(V_HeadXAxis,function(x){x/1000})
Greg_HeadY=lapply(V_HeadYAxis,function(x){x/1000})
Greg_HeadZ=lapply(V_HeadZAxis,function(x){x/1000})
TimeVec=lapply(Time,function(x){x/1000}) #fdaFUN will multiply this back to normal later
Size=lapply(Time,function(x){length(x)})

#correct shaft/head axis VM's to be RH (Align with Global CS)

#function to find the norm
norm_vec=function(x){sqrt(sum(x^2))}
#Calculate direction vectors for Greg's CS axis
Greg_Shaft_unitY=mapply(function(Greg_GripX,ShaftCS_origin){
  sweep(sweep(Greg_GripX,1,ShaftCS_origin),1,apply(sweep(Greg_GripX,1,ShaftCS_origin),1,norm_vec),"/")},
  Greg_GripX=Greg_GripX,ShaftCS_origin=ShaftCS_origin)
Greg_Shaft_unitX=mapply(function(Greg_GripY,ShaftCS_origin){
  sweep(sweep(Greg_GripY,1,ShaftCS_origin),1,apply(sweep(Greg_GripY,1,ShaftCS_origin),1,norm_vec),"/")},
  Greg_GripY=Greg_GripY,ShaftCS_origin=ShaftCS_origin)#Y axis should be X
Greg_Shaft_unitZ=mapply(function(Greg_GripZ,ShaftCS_origin){
  sweep(sweep(Greg_GripZ,1,ShaftCS_origin),1,apply(sweep(Greg_GripZ,1,ShaftCS_origin),1,norm_vec),"/")},
  Greg_GripZ=Greg_GripZ,ShaftCS_origin=ShaftCS_origin)

Greg_Head_unitY=mapply(function(Greg_HeadX,HeadCS_origin){
  sweep(sweep(Greg_HeadX,1,HeadCS_origin),1,apply(sweep(Greg_HeadX,1,HeadCS_origin),1,norm_vec),"/")},
  Greg_HeadX=Greg_HeadX,HeadCS_origin=HeadCS_origin)
Greg_Head_unitX=mapply(function(Greg_HeadY,HeadCS_origin){
  sweep(sweep(Greg_HeadY,1,HeadCS_origin),1,apply(sweep(Greg_HeadY,1,HeadCS_origin),1,norm_vec),"/")},
  Greg_HeadY=Greg_HeadY,HeadCS_origin=HeadCS_origin)
Greg_Head_unitZ=mapply(function(Greg_HeadZ,HeadCS_origin){
  sweep(sweep(Greg_HeadZ,1,HeadCS_origin),1,apply(sweep(Greg_HeadZ,1,HeadCS_origin),1,norm_vec),"/")},
  Greg_HeadZ=Greg_HeadZ,HeadCS_origin=HeadCS_origin)


#Transformation matrices for Greg's coordinate systems
SizeX=lapply(Greg_Shaft_unitY,function(x){nrow(x)})
R_Greg_Shaft=rep(list(NA),length(SizeX))
R_Greg_Head=rep(list(NA),length(SizeX))
for (j in 1:length(SizeX)){
  R_Greg_Shaft[[j]]=rep(list(NA),SizeX[[j]])
  R_Greg_Head[[j]]=rep(list(NA),SizeX[[j]])
  for (i in 1:SizeX[[j]]){
  R_Greg_Shaft[[j]][[i]]=rbind(Greg_Shaft_unitX[[j]][i,],Greg_Shaft_unitY[[j]][i,],Greg_Shaft_unitZ[[j]][i,])
  R_Greg_Head[[j]][[i]]=rbind(Greg_Head_unitX[[j]][i,],Greg_Head_unitY[[j]][i,],Greg_Head_unitZ[[j]][i,])
  }
}


#Create new VMs for corrected coordinate system (aligned with global)
GripXAxis_Correct=rep(list(NA),length(SizeX))
GripYAxis_Correct=rep(list(NA),length(SizeX))
GripZAxis_Correct=rep(list(NA),length(SizeX))
HeadXAxis_Correct=rep(list(NA),length(SizeX))
HeadYAxis_Correct=rep(list(NA),length(SizeX))
HeadZAxis_Correct=rep(list(NA),length(SizeX))

for (j in 1:length(SizeX)){
GripXAxis_Correct[[j]]=matrix(NA,nrow=SizeX[[j]],ncol=3)
GripYAxis_Correct[[j]]=matrix(NA,nrow=SizeX[[j]],ncol=3)
GripZAxis_Correct[[j]]=matrix(NA,nrow=SizeX[[j]],ncol=3)

HeadXAxis_Correct[[j]]=matrix(NA,nrow=SizeX[[j]],ncol=3)
HeadYAxis_Correct[[j]]=matrix(NA,nrow=SizeX[[j]],ncol=3)
HeadZAxis_Correct[[j]]=matrix(NA,nrow=SizeX[[j]],ncol=3)

for (i in 1:length(R_Greg_Shaft[[j]])){
  GripXAxis_Correct[[j]][i,]=ShaftCS_origin[[j]][i,]+t(t(R_Greg_Shaft[[j]][[i]])%*%rbind(-.3,0,0))
  GripYAxis_Correct[[j]][i,]=ShaftCS_origin[[j]][i,]+t(t(R_Greg_Shaft[[j]][[i]])%*%rbind(0,-.3,0))
  GripZAxis_Correct[[j]][i,]=Greg_GripZ[[j]][i,]
  
  HeadXAxis_Correct[[j]][i,]=HeadCS_origin[[j]][i,]+t(t(R_Greg_Head[[j]][[i]])%*%rbind(-.3,0,0))
  HeadYAxis_Correct[[j]][i,]=HeadCS_origin[[j]][i,]+t(t(R_Greg_Head[[j]][[i]])%*%rbind(0,-.3,0))
  HeadZAxis_Correct[[j]][i,]=Greg_HeadZ[[j]][i,]
}
}

#calculate direction vectors for correct cs axis
Correct_Shaft_unitX=mapply(function(GripXAxis_Correct,ShaftCS_origin){
  sweep(sweep(GripXAxis_Correct,1,ShaftCS_origin),1,apply(sweep(GripXAxis_Correct,1,ShaftCS_origin),1,norm_vec),"/")},
  GripXAxis_Correct=GripXAxis_Correct,ShaftCS_origin=ShaftCS_origin)
Correct_Shaft_unitY=mapply(function(GripYAxis_Correct,ShaftCS_origin){
  sweep(sweep(GripYAxis_Correct,1,ShaftCS_origin),1,apply(sweep(GripYAxis_Correct,1,ShaftCS_origin),1,norm_vec),"/")},
  GripYAxis_Correct=GripYAxis_Correct,ShaftCS_origin=ShaftCS_origin)
Correct_Shaft_unitZ=mapply(function(GripZAxis_Correct,ShaftCS_origin){
  sweep(sweep(GripZAxis_Correct,1,ShaftCS_origin),1,apply(sweep(GripZAxis_Correct,1,ShaftCS_origin),1,norm_vec),"/")},
  GripZAxis_Correct=GripZAxis_Correct,ShaftCS_origin=ShaftCS_origin)

Correct_Head_unitX=mapply(function(HeadXAxis_Correct,HeadCS_origin){
  sweep(sweep(HeadXAxis_Correct,1,HeadCS_origin),1,apply(sweep(HeadXAxis_Correct,1,HeadCS_origin),1,norm_vec),"/")},
  HeadXAxis_Correct=HeadXAxis_Correct,HeadCS_origin=HeadCS_origin)
Correct_Head_unitY=mapply(function(HeadYAxis_Correct,HeadCS_origin){
  sweep(sweep(HeadYAxis_Correct,1,HeadCS_origin),1,apply(sweep(HeadYAxis_Correct,1,HeadCS_origin),1,norm_vec),"/")},
  HeadYAxis_Correct=HeadYAxis_Correct,HeadCS_origin=HeadCS_origin)
Correct_Head_unitZ=mapply(function(HeadZAxis_Correct,HeadCS_origin){
  sweep(sweep(HeadZAxis_Correct,1,HeadCS_origin),1,apply(sweep(HeadZAxis_Correct,1,HeadCS_origin),1,norm_vec),"/")},
  HeadZAxis_Correct=HeadZAxis_Correct,HeadCS_origin=HeadCS_origin)

#transformation matrices for corrected coordinate systems
R_Correct_Shaft=rep(list(NA),length(SizeX))
R_Correct_Head=rep(list(NA),length(SizeX))

for (j in 1:length(SizeX)){
  R_Correct_Shaft[[j]]=rep(list(NA),SizeX[[j]])
  R_Correct_Head[[j]]=rep(list(NA),SizeX[[j]])
for (i in 1:SizeX[[j]]){
  R_Correct_Shaft[[j]][[i]]=rbind(Correct_Shaft_unitX[[j]][i,],Correct_Shaft_unitY[[j]][i,],Correct_Shaft_unitZ[[j]][i,])
  R_Correct_Head[[j]][[i]]=rbind(Correct_Head_unitX[[j]][i,],Correct_Head_unitY[[j]][i,],Correct_Head_unitZ[[j]][i,])
}}

#CALCULATE EULER ANGLES ASSOCIATED WITH A GIVEN [RMAT] ----
Rmat=R_Correct_Shaft
psi=rep(list(NA),length(R_Correct_Head))
th=rep(list(NA),length(R_Correct_Head))
phi=rep(list(NA),length(R_Correct_Head))

for (j in 1:length(R_Correct_Head)){
  psi[[j]]=rep(NA,length(R_Correct_Head[[j]]))
  th[[j]]=rep(NA,length(R_Correct_Head[[j]]))
  phi[[j]]=rep(NA,length(R_Correct_Head[[j]]))
for (i in 1:length(R_Correct_Head[[j]])){
  #ZYZ Body-Fixed (intrinsic) EulerAngles (phi,th,psi) -> Rz(phi)Ry(th)Rz(psi)
  sy=sqrt(Rmat[[j]][[i]][3,1]^2+Rmat[[j]][[i]][3,2]^2)
  if (sy<.Machine$double.neg.eps){
    psi[[j]][i]=atan2(-Rmat[[j]][[i]][3,2],Rmat[[j]][[i]][3,1])#rad
    th[[j]][i]=atan2(sy,Rmat[[j]][[i]][3,3])#rad
    phi[[j]][i]=0#rad
  }else{
    psi[[j]][i]=atan2(Rmat[[j]][[i]][3,2],Rmat[[j]][[i]][3,1])
    th[[j]][i]=atan2(sy,Rmat[[j]][[i]][3,3])
    phi[[j]][i]=atan2(Rmat[[j]][[i]][2,3],-Rmat[[j]][[i]][1,3])
  }
}}

#correct jumps between -pi and pi
library(signal)
psi=lapply(psi,function(x){unwrap(x,tol=2*pi/3)})
th=lapply(th,function(x){unwrap(x,tol=2*pi/3)})
phi=lapply(phi,function(x){unwrap(x,tol=2*pi/3)})

#estimate the euler rates
psi_smooth=rep(list(NA),length(psi))
th_smooth=rep(list(NA),length(th))
phi_smooth=rep(list(NA),length(phi))
psi_dot=rep(list(NA),length(psi))
th_dot=rep(list(NA),length(th))
phi_dot=rep(list(NA),length(phi))
w_Euler_body=rep(list(NA),length(th))
w_Euler_world=rep(list(NA),length(phi))
ShaftRotation=rep(list(NA),length(phi))
FaceClosure=rep(list(NA),length(phi))

for (i in 1:length(psi)){
psi_smooth[[i]]=FDAfun(psi[[i]],TimeVec[[i]],300,0)
th_smooth[[i]]=FDAfun(th[[i]],TimeVec[[i]],300,0)
phi_smooth[[i]]=FDAfun(phi[[i]],TimeVec[[i]],300,0)

psi_dot[[i]]=FDAfun(psi[[i]],TimeVec[[i]],300,1)
th_dot[[i]]=FDAfun(th[[i]],TimeVec[[i]],300,1)
phi_dot[[i]]=FDAfun(phi[[i]],TimeVec[[i]],300,1)

#calculate the angular velocity from euler angles/rates (body-frame) ---
w_Euler_body[[i]]=matrix(NA,nrow=length(psi[[i]]),ncol=3)
w_Euler_world[[i]]=matrix(NA,nrow=length(psi[[i]]),ncol=3)
for (j in 1:length(psi[[i]])){
  w_Euler_body[[i]][j,] = c(-psi_dot[[i]][j]*sin(th[[i]][j])*cos(phi[[i]][j]) + th_dot[[i]][j]*sin(phi[[i]][j]),
                       psi_dot[[i]][j]*sin(th[[i]][j])*sin(phi[[i]][j]) + th_dot[[i]][j]*cos(phi[[i]][j]),
                       psi_dot[[i]][j]*cos(th[[i]][j]) + phi_dot[[i]][j])
  
  x=data.frame(w_Euler_body[[i]][j,])
  x=data.matrix(x)
  w_Euler_world[[i]][j,] = t(t(Rmat[[i]][[j]])%*%x)
}

ShaftRotation[[i]]=w_Euler_body[[i]][,3]
FaceClosure[[i]]=w_Euler_world[[i]][,3]
}

#####
#Interpolate Impact 
#####

#sort variables by player/club
ShaftRotationALL=MoCapFUN(ShaftRotation,players,pswings)
FaceClosureALL=MoCapFUN(FaceClosure,players,pswings)
Size=MoCapFUN(Size,players,pswings)
TimeV=MoCapFUN(TimeVec,players,pswings)

#combine all into 1 matrix for each player
MatrixSize=lapply(Size,function(x){x[[which.max(x)]]})
ShaftRotationALL2=CombineListMatFUN(ShaftRotationALL,MatrixSize)
FaceClosureALL2=CombineListMatFUN(FaceClosureALL,MatrixSize)

#Register Curves
#sorted by player only
RegisteredShaftRot=registerFUN(ShaftRotationALL2,Size,MatrixSize,0)
RegisteredFaceClos=registerFUN(FaceClosureALL2,Size,MatrixSize,0)

#fit FDA curve to registered Swings
Rsize=1000
Rrange=c(1,1000)
Rtime=Rrange[1]:Rrange[2]

ShaftRotFDA=fdameanFUN(RegisteredShaftRot,Rsize,Rrange,Rtime,120)
FaceClosFDA=fdameanFUN(RegisteredFaceClos,Rsize,Rrange,Rtime,120)

#Avg of all smoothed swings
AVGShRot=lapply(ShaftRotFDA,function(x){mean(x)})
AVGFaClo=lapply(FaceClosFDA,function(x){mean(x)})

#plot gut check
plot.fd(ShaftRotFDA[[1]],col="red",lwd=1,xlim=c(0,1000),ylim=c(-5,15))
lines(AVGShRot[[1]],col="blue",lwd=2)



#All players
#ALL PLAYERS
par(mfrow=c(4,4))
for (i in 1:length(players)){
  plot.fd(ShaftRotFDA[[i]],col="coral2",xlim=c(200,1000),ylim=c(-5,15),main=players[i])
  lines(AVGShRot[[i]],col="red",lwd=4)
  lines(FaceClosFDA[[i]],col="cadetblue1")
  lines(AVGFaClo[[i]],col="blue",lwd=4)
}



#ANGULAR VELOCITY ABOUT A CERTAIN POINT
POTranslation=GripTop
#let's say the force is input 75 mm from grip top down the shaft
#POTranslation=lapply(POTranslation,function(x){x[,3]-75)})

POInterest=V_CenterFaceGR
ANGvel=w_Euler_world


R=mapply(function(POInterest,POTranslation){sweep(POInterest,1,POTranslation)},POInterest=POInterest,POTranslation=POTranslation)
#convert to meters
R=lapply(R,function(R){sweep(R,1,1000,"/")})


#find rotational velocity of POInterest
Vrot=rep(list(NA),length(R))
for (j in 1:length(R)){
Vrot[[j]]=matrix(NA,nrow=nrow(R[[j]]),ncol=ncol(R[[j]]))
for (i in 1:nrow(R[[j]])){
#V=omega x radius (cross product)
Vrot[[j]][i,]=c(ANGvel[[j]][i,2]*R[[j]][i,3]-ANGvel[[j]][i,3]*R[[j]][i,2],
                ANGvel[[j]][i,3]*R[[j]][i,1]-ANGvel[[j]][i,1]*R[[j]][i,3],
                ANGvel[[j]][i,1]*R[[j]][i,2]-ANGvel[[j]][i,2]*R[[j]][i,1])
}}


plot(Vrot[[1]][,1],ylim=c(-40,55))
lines(Vrot[[1]][,2],col="red",lwd=3)
lines(Vrot[[1]][,3],col="blue",lwd=3)


# VrotMAG=magFUN(Vrot)
# VrotMAGmph=VrotMAG*2.2369419
# plot(VrotMAGmph)

plot(w_Euler_body[[1]][,1],ylim=c(-60,80))
lines(w_Euler_body[[1]][,2],col="red",lwd=3)
lines(w_Euler_body[[1]][,3],col="blue",lwd=3)


#speed of the hands = speed of the pivot point
POTranslation=lapply(POTranslation,function(x){x/1000})
handspeedX=mapply(FDAfunVelNEW,swing=POTranslation,Time=TimeVec,nbasis=100,pos=1)
handspeedY=mapply(FDAfunVelNEW,swing=POTranslation,Time=TimeVec,nbasis=100,pos=2)
handspeedZ=mapply(FDAfunVelNEW,swing=POTranslation,Time=TimeVec,nbasis=100,pos=3)
handspeed=mapply(function(x,y,z){cbind(x,y,z)},x=handspeedX,y=handspeedY,z=handspeedZ)

#individual magnitudes
handspeedMAG=lapply(handspeed,function(x){magFUN(x)*2.2369419})
VrotMAG=lapply(Vrot,function(x){magFUN(x)*2.2369419})

Est.Speed.Tot=mapply(function(x,y){x+y},x=handspeed,y=Vrot)
Est.Speed.TotMAG=lapply(Est.Speed.Tot,function(x){magFUN(x)*2.2369419})


#Speed breakdown
plot(VrotMAG[[1]])
lines(VrotMAG[[1]])
lines(handspeedMAG[[1]])
lines(ShaftSpeedMAG[[1]],col="red")
lines(Velfda[[1]]*2.2369419)

#xyz breakdown
plot(Velfda[[1]]*2.2369419,ylim=c(-75,122),xlim=c(1000,1900))
lines(Vrot[[1]][,1]*2.2369419)
lines(Vrot[[1]][,2]*2.2369419,col="red",lwd=3)
lines(Vrot[[1]][,3]*2.2369419,col="blue",lwd=3)
lines(handspeedX[[1]]*2.2369419,col="azure4",lwd=3)
lines(handspeedY[[1]]*2.2369419,col="coral",lwd=3)
lines(handspeedZ[[1]]*2.2369419,col="cadetblue1",lwd=3)




#*****************SWING DIRECTION*********************
marker1=V_HoselCenter

#--------------------------------MANUAL CALCULATIONS-------------------------------

#get your swing direction at each point
Direction=lapply(marker1,function(marker1){diff(marker1,lag=2)})
#convert to unit vector
UnitVec=lapply(Direction,function(Direction){apply(Direction,1,function(x){sqrt(abs(x[1])^2+abs(x[2])^2+abs(x[3])^2)})})

#final unit vector for SWING DIRECTION VARIABLE
FinalUV=rep(list(NA),length(UnitVec))
for (i in 1:length(UnitVec)){
  FinalUV[[i]]=matrix(NA,nrow=length(UnitVec[[i]]),ncol=ncol(Direction[[i]]))
  for (j in 1:length(UnitVec[[i]])){
   FinalUV[[i]][j,]=Direction[[i]][j,]/UnitVec[[i]][j]
  }}

#update size for Time vector
TimeD=mapply(function(Time,FinalUV){tail(Time,n=nrow(FinalUV))},Time=Time,FinalUV=FinalUV)
SizeD=lapply(TimeD,function(x){length(x)})

#divide up into each direction
SwingXdir=lapply(FinalUV,function(x){x[,1]})
SwingYdir=lapply(FinalUV,function(x){x[,2]})
SwingZdir=lapply(FinalUV,function(x){x[,3]})


#-------------FDA------------
SwingXdir_smoothed=FDAtransFUN(SwingXdir,SizeD,players,pswings,1000,120)
SwingYdir_smoothed=FDAtransFUN(SwingYdir,SizeD,players,pswings,1000,120)
SwingZdir_smoothed=FDAtransFUN(SwingZdir,SizeD,players,pswings,1000,120)

AVGSwingXdir=lapply(SwingXdir_smoothed,function(x){mean(x)})
AVGSwingYdir=lapply(SwingYdir_smoothed,function(x){mean(x)})
AVGSwingZdir=lapply(SwingZdir_smoothed,function(x){mean(x)})

#plot gut check for 1 direction
plot(SwingZdir_smoothed[[1]],col="red")
lines(AVGSwingZdir[[1]],col="blue",lwd=3)

#plot all directions per player
plot.fd(SwingXdir_smoothed[[1]],col="coral2")
lines(AVGSwingXdir[[1]],col="red",lwd=4)
lines(SwingYdir_smoothed[[1]],col="cadetblue1")
lines(AVGSwingYdir[[1]],col="blue",lwd=4)
lines(SwingZdir_smoothed[[1]],col="darkolivegreen1")
lines(AVGSwingZdir[[1]],col="darkgreen",lwd=4)

#ALL PLAYERS
par(mfrow=c(4,4))
for (i in 1:length(players)){
  plot.fd(SwingXdir_smoothed[[i]],col="coral2",main=players[i])
  lines(AVGSwingXdir[[i]],col="red",lwd=4)
  lines(SwingYdir_smoothed[[i]],col="cadetblue1")
  lines(AVGSwingYdir[[i]],col="blue",lwd=4)
  lines(SwingZdir_smoothed[[i]],col="darkolivegreen1")
  lines(AVGSwingZdir[[i]],col="darkgreen",lwd=4)
}



#*****************MOI ENGAGEMENT*********************

#---------------------------------MANUAL CALCULATIONS-----------------------------------
#Continuation of Swing Direction code to find MOI XYZ contribution during swing
#This is a measure of face angle relative to swing direction to find which MOI axis is engaged
marker2=V_CenterFace
marker3=V_FaceNormal

#adjust size by deleting 1st and last
marker2=lapply(marker2,function(x){
  x=head(x,n=-1)
  x=tail(x,n=-1)
})
marker3=lapply(marker3,function(x){
  x=head(x,n=-1)
  x=tail(x,n=-1)
})

#get vector of CenterFace and FaceNormal
FaceAngle=mapply(function(marker2,marker3){sweep(marker3,1,marker2)},marker2=marker2,marker3=marker3)
#convert to unit vector
FAngleMAG=lapply(FaceAngle,function(y){apply(y,1,function(x){sqrt(abs(x[1])^2+abs(x[2])^2+abs(x[3])^2)})})
#final unit vector for direction
FAUnitVec=mapply(function(FaceAngle,FAngleMAG){sweep(FaceAngle,1,FAngleMAG,"/")},FaceAngle=FaceAngle,FAngleMAG=FAngleMAG)

#plot face angle direction
plot(FAUnitVec[[1]][,1],col="blue")
lines(FAUnitVec[[1]][,2],col="red")
lines(FAUnitVec[[1]][,3],col="green")

#dot product
dotMOIX=rep(list(NA),length(FAUnitVec))
for (j in 1:length(FAUnitVec)){
  dotMOIX[[j]]=matrix(NA,nrow=nrow(FAUnitVec[[j]]),ncol=1)
  for (i in 1:nrow(FAUnitVec[[j]])){
    dotMOIX[[j]][i]=sum(FAUnitVec[[j]][i,]*FinalUV[[j]][i,])
  }}

#take absolute value (swing directions negative)
dotMOIXa=lapply(dotMOIX,function(dotMOIX){abs(dotMOIX)})
dotMOIZ=lapply(dotMOIXa,function(x){1-x})

#gut check
plot(dotMOIXa[[1]],col="red")
lines(dotMOIZ[[1]],col="blue")

#-------FDA--------
dotMOIXa_smoothed=FDAtransFUN(dotMOIXa,SizeD,players,pswings,1000,120)
dotMOIZ_smoothed=FDAtransFUN(dotMOIZ,SizeD,players,pswings,1000,120)

AVGdotMOIXa=lapply(dotMOIXa_smoothed,function(x){mean(x)})
AVGdotMOIZ=lapply(dotMOIZ_smoothed,function(x){mean(x)})

#plot all directions per player
plot.fd(dotMOIXa_smoothed[[1]],col="coral2",xlim=c(180,1000))
lines(AVGdotMOIXa[[1]],col="red",lwd=4)
lines(dotMOIZ_smoothed[[1]],col="cadetblue1")
lines(AVGdotMOIZ[[1]],col="blue",lwd=4)


#ALL PLAYERS
par(mfrow=c(4,4))
for (i in 1:length(players)){
  plot.fd(dotMOIXa_smoothed[[i]],col="coral2",main=players[i])
  lines(AVGdotMOIXa[[i]],col="red",lwd=4)
  lines(dotMOIZ_smoothed[[i]],col="cadetblue1")
  lines(AVGdotMOIZ[[i]],col="blue",lwd=4)
}


#*****************SHAFT DEFLECTION/SHAFT SPEED*********************
marker2a=V_CenterFace
marker2b=V_CenterFaceGR
Time2=Time


#function to filter out outliers
outliersFUN=function(list){
  quantiles=quantile(list,probs=c(.001,.99))
  list=subset(list,list>quantiles[1] & list<quantiles[2])
  return(list)
}

#find xyz coordinates between both points
deflection=mapply(function(marker1,marker2){marker2-marker1},marker1=marker2a,marker2=marker2b)
#calculate deflection in body coordinate system
deflection_body=rep(list(NA),length(deflection))
for (i in 1:length(deflection)){
  deflection_body[[i]]=matrix(NA,nrow=length(Rmat[[i]]),ncol=3)
  for (j in 1:length(Rmat[[i]])){
    x=data.matrix(deflection[[i]][j,])
  deflection_body[[i]][j,]=t(Rmat[[i]][[j]]%*%data.matrix(deflection[[i]][j,]))
}}
#input names of clubs
names(deflection_body)=paste(temp3)

deflectionX=lapply(deflection_body,function(x){x[,1]})
deflectionY=lapply(deflection_body,function(x){x[,2]})
deflectionZ=lapply(deflection_body,function(x){x[,3]})
#calculate magnitude of this distance
deflectionMAG=lapply(deflection,magFUN)

#filter out outliers on each swing
deflectionMAG=lapply(deflectionMAG,outliersFUN)
deflectionX=lapply(deflectionX,outliersFUN)
deflectionY=lapply(deflectionY,outliersFUN)
deflectionZ=lapply(deflectionZ,outliersFUN)


#create new time vector bc we changed the sizes due to outliers
Time2=mapply(function(time,list){time=head(time,-(length(time)-length(list)))},time=Time2,list=deflectionMAG)
ShaftSpeedMAG=mapply(FDAfunVelNEW,swing=deflection,Time=Time,nbasis=100,pos=4)
ShaftSpeedX=mapply(FDAfunVelNEW,swing=deflection,Time=Time,nbasis=100,pos=1)
ShaftSpeedY=mapply(FDAfunVelNEW,swing=deflection,Time=Time,nbasis=100,pos=2)
ShaftSpeedZ=mapply(FDAfunVelNEW,swing=deflection,Time=Time,nbasis=100,pos=3)


#Time lists/vectors prep for FDA --> DEFLECTION TIME VECTORS
TimeMS=lapply(Time2,function(x){x*1000})
Size=lapply(TimeMS,function(x){length(x)})

#Time lists/vectors prep for FDA --> SHAFT SPEED TIME VECTORS
TimeMSxs=lapply(Time,function(x){x*1000})
Sizexs=lapply(ShaftSpeedMAG,function(x){length(x)})


#----------FDA------------
DeflectionX_smoothed=FDAtransFUN(deflectionX,Size,players,pswings,1000,120)
DeflectionY_smoothed=FDAtransFUN(deflectionY,Size,players,pswings,1000,120)
DeflectionZ_smoothed=FDAtransFUN(deflectionZ,Size,players,pswings,1000,120)
Deflection_smoothed=FDAtransFUN(deflectionMAG,Size,players,pswings,1000,120)

ShaftSpeed_smoothed=FDAtransFUN(ShaftSpeedMAG,Sizexs,players,pswings,300,120)
ShaftSpeedX_smoothed=FDAtransFUN(ShaftSpeedX,Sizexs,players,pswings,300,120)
ShaftSpeedY_smoothed=FDAtransFUN(ShaftSpeedY,Sizexs,players,pswings,300,120)
ShaftSpeedZ_smoothed=FDAtransFUN(ShaftSpeedZ,Sizexs,players,pswings,300,120)


AVGDeflectionX=lapply(DeflectionX_smoothed,function(x){mean(x)})
AVGDeflectionY=lapply(DeflectionY_smoothed,function(x){mean(x)})
AVGDeflectionZ=lapply(DeflectionZ_smoothed,function(x){mean(x)})
AVGDeflection=lapply(Deflection_smoothed,function(x){mean(x)})
AVGShaftSpeed=lapply(ShaftSpeed_smoothed,function(x){mean(x)})
AVGShaftSpeedX=lapply(ShaftSpeedX_smoothed,function(x){mean(x)})
AVGShaftSpeedY=lapply(ShaftSpeedY_smoothed,function(x){mean(x)})
AVGShaftSpeedZ=lapply(ShaftSpeedZ_smoothed,function(x){mean(x)})


#plot all directions per player
plot.fd(Deflection_smoothed[[1]],col="coral2",ylim=c(-120,160))
lines(AVGDeflection[[1]],col="red",lwd=4)
lines(DeflectionX_smoothed[[1]],col="cadetblue")
lines(AVGDeflectionX[[1]],col="blue",lwd=4)
lines(DeflectionY_smoothed[[1]],col="darkorchid1")
lines(AVGDeflectionY[[1]],col="darkviolet",lwd=4)
lines(DeflectionZ_smoothed[[1]],col="darkolivegreen1")
lines(AVGDeflectionZ[[1]],col="darkgreen",lwd=4)
legend(1,120,legend=c("DeflectionMAG","Deflection X","Deflection Y","Deflection Z"),col=c("red","blue","darkviolet","darkgreen"),lwd=4,cex=1)

#plot all shaft speeds per player
plot.fd(ShaftSpeed_smoothed[[1]],col="coral2",ylim=c(-5,5))
lines(AVGShaftSpeed[[1]],col="red",lwd=4)
lines(ShaftSpeedX_smoothed[[1]],col="cadetblue")
lines(AVGShaftSpeedX[[1]],col="blue",lwd=4)
lines(ShaftSpeedY_smoothed[[1]],col="darkorchid1")
lines(AVGShaftSpeedY[[1]],col="darkviolet",lwd=4)
lines(ShaftSpeedZ_smoothed[[1]],col="darkolivegreen1")
lines(AVGShaftSpeedZ[[1]],col="darkgreen",lwd=4)
legend(1,5,legend=c("ShaftSpeedMAG","ShaftSpeed X","ShaftSpeed Y","ShaftSpeed Z"),col=c("red","blue","darkviolet","darkgreen"),lwd=4,cex=1)


plot.fd(ShaftSpeed_smoothed[[1]],col="cadetblue1")
lines(AVGShaftSpeed[[1]],col="blue",lwd=4)

#ALL PLAYERS
par(mfrow=c(4,4))
for (i in 1:length(players)){
  plot.fd(Deflection_smoothed[[i]],col="coral2",main=players[i],ylim=c(-100,200))
  lines(AVGDeflection[[i]],col="red",lwd=4)
  lines(DeflectionX_smoothed[[i]],col="cadetblue")
  lines(AVGDeflectionX[[i]],col="blue",lwd=4)
  lines(DeflectionY_smoothed[[i]],col="darkorchid1")
  lines(AVGDeflectionY[[i]],col="darkviolet",lwd=4)
  lines(DeflectionZ_smoothed[[i]],col="darkolivegreen1")
  lines(AVGDeflectionZ[[i]],col="darkgreen",lwd=4)
  #legend(1,120,legend=c("DeflectionMAG","Deflection X","Deflection Y","Deflection Z"),col=c("red","blue","darkviolet","darkgreen"),lwd=4,cex=.3)
}

par(mfrow=c(4,4))
for (i in 1:length(players)){
  plot.fd(ShaftSpeed_smoothed[[i]],col="cadetblue1",main=players[i],ylim=c(-4,6))
  lines(AVGShaftSpeed[[i]],col="blue",lwd=4)
  lines(ShaftSpeedX_smoothed[[i]],col="cadetblue")
  lines(AVGShaftSpeedX[[i]],col="blue",lwd=4)
  lines(ShaftSpeedY_smoothed[[i]],col="darkorchid1")
  lines(AVGShaftSpeedY[[i]],col="darkviolet",lwd=4)
  lines(ShaftSpeedZ_smoothed[[i]],col="darkolivegreen1")
  lines(AVGShaftSpeedZ[[i]],col="darkgreen",lwd=4)
}



#SORTED BY PLAYER AND CLUB
ANALYSISdeflection=clubanalysisFUN(deflectionMAG,Size,clubs,120,1,2)
ANALYSISshaftspeed=clubanalysisFUN(ShaftSpeedMAG,Sizexs,clubs,120,1,2)


#**********SPEED/ACCELERATION**********

#Markers to be used (error is high on grip markers for now)
Marker1=V_CenterFace
TimeN=lapply(Time,function(x){x*1000})
SizeS=lapply(TimeN,function(x){length(x)})

#2 methods to find impact speed->poly fit last few frames/fda fit raw derivative

#calculate velocity of specified marker for every swing
VTot=mapply(function(Marker1,Time){cbind(diff(Marker1[,1])/diff(Time),diff(Marker1[,2])/diff(Time),diff(Marker1[,3])/diff(Time))},
            Marker1=Marker1,Time=TimeN)
VTot=lapply(VTot,function(VTot){magFUN(VTot)})

IMPACT=mapply(time2impact,swing=MoCapALL,marker=Marker1,time=Time,frames=15,order=2,vel=1,SIMPLIFY = FALSE)

#Add new impact values to the velocity/time vectors
VTot=mapply(function(x,y){x=append(x,y[2],after=length(x))},x=VTot,y=IMPACT)
TimeN=mapply(function(x,y){x=append(x,y[1],after=(length(x)-1))},x=TimeN,y=IMPACT)
TimeN=lapply(TimeN,function(x){x=x[-length(x)]})


#FILTER/CLEAN THE DATA#################################################################
#filter out the worst outliers easy to pick out
TimeN=mapply(SimpleScrubFUN,V=VTot,TimeN=TimeN,num=2)
VTot=mapply(SimpleScrubFUN,V=VTot,TimeN=TimeN,num=1)

#acceleration filter method
ATot=mapply(function(VTot,Time){diff(VTot)/diff(Time)},VTot=VTot,Time=TimeN)
#TimeN2=lapply(TimeN,function(x){x=x[-1]})

TimeN=mapply(AccelScrubFUN,ATot=ATot,V=VTot,TimeN=TimeN,num=2)
VTot=mapply(AccelScrubFUN,ATot=ATot,V=VTot,TimeN=TimeN,num=1)

ATot=mapply(function(VTot,Time){diff(VTot)/diff(Time)},VTot=VTot,Time=TimeN)

#if rerunning code reinitialize ManualScrub
#****ManualScrub=0
#now run manual function

#*****ManualScrub=mapply(ManualScrubFUN,ATot=ATot,VTot=VTot)

#VTot=mapply(ManualReplaceFUN,List=VTot,ScrubValList=ManualScrub)
#TimeN=mapply(ManualReplaceFUN,List=TimeN,ScrubValList=ManualScrub)


#resize
SizeS=lapply(VTot,function(x){length(x)})

#Impact Speed
polyfitspeed=lapply(VTot,function(x){tail(x,n=1)*2.23694})
polyfitspeed=t(data.frame(polyfitspeed))


#SORTED BY PLAYER AND CLUB
#ANALYSISspeed=clubanalysisFUNavg(VTot,SizeS,clubs,2,1)




#--------FDA--------

#produces average for each player ALL SWINGS
Size=lapply(TimeN,function(x){length(x)})
#try to incorporate VTot instead of putting data through 2 smoothing functions
Velocity_smoothed=FDAtransFUN(VTot,Size,players,pswings,1000,101)
Acceleration_smoothed=FDAtransFUN(Accfda,Size,players,pswings,400,120)

#stdev
#sdavg=lapply(Velocity_smoothed,function(x){stdev.fd(x*2.23694)})
#plot(stdev.fd(Velocity_smoothed[[1]]*2.23694))


#AVG of all smoothed swings
AVGvel=lapply(Velocity_smoothed,function(x){mean(x)})
AVGaccel=lapply(Acceleration_smoothed,function(x){mean(x)})

#Velocity
plot(Velocity_smoothed[[1]])
lines(AVGvel[[1]],col="red",lwd=3)
AVGvel[[1]]$coefs[120]
#Acceleration
plot(Acceleration_smoothed[[1]])
lines(AVGaccel[[1]],col="red",lwd=3)

plot(Velocity_smoothed[[11]],xlim=c(980,1000),ylim=c(40,45))
lines(AVGvel[[11]],col="red",lwd=3)
AVGvel[[11]]$coefs[120]


for (i in 1:15){
plot(Velocity_smoothed[[i]],xlim=c(980,1000),main=players[[i]])
lines(AVGvel[[i]],col="red",lwd=3)
AVGvel[[i]]$coefs[101]
}

i=11
j=1
plot(Velocity_smoothed[[i]][j],xlim=c(980,1000))
lines(AVGvel[[i]][j],col="red",lwd=3)
AVGvel[[i]][j]$coefs[101]



#ERROR########### fda compared to our most accurate speed calculation (polyfitspeed)
ErrorVelocitySmoothed=rep(NA,length(polyfitspeed))
k=1
for (j in 1:length(Velocity_smoothed)){
  for (i in 1:ncol(Velocity_smoothed[[j]]$coefs)){
  ErrorVelocitySmoothed[k]=abs(Velocity_smoothed[[j]][i]$coefs[101]*2.23694-polyfitspeed[k])
k=k+1
  }
}


AVGErrorVelSmooth=mean(ErrorVelocitySmoothed)
maxerror=ErrorVelocitySmoothed[which.max(ErrorVelocitySmoothed)]




#TIME VARIABLES - TIME RATIOS OF THE SWING

#1: Ratio of Back Swing Time to Down Swing Time
AVGtest=lapply(AVGvel,function(x){return(x$coef)})

#identify 2 pts
#1st pt represents start of swing
#2nd pt represents start of downswing
bad2=lapply(AVGtest,
           function(x){
             plot(x)
             lines(x)
             return(identify(x))})

#ratio of back swing time to down
Swing_Time_Ratio = mapply(function(bad2,AVGtest){(bad2[2]-bad2[1])/(nrow(AVGtest)-bad2[2])},bad2=bad2,AVGtest=AVGtest)
hist(Swing_Time_Ratio)
pl=data.frame(Swing_Time_Ratio)

##############################################-----CLUSTER FUN-----###################################################


#import
allswings=dotMOIZ_smoothed
AVGswings=AVGdotMOIZ

allswings=lapply(allswings,function(x){return(x$coefs)})
AVGswings=lapply(AVGswings,function(x){return(x$coefs)})
names(allswings)=paste(players)
names(AVGswings)=paste(players)

#identify cut points
bad=lapply(AVGswings,
           function(x){
             plot(x)
             lines(x)
             return(identify(x))})

#ALL SWINGS
allswings=mapply(function(list,bad){list=list[bad:nrow(list),]},list=allswings,bad=bad)
#AVERAGE
AVGswings=mapply(function(list,bad){list=list[bad:nrow(list),]},list=AVGswings,bad=bad)

fdata=clusterFUN2(AVGswings, bad, 1)

# plot all the smoothed curves
plot.fd(fdata, col = "red")


# Cluster curves into 3 groups (trying to match male and female)
k = 2  #number of clusters
results = funclust(fdata, K = k) #function that performs the actual clustering
pred_cluster = results$cls #predicted group from clustering
pred_cluster
pl=data.frame(pred_cluster)

plot.fd(fdata,col=pred_cluster)










clusterFUN2 = function(allswings, bad, val){
  
  if (val==1){
    list_lengths=lapply(allswings,length)
  }else{
    list_lengths=lapply(allswings,nrow)}
  
  thresh=unlist(list_lengths[which.max(list_lengths)])
  thresh=thresh+1
  
  #find the position in the array/list with the largest number of points
  sizebad=lapply(bad,function(bad){return(thresh-bad)})
  maxval=unlist(sizebad[which.max(unlist(sizebad))])
  
  tvec2=mapply(function(list,maxval){
    if (val==0){
      test=list[,1]
    }else{
      test=list}
    tstep=(maxval-length(test))/length(test)+1
    tvec=seq(from=tstep,to=maxval,by=tstep)
    return(tvec)
  },list=allswings,maxval=maxval)
  
  
  #this function takes in a feature/time with the same size and outputs the new values with a new index from 1 to specified size
  interp=rep(list(NA),length(allswings))
  
  if (val==0){
    for (j in 1:length(allswings)){ 
      interp[[j]]=matrix(NA,nrow=maxval-1,ncol=ncol(allswings[[j]]))
      for (i in 1:ncol(allswings[[j]])){
        swing=allswings[[j]][,i]
        oldtime=tvec2[[j]]
        newtime=seq(from=1,to=maxval,by=1)
        
        interp[[j]][,i]=approx(oldtime,swing,xout=newtime,method="linear")$y[2:maxval]
        interp[[j]][maxval-1,i]=swing[length(swing)]
      }}
  }else{
    for (j in 1:length(allswings)){ 
      swing=allswings[[j]]
      oldtime=tvec2[[j]]
      newtime=seq(from=1,to=maxval,by=1)
      
      interp[[j]]=approx(oldtime,swing,xout=newtime,method="linear")$y[2:maxval]
      interp[[j]][maxval-1]=swing[length(swing)]
    }
  }
  
  allswings2=list.cbind(interp)
  tvec3=seq(from=1,to=maxval-1,by=1)
  
  # create spline smoothed representations of the curves
  splines = create.bspline.basis(rangeval = c(tvec3[1], max(tvec3)) , nbasis = 25, norder = 6)
  fdata = Data2fd(y = allswings2, argvals = tvec3, basisobj = splines)
  return(fdata)
}





# #find the longest swing in the list of swings
# list_lengths=lapply(allswings,nrow)
# thresh=unlist(list_lengths[which.max(list_lengths)])
# thresh=thresh+1
# 
# #find the position in the array/list with the largest number of points
# sizebad=lapply(bad,function(bad){return(thresh-bad)})
# maxval=unlist(sizebad[which.max(unlist(sizebad))])
# 
# tvec2=mapply(function(list,maxval){
#   test=list[,1]
#   tstep=(maxval-length(test))/length(test)+1
#   tvec=seq(from=tstep,to=maxval,by=tstep)
#   return(tvec)
# },list=allswings,maxval=maxval)
# 
# 
# 
# #allswings2=mapply(interpolationlistFUN(allswings,tvec2,size),allswings=allswings,tvec2=tvec2,size=120)
# 
# 
# #this function takes in a feature/time with the same size and outputs the new values with a new index from 1 to specified size
# interp=rep(list(NA),length(allswings))
# 
# for (j in 1:length(allswings)){ 
#   interp[[j]]=matrix(NA,nrow=maxval-1,ncol=ncol(allswings[[j]]))
#   for (i in 1:ncol(allswings[[j]])){
#     swing=allswings[[j]][,i]
#     oldtime=tvec2[[j]]
#     newtime=seq(from=1,to=maxval,by=1)
#     
#     interp[[j]][,i]=approx(oldtime,swing,xout=newtime,method="linear")$y[2:maxval]
#     interp[[j]][maxval-1,i]=swing[length(swing)]
#   }
# }
# 
# allswings2=list.cbind(interp)
# tvec3=seq(from=1,to=maxval-1,by=1)
# 
# # create spline smoothed representations of the curves
# splines = create.bspline.basis(rangeval = c(tvec3[1], max(tvec3)) , nbasis = 25, norder = 6)
# fdata = Data2fd(y = allswings2, argvals = tvec3, basisobj = splines)
# 
# # plot all the smoothed curves
# plot.fd(fdata, col = "red")
# 
# 
# 
# # Cluster curves into 3 groups (trying to match male and female)
# k = 2   #number of clusters
# results = funclust(fdata, K = k) #function that performs the actual clustering
# pred_cluster = results$cls #predicted group from clustering
# pred_cluster
# 
# plot.fd(fdata,col=pred_cluster)
























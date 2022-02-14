#library("ggplot2")
library("getopt")
library("data.table")
library(Cairo)
spec <- matrix(
  c("jiange","k",2,"numeric", "The numer of bins!",
    "len","l",2,"numeric", "Use len-bin!",
    "name","n",1,"character", "genome",
    "out","o",1,"character" ,"out_dir",
    "images_out","g",1,"character" ,"images_out_dir"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
len=opt$len
jiange=opt$jiange
name=opt$name
out=opt$out
images_out=opt$images_out
path=paste(out,name,sep="/")

images_path=paste(images_out,name,sep="/")
data=as.matrix(read.table(paste(path,".site.txt",sep=""),header = F,sep="\t",stringsAsFactors = F))
info=read.table(paste(path,".genome_info.txt",sep=""),header=F,sep = "\t",stringsAsFactors = F)
iden=as.matrix(fread(paste(path,".inentity.txt",sep=""),header = F,sep="\t"))
if(is.null(jiange)){
  if(as.matrix(info[2])>len){
    qujian=seq(0,info[2]$V2,by=len) 
  }
}else{
  d=info[2]$V2/jiange
  qujian=seq(0,info[2]$V2,by=d)
}


if(qujian[length(qujian)]<as.matrix(info[2])){
  qujian[length(qujian)+1]=qujian[length(qujian)]+qujian[2]
}
uninum=hist(data[data[,2]==1,1],breaks = qujian)
mutinum=hist(data[data[,2]==2,1],breaks = qujian)


tongji=c()
for(i in 1:(length(qujian)-1)){
  big=iden[which(iden[,1]>qujian[i]),]
  small=big[which(big[,1]<qujian[i+1]),2]
  temp=c((qujian[i+1]+qujian[i])/2,mean(small))
  if(length(small)==0){
    temp=c((qujian[i+1]+qujian[i])/2,0)
  }
  tongji=rbind(tongji,temp)
}
covearge=round(info[4],5)*100
m_title=paste(gsub("_"," ",info[5]),"\nCoverage: ",covearge,"%",paste="")
tongji[,2]=tongji[,2]*100

getmax=function(x){
  i=0
  while (x>10) {
    x=x/10
    i=i+1
  }
  if(round(x)<x){
    x=round(x)+1
  }else{
    x=round(x)
  }
  x=x*10^i
  return(c(x,i))
}

maxy=getmax(max(max(uninum$counts),max(mutinum$counts)))[1]
i=getmax(max(max(uninum$counts),max(mutinum$counts)))[2]

if(info[2]>1000000){
  ticknum=50
  attick=seq(0,qujian[length(qujian)],by=qujian[length(qujian)]/ticknum)
  ticklabel=0
  newat=0
  for (i in 1:10) {
    ticklabel=c(ticklabel,paste(round(attick[i*5+1]/1000000,1),"M",sep=""))
    newat=c(newat,attick[i*5+1])
  }
}else{
  ticknum=50
  attick=seq(0,qujian[length(qujian)],by=qujian[length(qujian)]/ticknum)
  ticklabel=0
  newat=0
  for (i in 1:10) {
    ticklabel=c(ticklabel,paste(round(attick[i*5+1]/1000,1),"K",sep=""))
    newat=c(newat,attick[i*5+1])
  }
}

CairoPDF(file=paste(images_path,"_identity.pdf",sep=""),width=12,height=6)
#pdf(file=paste(images_path,"_identity.pdf",sep=""),width = 12)
#map<-dev.cur()
#png(file=paste(images_path,"_identity.png",sep=""),type = "cairo",width=750*3,height=3*500,res=72*3)
#dev.control("enable")

par(mai=c(0.5,0.5,0.8,0.5),xpd=F)
par(xaxs = "i", yaxs = "i",oma=c(4,4,4,4))
#Identity Y.axis

plot(tongji,type="l",ylim=c(0,100),col="#B34BF3",yaxt="n",xaxt="n",xlim=c(0,max(qujian)),xlab="",ylab = "",lwd=2.5,bty="n")
#axis(side=4,las=1,cex.lab=1.5,cex.axis=1.5,col="#B34BF3",lwd=2,col.ticks = "#B34BF3",col.axis="#B34BF3",col.lab="#B34BF3")
axis(side=4,las=1,cex.lab=1.5,cex.axis=1.5,col="#B34BF3",col.axis="#B34BF3",lwd=2)
for (i in seq(0,100,by=20)) {
  abline(h=i,lty=2,xpd=F,lwd=1.2,col="orange")
}
mtext(4,text="Identity",cex = 1.5,col="#B34BF3",outer = T,line=1.5)

par(new=T,xaxs = "i", yaxs = "i",xpd=F)
#Mapped reads  Y.axis
plot(c(tongji[1,1],tongji[1:length(uninum$counts),1],rev(tongji[1:length(uninum$counts),1]),tongji[1,1]),c(0,uninum$counts,rep(0,length(uninum$counts)),0),type="l",ylim=c(0,maxy),xlim=c(0,max(qujian)),col=rgb(0,138,0,alpha=150,max=255),yaxt="n",xaxt="n",cex.lab=1.5,cex.axis=1.5,bty="n",lwd=2.5,ylab="",xlab="Genomic position")
axis(2,col="#008A00",at=seq(0,maxy,by=10^getmax(max(max(uninum$counts),max(mutinum$counts)))[2]),las=1,col.ticks = "#008A00",col.axis="#008A00",cex.lab=1.2,cex.axis=1.2,lwd=2,col.lab="#008A00")
mtext(2,text="Mapped reads",cex = 1.5,col="#008A00",outer = T,line=1.4)
polygon(c(tongji[1:length(uninum$counts),1],rev(tongji[1:length(uninum$counts),1])),c(uninum$counts,rep(0,length(uninum$counts))),col=rgb(0,138,0,alpha=50,max=255),border = F,bty="u")
#grid(col="orange",lwd=1.5,lty=2)
#X.axis
par(new=T,xaxs = "i", yaxs = "i",xpd=F)
plot(c(tongji[1,1],tongji[1:length(mutinum$counts),1],rev(tongji[1:length(mutinum$counts),1]),tongji[1,1]),c(0,mutinum$counts,rep(0,length(mutinum$counts)),0),type="l",ylim=c(0,maxy),xaxt="n",xlim=c(0,max(qujian)),col=rgb(52,152,219,alpha=150,max=255),cex.lab=1.5,yaxt="n",cex.axis=1.5,bty="n",lwd=2.5,ylab="",xlab="Genomic position")

#axis(1,at=attick,labels = F,cex.lab=1.2,cex.axis=1.2,lwd=2)
axis(1,at=newat,labels = ticklabel,cex.lab=1.2,cex.axis=1.2,lwd=2)
for (i in 1:10) {
  abline(v=newat[i+1],lty=2,xpd=F,lwd=1.2,col="orange")
}
mtext(1,text="Genomic position",cex = 1.5,line=1,outer = T)
polygon(c(tongji[1:length(mutinum$counts),1],rev(tongji[1:length(mutinum$counts),1])),c(mutinum$counts,rep(0,length(mutinum$counts))),col=rgb(52,152,219,alpha=50,max=255),border = F,bty="u")
lines(x=c(0,info[2]$V2),y=c(0,0),lwd=2)
#lines(x=c(0,max(qujian)),y=c(0,0),lwd=2)
legend(x=max(qujian)/7,y=(maxy*7/6),xpd = T,bty="n",horiz=T,legend =c("Multiple Mapping","Unique Mapping","Identity"),cex=1.2,lwd=c(2.6,2.6,2.6),seg.len=c(2,2,2),lty=c(1,1,1),col=c(rgb(52,152,219,max=255),"#008A00","#B34BF3"),text.col =c(rgb(52,152,219,max=255),"#008A00","#B34BF3") )
mtext(m_title,side=3,outer = T,cex=1.5)
#dev.copy(which=map)
#dev.off()
dev.off()
result_data=cbind(tongji[1:length(uninum$counts),1],uninum$counts,mutinum$counts)
write.table(result_data,paste(path,"_counts.txt",sep=""),quote = F,sep="\t",col.names = c("breakpoints","unique_map_num","mutiple_map_num"),row.names = F)

CairoPNG(file=paste(images_path,"_identity.png",sep=""),width=960,height=480,pointsize=13)
#png(file=paste(images_path,"_identity.png",sep=""),type = "cairo",width=700*3,height=3*500,res=72*3)

par(mai=c(0.5,0.5,0.8,0.5),xpd=F)
par(xaxs = "i", yaxs = "i",oma=c(4,4,4,4))
#Identity Y.axis
#plot(tongji,type="l",ylim=c(0,100),col="#B34BF3",yaxt="n",xaxt="n",xlim=c(1,info[2]$V2)),xlab="",ylab = "",lwd=2.5,bty="n")

plot(tongji,type="l",ylim=c(0,100),col="#B34BF3",yaxt="n",xaxt="n",xlim=c(0,max(qujian)),xlab="",ylab = "",lwd=2.5,bty="n")
#axis(side=4,las=1,cex.lab=1.5,cex.axis=1.5,col="#B34BF3",lwd=2,col.ticks = "#B34BF3",col.axis="#B34BF3",col.lab="#B34BF3")
axis(side=4,las=1,cex.lab=1.5,cex.axis=1.5,col="#B34BF3",col.axis="#B34BF3",lwd=2)
for (i in seq(0,100,by=20)) {
  abline(h=i,lty=2,xpd=F,lwd=1.2,col="orange")
}
mtext(4,text="Identity",cex = 1.5,col="#B34BF3",outer = T,line=1.5)

par(new=T,xaxs = "i", yaxs = "i",xpd=F)
#Mapped reads  Y.axis
plot(c(tongji[1,1],tongji[1:length(uninum$counts),1],rev(tongji[1:length(uninum$counts),1]),tongji[1,1]),c(0,uninum$counts,rep(0,length(uninum$counts)),0),type="l",ylim=c(0,maxy),xlim=c(0,max(qujian)),col=rgb(0,138,0,alpha=150,max=255),yaxt="n",xaxt="n",cex.lab=1.5,cex.axis=1.5,bty="n",lwd=2.5,ylab="",xlab="Genomic position")
axis(2,col="#008A00",at=seq(0,maxy,by=10^getmax(max(max(uninum$counts),max(mutinum$counts)))[2]),las=1,col.ticks = "#008A00",col.axis="#008A00",cex.lab=1.2,cex.axis=1.2,lwd=2,col.lab="#008A00")
mtext(2,text="Mapped reads",cex = 1.5,col="#008A00",outer = T,line=1.4)
polygon(c(tongji[1:length(uninum$counts),1],rev(tongji[1:length(uninum$counts),1])),c(uninum$counts,rep(0,length(uninum$counts))),col=rgb(0,138,0,alpha=50,max=255),border = F,bty="u")
#grid(col="orange",lwd=1.5,lty=2)
#X.axis
par(new=T,xaxs = "i", yaxs = "i",xpd=F)
plot(c(tongji[1,1],tongji[1:length(mutinum$counts),1],rev(tongji[1:length(mutinum$counts),1]),tongji[1,1]),c(0,mutinum$counts,rep(0,length(mutinum$counts)),0),type="l",ylim=c(0,maxy),xaxt="n",xlim=c(0,max(qujian)),col=rgb(52,152,219,alpha=150,max=255),cex.lab=1.5,yaxt="n",cex.axis=1.5,bty="n",lwd=2.5,ylab="",xlab="Genomic position")

#axis(1,at=attick,labels = F,cex.lab=1.2,cex.axis=1.2,lwd=2)
axis(1,at=newat,labels = ticklabel,cex.lab=1.2,cex.axis=1.2,lwd=2)
for (i in 1:10) {
  abline(v=newat[i+1],lty=2,xpd=F,lwd=1.2,col="orange")
}
mtext(1,text="Genomic position",cex = 1.5,line=1,outer = T)
polygon(c(tongji[1:length(mutinum$counts),1],rev(tongji[1:length(mutinum$counts),1])),c(mutinum$counts,rep(0,length(mutinum$counts))),col=rgb(52,152,219,alpha=50,max=255),border = F,bty="u")
lines(x=c(0,info[2]$V2),y=c(0,0),lwd=2)
#lines(x=c(0,max(qujian)),y=c(0,0),lwd=2)
legend(x=max(qujian)/7,y=(maxy*7/6),xpd = T,bty="n",horiz=T,legend =c("Multiple Mapping","Unique Mapping","Identity"),cex=1.2,lwd=c(2.6,2.6,2.6),seg.len=c(2,2,2),lty=c(1,1,1),col=c(rgb(52,152,219,max=255),"#008A00","#B34BF3"),text.col =c(rgb(52,152,219,max=255),"#008A00","#B34BF3") )
mtext(m_title,side=3,outer = T,cex=1.5)
dev.off()

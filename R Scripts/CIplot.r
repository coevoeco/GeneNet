CIplot<-function(data,level,out=F,outname="",xlabel="x",ylabel="y")
{
	colnames(data)=c('x','y')
	mod=lm(y~x,data=data)
	newx=seq(min(data$x)-10,max(data$x)+10,length.out=1000)
	pred=predict(mod,newdata=data.frame(x=newx),interval='confidence',level=level)
	realpred=predict(mod,newdata=data.frame(x=data$x),interval='confidence',level=level)
	outliervec=rep(0,length(data$x))
	for(i in 1:length(data$x)){if(data$y[i]>realpred[i,3]|data$y[i]<realpred[i,2]){outliervec[i]=1}}
	outliers=which(outliervec==1)
	plot(y~x,data=data,pch=16,cex=.65,xlab=xlabel,ylab=ylabel)
	polygon(c(rev(newx),newx),c(rev(pred[,3]),pred[,2]),col=rgb(.1,.1,.1,.3),border=NA)
	lines(newx,pred[,2],lty=2,lwd=.5,col=2)
	lines(newx,pred[,3],lty=2,lwd=.5,col=2)
	lines(newx,pred[,1],lty=2,lwd=.5)
	points(data$x[outliers],data$y[outliers],pch=16,col=2,cex=.65)

	if(out==T)
	{
		pdf(outname,height=6,width=6)
		plot(y~x,data=data,pch=16,cex=.65,xlab=xlabel,ylab=ylabel)
		polygon(c(rev(newx),newx),c(rev(pred[,3]),pred[,2]),col=rgb(.1,.1,.1,.3),border=NA)
		lines(newx,pred[,2],lty=2,lwd=.5,col=2)
		lines(newx,pred[,3],lty=2,lwd=.5,col=2)
		lines(newx,pred[,1],lty=2,lwd=.5)
		points(data$x[outliers],data$y[outliers],pch=16,col=2,cex=.65)
		dev.off()
	}

	results=NULL
	results$pred=pred
	results$realpred=realpred
	results$outliers=outliers
	results$mod=mod

	return(results)
}

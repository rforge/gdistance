WeitzmanCalc <- function(binarySeq,distObject) 
{
	n <- attr(distObject, "Size")
	distValue <- 0
	for(i in 1:(length(binarySeq)))
	{
		distValue <- distValue + min(distObject)
		if(i<(length(binarySeq)))
		{
			Pair <- matrIndex(which(distObject==min(distObject)),n+1-i)
			Link <- Pair[binarySeq[i] +1]
			if(Link>1){index1 <- distIndex(Link,1:(Link-1),n+1-i)}else(index1<-NA)
			if(Link<(n+1-i)){index2 <- distIndex((Link+1):(n+1-i),Link,n+1-i)}else(index2<-NA)
			index <- na.omit(c(index1,index2))
			distObject <- distObject[-index]
		}
	}
	return(distValue)
}

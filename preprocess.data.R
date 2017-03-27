preprocess.data <-function(Xp,X,Preprocess.net){
	if(class(Preprocess.net)!='rarb.preprocess.net') stop('class(Preprocess.net) is wrong, use function preprocess.net to generate it.')
	if(class(X)!='matrix'|class(Xp)!='matrix') stop('class(X) and class(Xp) should be \'matrix\'.')

	Edge.tree=Preprocess.net$Edge.tree
	net.vertex=names(Edge.tree)
	Xp=Xp[,net.vertex]
	X=X[,net.vertex]

	Sigma=lapply(Edge.tree,function(x)cov(Xp[,x,drop=FALSE],use="pairwise.complete.obs")+cov(X[,x,drop=FALSE],use="pairwise.complete.obs"))
	mu=colMeans(Xp,na.rm=TRUE)-colMeans(X,na.rm=TRUE)
	MU=lapply(Edge.tree,function(x)mu[unlist(x)])
	print('finish computing Sigma and MU.',quote=FALSE)

	robj=list(Sigma=Sigma,MU=MU,Xp=Xp,X=X)
	class(robj)='rarb.preprocess.data'
	return(robj)
}

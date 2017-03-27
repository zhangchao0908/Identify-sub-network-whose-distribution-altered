preprocess.net <-function(edgem,net.vertex){
	if(class(edgem)!='matrix') stop('class(edgem) should be \'matrix\'.')
	if(mode(edgem)!='character') stop('mode(edgem) should be \'character\'.')
	if(ncol(edgem)!=2) stop('ncol(edgem) should be 2.')

	Edge.matrix=rbind(edgem, edgem[,2:1], cbind(net.vertex,net.vertex))
	dimnames(Edge.matrix)=NULL
	Edge.tree=lapply(net.vertex,function(v)Edge.matrix[Edge.matrix[,1]==v,2])
	names(Edge.tree)=net.vertex
	
	robj=list(edgem=edgem,Edge.tree=Edge.tree)
	class(robj)='rarb.preprocess.net'
	return(robj)
}

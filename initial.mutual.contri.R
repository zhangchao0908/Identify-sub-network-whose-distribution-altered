initial.mutual.contri <-function(Preprocess.net,Preprocess.data,shrink){
  if(class(Preprocess.net)!='rarb.preprocess.net') stop('class(Preprocess.net) is wrong, use function preprocess.net to generate it.')
  if(class(Preprocess.data)!='rarb.preprocess.data') stop('class(Preprocess.data) is wrong, use function preprocess.data to generate it.')
  
  ####################################################################
  ## assistant function: compute signal to noise ratio at vertex v. ##
  ####################################################################
  s2n <-function(v){
    S=Sigma[[v]]
    M=MU[[v]]
    if(ncol(S)==1) return(as.vector(M^2/S))
    else{
      robj=shrink*S+(1-shrink)*diag(diag(S))
      if(any(eigen(robj)$values<0))stop('\'shrink\' should be set smaller.')
      robj=(M%*%solve(robj)%*%M)[1]
      return(robj)
    }
  }
  
  #################################################################################
  ## assistant function: compute the contributions of the neighbors of vertex v. ##
  #################################################################################
  contri <-function(v){
    S=Sigma[[v]]
    M=MU[[v]]
    D=diag(S)
    n=length(M)
    
    if (n>2){
      robj=sapply(1:(n-1),function(j)S2N[v]-M[-j]%*%solve(shrink*S[-j,-j]+(1-shrink)*diag(D[-j]))%*%M[-j])
      names(robj)=names(M)[-n]
      return(robj)
    }
    if (n==2){
      robj=S2N[v]-M[2]^2/S[2,2]
      names(robj)=names(M)[1]
      return(robj)
    }
    if (n==1) return(S2N[v])
  }
  #################################################################################
  Sigma=Preprocess.data$Sigma
  MU=Preprocess.data$MU
  net.vertex=names(Sigma)
  
  S2N=sapply(net.vertex,s2n);
  print('finish computing S2N...',quote=FALSE)
  
  Contri=lapply(net.vertex,function(v)contri(v))
  names(Contri)=net.vertex
  print('finish computing Contri...',quote=FALSE)
  
  edgem=Preprocess.net$edgem
  Degree=degree(graph.data.frame(edgem,directed=FALSE))
  names(Degree)=net.vertex
  
  mC=apply(edgem,1,function(x){
    d1=Degree[x[1]]
    d2=Degree[x[2]]
    if(d1>1&d2>1)
      return(mean(Contri[[x[2]]][x[1]],Contri[[x[1]]][x[2]]))
    else
      return(min(Contri[[x[1]]][x[2]],Contri[[x[2]]][x[1]]))
  })
  print('finish computing mC.',quote=FALSE)
  
  robj=list(S2N=S2N,Contri=Contri,mC=mC,shrink=shrink)
  class(robj)='rarb.initial.mutual.contri'
  
  return(robj)
}

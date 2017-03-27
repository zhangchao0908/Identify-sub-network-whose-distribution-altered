rank.edge <-function(file,Preprocess.net,Preprocess.data,Initial.mutual.contri){
  if(class(Preprocess.net)!='rarb.preprocess.net') stop('class(Preprocess.net) is wrong, use function preprocess.net to generate it.')
  if(class(Preprocess.data)!='rarb.preprocess.data') stop('class(Preprocess.data) is wrong, use function preprocess.data to generate it.')
  if(class(Initial.mutual.contri)!='rarb.initial.mutual.contri') stop('class(Initial.mutual.contri) is wrong, use function initial.mutual.contri to generate it.')
  
  
  Edge.tree=Preprocess.net$Edge.tree
  edgem=Preprocess.net$edgem
  Sigma=Preprocess.data$Sigma
  MU=Preprocess.data$MU
  Xp=Preprocess.data$Xp
  X=Preprocess.data$X
  S2N=Initial.mutual.contri$S2N
  Contri=Initial.mutual.contri$Contri
  mC=Initial.mutual.contri$mC
  shrink=Initial.mutual.contri$shrink
  
  ####################################################################
  ## assistant function: compute signal to noise ratio at vertex v. ##
  ####################################################################
  s2n <-function(v){
    S=Sigma[[v]]
    M=MU[[v]]
    if(ncol(S)==1) return(as.vector(M^2/S))
    else return((M%*%solve(shrink*S+(1-shrink)*diag(diag(S)))%*%M)[1])
  }
  ####################################################################
  
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
  
  ###############################################################
  ## assistant function: find the minimum mutual contribution. ##
  ###############################################################
  find.min <-function(){
    if(any(is.na(mC))) stop('NA 2')
    w=which.min(mC)
    V=edgem[w,1]
    v=edgem[w,2]
    print(paste(length(mC),'::','minimum contribution is',mC[w]))
    cat(V,v,mC[w],sum(S2N),file=file,append=TRUE,sep='\t',fill=TRUE)
    edgem<<-edgem[-w,,drop=FALSE]
    mC<<-mC[-w]
    return(c(V,v))
  }
  ###############################################################
  
  #############################################################################
  ## assistant function: updata the structure of network and its attributes. ##
  #############################################################################
  renew.Sigma.Contri <-function(V,v){
    EV=Edge.tree[[V]]
    Ev=Edge.tree[[v]]
    not.v=v!=EV
    not.V=V!=Ev
    
    Edge.tree[[v]]<<-Ev[not.V]
    Edge.tree[[V]]<<-EV[not.v]
    Sigma[[v]]<<-Sigma[[v]][not.V,not.V,drop=FALSE] 
    Sigma[[V]]<<-Sigma[[V]][not.v,not.v,drop=FALSE]
    MU[[v]]<<-MU[[v]][not.V]
    MU[[V]]<<-MU[[V]][not.v]
    S2N[[v]]<<-s2n(v)
    S2N[[V]]<<-s2n(V)
    Contri[[v]]<<-contri(v)
    Contri[[V]]<<-contri(V)
    
    if(length(EV)>2|length(Ev)>2){
      if((length(EV)-1)!=length(EV[not.v]))stop('hehe')
      cham=rbind(cbind(rep(V,length(EV)-1),EV[not.v]),cbind(rep(v,length(Ev)-1),Ev[not.V]))
      cham=rbind(cham,cbind(cham[,2],cham[,1]))
      cham=unique(cham[cham[,1]>cham[,2],,drop=FALSE])
      
      apply(cham,1,function(x){
        d1=length(Edge.tree[[x[1]]])-1
        d2=length(Edge.tree[[x[2]]])-1
        if(d1>1&d2>1)
          mC[edgem[,1]==x[1]&edgem[,2]==x[2]]<<-mean(Contri[[x[2]]][x[1]],Contri[[x[1]]][x[2]])
        else
          mC[edgem[,1]==x[1]&edgem[,2]==x[2]]<<-min(Contri[[x[2]]][x[1]],Contri[[x[1]]][x[2]])
        
      })
    }
    return(NULL)
  }
  #############################################################################
  
  lapply(1:nrow(edgem),function(i){
    Vv=find.min()
    renew.Sigma.Contri(Vv[1],Vv[2])
  })
  return(NULL)
}


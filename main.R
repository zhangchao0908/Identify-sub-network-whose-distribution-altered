#### helper functions ####
within_boundary = function(vertex){
  vertex = as.numeric(unlist(strsplit(vertex, '_')))
  x = vertex[1]
  y = vertex[2]
  return(-2 <= x & x <= 2 & -2 <= y & y <= 2)
}
up = function(vertex){
  vertex = as.numeric(unlist(strsplit(vertex, '_')))
  x = vertex[1]
  y = vertex[2] + 1
  return(paste(x, y, sep = '_'))
}
down = function(vertex){
  vertex = as.numeric(unlist(strsplit(vertex, '_')))
  x = vertex[1]
  y = vertex[2] - 1
  return(paste(x, y, sep = '_'))
}
left = function(vertex){
  vertex = as.numeric(unlist(strsplit(vertex, '_')))
  x = vertex[1] - 1
  y = vertex[2]
  return(paste(x, y, sep = '_'))
}
right = function(vertex){
  vertex = as.numeric(unlist(strsplit(vertex, '_')))
  x = vertex[1] + 1
  y = vertex[2]
  return(paste(x, y, sep = '_'))
}
plot.edge = function(edges, ...){
  for(i in 1:nrow(edges)){
    v1 = as.numeric(unlist(strsplit(edges[i, 1], split = '_')))
    v2 = as.numeric(unlist(strsplit(edges[i, 2], split = '_')))
    segments(v1[1], v1[2], v2[1], v2[2], ...)
  }
}
plot.vertex = function(vertices, ...){
  for(vertex in vertices){
    vertex = as.numeric(unlist(strsplit(vertex, split = '_')))
    points(vertex[1], vertex[2], ...)
  }
}
source('preprocess.net.R')
source('preprocess.data.R')
source('initial.mutual.contri.R')
source('rank.edge.R')

#### build network ####
net = outer((-5):5,5:(-5),paste,sep='_')
sub_net = net[4:8,4:8]
vertex_with_changed_mean = outer((-2):2,(-2):2,paste,sep='_')[seq(1, 25, 2)]
vertex_with_changed_cor = outer((-2):2,(-2):2,paste,sep='_')[seq(2, 25, 2)]
edgem = lapply(1:10, function(i)rbind(net[,i:(i+1)], t(net[i:(i+1),])))
edgem = do.call(rbind, edgem)
sub_edgem = lapply(1:4, function(i)rbind(sub_net[,i:(i+1)], t(sub_net[i:(i+1),])))
sub_edgem = do.call(rbind, sub_edgem)
net.vertex = as.vector(net)

#### plot network ####
par(bg = 'gray', bty = 'n')
plot(c(-5, 5), c(-5, 5), type='n', axes=FALSE, xlab='', ylab='')
axis(side=1, at=(-5):5, labels=(-5):5, tick = FALSE)
axis(side=2, at=(-5):5, labels=(-5):5, tick = FALSE, las = 2)
plot.edge(edgem)
plot.edge(sub_edgem, lwd = 2)
plot.vertex(net.vertex, pch=19)
plot.vertex(vertex_with_changed_cor, col='yellow', cex=1.2, pch=19)
plot.vertex(vertex_with_changed_mean, col='red', cex=1.2, pch=19)

#### generate data ####
n_subject = 20
set.seed(0)
X = matrix(rnorm(11*11*n_subject), nrow = n_subject)
Xp = matrix(rnorm(11*11*n_subject), nrow = n_subject)
colnames(X) = colnames(Xp) = net.vertex

for(vertex in vertex_with_changed_mean){
  if(within_boundary(up(vertex))){
    Xp[, vertex] = 0.4*Xp[, vertex] + 0.9*Xp[, up(vertex)]
  }
  if(within_boundary(down(vertex))){
    Xp[, vertex] = 0.9*Xp[, vertex] + 0.4*Xp[, down(vertex)]
  }
  if(within_boundary(left(vertex))){
    Xp[, vertex] = 0.4*Xp[, vertex] + 0.9*Xp[, left(vertex)]
  }
  if(within_boundary(right(vertex))){
    Xp[, vertex] = 0.9*Xp[, vertex] + 0.4*Xp[, right(vertex)]
  }
  Xp[, vertex] = Xp[, vertex] + 1
}

#### rank edges ####
Preprocess.net <- preprocess.net(edgem,net.vertex)
Preprocess.data <- preprocess.data(Xp,X,Preprocess.net)
Initial.mutual.contri <- initial.mutual.contri(Preprocess.net,Preprocess.data,shrink=0.8)
Rank.edge = rank.edge('test.txt',Preprocess.net,Preprocess.data,Initial.mutual.contri)
result = read.table('test.txt',sep='\t', stringsAsFactors = FALSE)

vertices=c()
i = nrow(result)
while(length(vertices)<25){
  vertices = unique(c(vertices, unlist(result[i,1:2])))
  i = i-1
}
plot.vertex(vertices, cex = 1.4, col='green', pch = 19)

#### compare with T-test ####
# Tresult = sapply(net.vertex, function(vertex)t.test(X[,vertex],Xp[,vertex])$p.value)
# Tresult = names(sort(Tresult))[1:25]
# 
# plot.vertex(Tresult, cex = 1.4, col='blue', pch = 19)




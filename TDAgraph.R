library(igraph)
library(TDA)
library(TDAvec)
  
graphToPD <- function(g,maxSimplxDim=2,maxHomDim=2,nodeActFunction='degree',nodeAct=NULL){
  # g: an igraph graph 
  # maxSimplxDim: maximum dimension for simplices; choose from 1,2 or 3
  # maxHomDim: maximum homological dimension; set maxHomDim = maxSimpDim or maxSimpDim-1
  # nodeActFunction: type of node activation function; choose from "degree", "closeness", "betweenness" or "custom"
  # if nodeActFunction="custom", additionally pass a vector of quantitative node attributes to nodeAct
  
  # construct simplices of dimensions 0,...,maxHomDim
  if (maxSimplxDim==1){
    zeroSimplx <- as.list(V(g))
    oneSimplx<-data.frame(t(as_edgelist(g)))
    cmplx<-c(zeroSimplx,oneSimplx)
  } else 
    if (maxSimplxDim==2){
      zeroSimplx <- as.list(V(g))
      oneSimplx <- data.frame(t(as_edgelist(g)))
      twoSimplx <- data.frame(matrix(triangles(g),nrow = 3))
      cmplx<-c(zeroSimplx,oneSimplx,twoSimplx)
    } else 
      if (maxSimplxDim==3){ # use 3 for small/moderate sized graphs
        cmplx <- cliques(g,min=1,max=maxSimplxDim+1) # could be slow to run
      } else cat("The maximum dimension of simplices must not exceed 3")
      
  # compute node activation values
  if (nodeActFunction!="custom") nodeAct <- get(nodeActFunction)(g) 
  nodeAct <- nodeAct/max(nodeAct) # normalize to [0,1]
  
  # construct sublevel set filtration based on normalized node activation values
  fltn <- funFiltration(FUNvalues = nodeAct,cmplx = cmplx,sublevel = T)
  
  # compute the PD of the filtration
  D <- filtrationDiag(filtration = fltn,
                      maxdimension = maxHomDim,
                      library = 'Dionysus',
                      location = T,
                      diagLimit =1.01)$diagram
return(D)
}

# BODY
nNode <- 100 # number of nodes

# sample from a Dirichlet distribution
probs <- sample_dirichlet(nNode,alpha=c(1.5,1.5,1.5)) 

# generate a graph from the random dot product model
g <- sample_dot_product(probs)

vecDim <- 10 # dimension of vector representation
V0 <- V1 <- V2 <- matrix(NA,nrow = 4,ncol = vecDim)
rownames(V0) <- rownames(V1) <- rownames(V2) <- c("degree","closeness","betweenness","custom")

for (nodeActFunction in rownames(V0)){
  # compute PD of g using the given node activation function
  if (nodeActFunction!="custom")
    D <- graphToPD(g,maxSimplxDim = 2,maxHomDim = 2,nodeActFunction)
  else{
    nodeAct <- colSums(-probs*log2(probs)) # custom node activation values
    D <- graphToPD(g,maxSimplxDim = 2,maxHomDim = 2,nodeActFunction = "custom",nodeAct=nodeAct)
  }
  # compute VAB (vector of averaged Bettis) of D for homological dimensions 0,1 and 2
  V0[nodeActFunction,] <- computeVAB(D,homDim = 0,scaleSeq = seq(0,1,length.out=vecDim+1))
  V1[nodeActFunction,] <- computeVAB(D,homDim = 1,scaleSeq = seq(0,1,length.out=vecDim+1))
  V2[nodeActFunction,] <- computeVAB(D,homDim = 2,scaleSeq = seq(0,1,length.out=vecDim+1))
}
print(V0,digits = 3)
print(V1,digits = 3)
print(V2,digits = 3)










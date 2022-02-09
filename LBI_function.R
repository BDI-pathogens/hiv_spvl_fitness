library(ape)

lbi_o<-function(tree,l){
  treetest<-reorder.phylo(tree,order="postorder")
  ntips<-length(treetest$tip.label)
  Inodedown<-rep(0,ntips+treetest$Nnode)
  Inodedownleft<-rep(0,ntips+treetest$Nnode)
  Inodedownright<-rep(0,ntips+treetest$Nnode)
  left_done<-rep(0,ntips+treetest$Nnode)
  for(i in 1:length(treetest$edge[,1])){
    Inodedown[treetest$edge[i,1]]<-Inodedown[treetest$edge[i,1]]+Inodedown[treetest$edge[i,2]]*exp(-treetest$edge.length[i]/l)+(1-exp(-treetest$edge.length[i]/l))*l
    if(left_done[treetest$edge[i,1]]==1){
      Inodedownright[treetest$edge[i,1]]<-Inodedown[treetest$edge[i,2]]*exp(-treetest$edge.length[i]/l)+(1-exp(-treetest$edge.length[i]/l))*l    
    } else {
      Inodedownleft[treetest$edge[i,1]]<-Inodedown[treetest$edge[i,2]]*exp(-treetest$edge.length[i]/l)+(1-exp(-treetest$edge.length[i]/l))*l    
      left_done[treetest$edge[i,1]]<-1
    }
  }
  fitness_lbi<-rep(NA,ntips)
  for(i in 1:ntips){
    fitness<-0
    nodes2tip<-nodepath(treetest,from=ntips+1,to=i)
    nnodes2tip<-length(nodes2tip)
    for(n in 2:nnodes2tip){
      children_indices<-which(treetest$edge[,1]==nodes2tip[n-1])
      children<-treetest$edge[children_indices,2]
      if(children[1]==nodes2tip[n]){
        fitness<-(fitness+Inodedownright[nodes2tip[n-1]])*exp(-treetest$edge.length[children_indices[1]]/l)+(1-exp(-treetest$edge.length[children_indices[1]]/l))*l
      } else {
        fitness<-(fitness+Inodedownleft[nodes2tip[n-1]])*exp(-treetest$edge.length[children_indices[2]]/l)+(1-exp(-treetest$edge.length[children_indices[2]]/l))*l
      }
    }
    fitness_lbi[i]<-fitness
  }
  return(fitness_lbi)
}

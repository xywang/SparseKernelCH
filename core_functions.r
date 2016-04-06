# Core functions of SparseKernekHC =========================================

lw_update_rule = function(D,k,l,method,c)
{
  #Implement the LW updated formula
  #INPUT
  #D, mat (nxn): distances matrix 
  #k, l, scalars: cluster ids to be merged
  #method, string: name of the method
  #c, vec (nx1): number of points in each cluster
  #OUTPUT
  #D, mat ((n-1)x(n-1): updated distances matrix
  #c, vec ((n-1)x1): update number of points in each cluster
  list_id=rownames(D)
  k=which(list_id==k)
  l=which(list_id==l)
  Dk=D[k,]
  Dl=D[l,]
  ck=c[k]
  cl=c[l]
  switch(method, # like LW 
         "single"={w=c(0.5,0.5,0,-0.5)},
         "complete"={w=c(0.5,0.5,0,0.5)},
         "average"={w=c(ck/(ck+cl),cl/(ck+cl),0,0)},
         "mcquitty"={w=c(0.5,0.5,0,0)},
         "median"={w=c(0.5,0.5,-0.25,0)},
         "centroid"={w=c(ck/(ck+cl),cl/(ck+cl),-(ck*cl)/(ck+cl)^2,0)}
  )
  if (method=="ward")
  {
    Dn=((c+ck)*Dk+(c+cl)*Dl-c*D[k,l])/(ck+cl+c)
  }
  else
  {
    Dn=w[1]*Dk+w[2]*Dl+w[3]*D[k,l]+w[4]*abs(Dk-Dl)   
  }
  D=D[-l,-l]
  Dn=Dn[-l]
  if (length(Dn)==1)
  {
    D[k]=Dn
  }
  else
  {
    D[k,]=Dn
    D[,k]=Dn
  }
  c[k]=c[k]+c[l]
  c=as.matrix(c[-l],ncol=1)
  rownames(c)=list_id[-l]
  return(list(D=D,c=c))
}

lw_merge_rule = function(D)
{
  #Implement the AHC merging rule
  #INPUT
  #D, mat (nxn): distances matrix
  #OUTPUT
  #cluster ids to merge
  #the distance measure (height in the dendrogram)
  list_id=as.numeric(rownames(D))
  D=D+diag(rep(Inf,nrow(D))) # inf on the diag for searching the min outside diag
  m=min(D) # search for the minimum
  min_ind=which(D==m,arr.ind = TRUE) # search for the index
  return(c(sort(list_id[min_ind[1,]]),m))
}

ahc = function(D,method)
{
  #Implement the AHC algorithm (tree building)
  #INPUT
  #D, mat (nxn): distances matrix
  #method, string: name of the method
  #OUTPUT
  #dend, mat ((n-1)x3): two clusters ids that are merged, height
  list_id=as.numeric(rownames(D))
  n_ind=length(list_id) # number of data points
  dend=matrix(0,ncol=3,nrow=n_ind-1) #initialisation of dendogram
  c=matrix(1,nrow=n_ind,ncol=1) #initialisation of the size of the clusters
  rownames(c)=list_id
  cpt=1
  while (cpt<(n_ind)) # n-1 iterations
  {
    cl_to_merge=lw_merge_rule(D) # id of grouped clusters
    res=lw_update_rule(D,cl_to_merge[1],cl_to_merge[2],method,c) # matrix distance
    D=res$D;
    c=res$c;
    dend[cpt,]=c(cl_to_merge) # the id of grouped clusters, the ultrametric distance and the size of grouped cluster are stored
    cpt=cpt+1
  }
  return(dend)
}

#[150215] KERNELISATION DE LW ET AHC--------------------------------------------------------------------------------------------------------------------------------------------------------

repvec = function(v,d,n)
{
  #R equivalent of repmat (matlab)
  #INPUT 
  #v, vec (nx1): a vector to replicate
  #d, scalar (1 or 2): if 1 replicate rows, if 2 replicate columns
  #n, scalar: nb of times to replicate
  #OUTPUT
  #mat (nxd or dxn): replicated vectors
  if (d==1)
  {
    return(matrix(v,ncol=length(v),nrow=n,byrow=TRUE))
  }
  else
  {
    return(matrix(v,nrow=length(v),ncol=n,byrow=FALSE))
  }
}

kernel_lw_update_rule = function(S,k,l,method,c)
{
  #Implement the LW updated formula using kernel matrix
  #INPUT
  #S, mat (nxn): similarities matrix 
  #k, l, scalars: cluster ids to be merged
  #method, string: name of the method
  #c, vec (nx1): number of points in each cluster
  #OUTPUT
  #S, mat ((n-1)x(n-1): updated similarities matrix
  #c, vec ((n-1)x1): update number of points in each cluster  
  list_id=rownames(S) # list of id of active clusters
  k=which(list_id==k) # map k in D
  l=which(list_id==l) # map l in D
  n=min(k,l) # choose the id of cluster
  Sk=S[k,]
  Sl=S[l,]
  ck=c[k]
  cl=c[l]
  switch(method, # function as LW formula, but with pairwise similarities
         "single"={w=c(0.5,0.5,0.5,0.5,0.5,0)},
         "complete"={w=c(0.5,0.5,-0.5,0.5,0.5,0)},
         "average"={w=c(ck/(ck+cl),cl/(ck+cl),0,ck/(ck+cl),cl/(ck+cl),0)},
         "mcquitty"={w=c(0.5,0.5,0,0.5,0.5,0)},
         "median"={w=c(0.5,0.5,0,0.25,0.25,0.5)},
         "centroid"={w=c(ck/(ck+cl),cl/(ck+cl),0,(ck/(ck+cl))^2,(cl/(ck+cl))^2,2*(ck*cl)/(ck+cl)^2)},
         "ward"={w=c(ck/(ck+cl),cl/(ck+cl),0,(ck/(ck+cl))^2,(cl/(ck+cl))^2,2*(ck*cl)/(ck+cl)^2)}
  )
  Sn=w[1]*Sk+w[2]*Sl+w[3]*abs(Sk-Sl) #formula for the similarities (1)
  Snn=w[4]*S[k,k]+w[5]*S[l,l]+w[6]*S[k,l] #formula for the norms (2)
  S=S[-l,-l] # replecement in the matrix
  Sn=Sn[-l]
  if (length(Sn)==1) # in case that only 1 element remains
  {
    S[k]=Snn
  }
  else
  {
    S[k,]=Sn
    S[,k]=Sn
    S[k,k]=Snn
  }
  c[k]=c[k]+c[l]
  c=as.matrix(c[-l],ncol=1)#maj de c
  rownames(c)=list_id[-l]#bizarerie de r
  return(list(S=S,c=c))
}

kernel_lw_merge_rule = function(S,method,c)
{
  #Implement the AHC merging rule when using kernel matrix
  #INPUT
  #S, mat (nxn): similarities matrix
  #method, string: name of the method
  #c, vec (nx1): number of points in each cluster  
  #OUTPUT
  #cluster ids to merge
  #the distance measure (height in the dendrogram)
  list_id=as.numeric(rownames(S))#liste des id actifs
  Sd=diag(S)
  switch(method,
         "single"={p=2;m=-1},
         "complete"={p=2;m=-1},
         "average"={p=2;m=-1},
         "mcquitty"={p=2;m=-1},
         "median"={p=2;m=-1},
         "centroid"={p=2;m=-1},
         "ward"={
           D=(repvec(c,1,length(c))*repvec(c,2,length(c)))/(repvec(c,1,length(c))+repvec(c,2,length(c)))
           p=2*D;m=-D
         }
  )
  S=p*S+m*(repvec(Sd,1,length(Sd))+repvec(Sd,2,length(Sd)))  
  S=S+diag(rep(-Inf,nrow(S)))#on ajoute des -Inf sur la diag pour la recherche du max hors diag
  m=max(S)#recherche de maximum
  max_ind=which(S==m,arr.ind = TRUE)    
  return(c(sort(list_id[max_ind[1,]]),m))  
}

kernel_ahc = function(S,method)
{
  #Implement the AHC algorithm (tree building) when using the kernel matrix
  #INPUT
  #S, mat (nxn): similarities matrix
  #method, string: name of the method
  #OUTPUT
  #dend, mat ((n-1)x3): two clusters ids that are merged, height
  list_id=as.numeric(rownames(S)) #list of id of active clusters
  n_ind=length(list_id) # number of individuls
  dend=matrix(0,ncol=3,nrow=n_ind-1) # initialisation of dendogram
  c=matrix(1,nrow=n_ind,ncol=1) # initialisation of the size of clusters
  rownames(c)=list_id # the same list of id
  cpt=1
  while (cpt<(n_ind))# there are n-1 iterations
  {
    cl_to_merge=kernel_lw_merge_rule(S,method,c)# search id of clusters, group them
    res=kernel_lw_update_rule(S,cl_to_merge[1],cl_to_merge[2],method,c)#maj matrice distance
    S=res$S;
    c=res$c;
    dend[cpt,]=c(cl_to_merge)# the id of grouped clusters, the ultrametric distance and the size of the grouped cluster are stored
    cpt=cpt+1
  }
  return(dend)
}

P_matrices = function(X,center,scale,normalize)
{
  #Center and/or scale and/or normalize similarities and distances matrices 
  #INPUT
  #X, mat (nxp): a feature matrix
  #center (TRUE or FALSE): if TRUE center wrt the mean vector
  #scale (TRUE or FALSE): if TRUE scale each feature wrt standard deviation
  #normalize (TRUE or FALSE): if TRUE normalize the vectors lengths to 1
  #OUTPUT
  #S, mat (nxn): dot products matrix
  #D, mat (nxn): euclidean distances matrix
  X=scale(X,center = center,scale = scale)
  S=X%*%t(X)
  if (normalize)
  {
    S=S/sqrt((repvec(diag(S),1,nrow(S))*(repvec(diag(S),2,nrow(S)))))#cosine normalization
  }
  D=repvec(diag(S),1,nrow(S))+repvec(diag(S),2,nrow(S))-2*S
  rownames(D)=1:nrow(D)
  colnames(D)=1:nrow(D)
  rownames(S)=1:nrow(S)
  colnames(S)=1:nrow(S)
  return(list(S=S,D=D))
}

p_rbf_mat = function(X,center,scale,normalize,r) ####### added by xywang 15 jan 2016 ######
{
  # output S and D given a specified sigma in rbf kernel
  X=scale(X,center = center,scale = scale)
  kernel=rbfdot(r)
  S=kernelMatrix(kernel,X)
  if (normalize)
  {
    S=S/sqrt((repvec(diag(S),1,nrow(S))*(repvec(diag(S),2,nrow(S)))))#cosine normalization
  }
  D=repvec(diag(S),1,nrow(S))+repvec(diag(S),2,nrow(S))-2*S
  rownames(D)=1:nrow(D)
  colnames(D)=1:nrow(D)
  rownames(S)=1:nrow(S)
  colnames(S)=1:nrow(S)
  return(list(S=S,D=D))
}



# COPHENETIC MATRICES --------------------------------------------------------------------------------------------------------------------------------

cophenetic_matrices = function (D)
  #Implement the computation of a cophenetic matrix
  #INPUT
  #D, mat (nx3): a dendogram provided by ahc or kernel_ahc
  #OUTPUT
  #C, mat (nxn): a cophenetic matrix
{
  n=nrow(D)+1
  I=matrix(0,nrow=n,ncol=n)
  I=diag(rep(1,n))
  C=matrix(0,nrow=n,ncol=n)
  for (i in 1:(n-1))
  {
    row=D[i,1]
    col=D[i,2]
    val=D[i,3]
    row_set=which(I[row,]==1)
    col_set=which(I[col,]==1)
    C[row_set,col_set]=val
    C[col_set,row_set]=val
    I[row,col_set]=1    
  }
  return(C)
}

cophenetic_correlation = function(C1,C2)
{
  #Implement the computation of a cophenetic measure
  #INPUT
  #C1, C2 mat (nxn): two cophenetic matrices
  #OUTPUT
  #scalar: the cophenetic measure
  return(cor(C1[upper.tri(C1)],C2[upper.tri(C2)]))
}


convert_dend = function(dend,list_id,k)
  # convert dendrogram to plottable format
  # dendrogram = output of ahc() or kernel_ahc
  # list_id = as.numeric(rownames(D)) or list_id = as.numeric(rownames(S))
  # k = user input, num of expected clusters
  # res will be used in "plot(dataset, col=res)"
  {
    n_ind = length(list_id)
    res <- matrix(c(list_id),nrow=n_ind,ncol=1)
    rownames(res)<-list_id
    dend = dend[,-3]
    n = nrow(dend)+1
    I = diag(rep(1,n))
    for (i in 1:(n-1))
    {row = dend[i,1]
     col = dend[i,2]
     col_set = which(I[col,]==1)
     I[row,col_set] =1
     res[col_set] <- row
     if (k == length(unique(res))) break
    }
    return(res)
  } 


get_ker_D = function(s_dot)
  # s_dot is the kernal similarity matrix, 
  # if s_dot is a dot product, then ker_D is the euclidean distance
{
  dia_s_dot <- diag(s_dot)
  mat_1 <- repvec(dia_s_dot,1,nrow(s_dot))
  mat_2 <- repvec(dia_s_dot,2,nrow(s_dot))
  mat_3 <- 2*s_dot  
  dist_mat <- mat_1 + mat_2 - mat_3
  return(sqrt(dist_mat))
}

# ----------------------------
ARI_kerAHC = function(x,m)
{ 
  dend = kernel_ahc(x,m)
  list_id = as.numeric(rownames(x))
  pre_tags = convert_dend(dend,list_id,15)
  #tags <- factor(pre_tags, labels = c(1:length(unique(pre_tags))))
  #plot(data, col=tags)
  ARI = adjustedRandIndex(label, pre_tags)
  return(ARI)
}


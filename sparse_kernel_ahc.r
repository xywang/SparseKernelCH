# Implementation of SparseKernelHC
# similarity matrix is formatted in sparse way: i.e
# each value is represented by a triple (ind.row, ind.vol, val)

######################################################## 

sp_lw_update_rule = function(Di,Dj,Dx,k,l,method,c)
{
  #Implement the LW updated formula
  #INPUT
  #Di, Dj , Dx, vecs (nx1): distance matrix in sparse triple format
  #k, l, scalars: cluster ids to be merged
  #method, string: name of the method
  #c, vec (nx1): number of points in each cluster
  #OUTPUT
  #Di, Dj, Dx; vecs ((n-1)x1): updated distances matrix
  #c; vec ((n-1)x1): update number of points in each cluster
  list_id=rownames(c)
  id.ik=which(Di==k)#k is involved as a row
  id.jk=which(Dj==k)#k is involved as a column

  Dk=cbind(c(Dj[id.ik],Di[id.jk]),c(Dx[id.ik],Dx[id.jk]))#distance row of k
  id.il=which(Di==l)#l is involved as a row
  id.jl=which(Dj==l)#l is involved as a column

  max.id=max(as.integer(list_id))
  Dl=cbind(c(Dj[id.il],Di[id.jl]),c(Dx[id.il],Dx[id.jl]))#distance row of k
  Dkl=Dx[intersect(id.ik,id.jl)]#Dkl value to be used afterwards
  all.id=union(union(id.ik,id.jk),union(id.il,id.jl))#all id involving k and l
  Dk=spMatrix(nrow=1, ncol=max.id+1, i=rep(1,nrow(Dk)), j=as.integer(Dk[,1])+1, x=Dk[,2])
  Dl=spMatrix(nrow=1, ncol=max.id+1, i=rep(1,nrow(Dl)), j=as.integer(Dl[,1])+1, x=Dl[,2])
  kk=which(list_id==k)
  ll=which(list_id==l)
  ck=c[kk]
  cl=c[ll]
  switch(method,
         "single"={w=c(0.5,0.5,0,-0.5)},
         "complete"={w=c(0.5,0.5,0,0.5)},
         "average"={w=c(ck/(ck+cl),cl/(ck+cl),0,0)},
         "mcquitty"={w=c(0.5,0.5,0,0)},
         "median"={w=c(0.5,0.5,-0.25,0)},
         "centroid"={w=c(ck/(ck+cl),cl/(ck+cl),-(ck*cl)/(ck+cl)^2,0)}
  )
  #compute the new row
  if (method=="ward")
  {
    cck=spMatrix(nrow=1, ncol=max.id+1, i=rep(1,nrow(c)), j=as.integer(list_id)+1, x=(c+ck)/(c+ck+cl))
    ccl=spMatrix(nrow=1, ncol=max.id+1, i=rep(1,nrow(c)), j=as.integer(list_id)+1, x=(c+cl)/(c+ck+cl))
    cckcl=spMatrix(nrow=1, ncol=max.id+1, i=rep(1,nrow(c)), j=as.integer(list_id)+1, x=-c/(c+ck+cl))
    Dn=cck*Dk+ccl*Dl+cckcl*Dkl
  }
  else#cas des autres méthodes
  {
    Dn=w[1]*Dk+w[2]*Dl+w[3]*Dkl+w[4]*abs(Dk-Dl)   
  }
  #remove from previous values the obsolete row
  Di=Di[-all.id]
  Dj=Dj[-all.id]
  Dx=Dx[-all.id]
  #format the new row in sparse format
  Dn=as(Dn,"dgTMatrix")
  Dnj=Dn@j
  Dnx=Dn@x
  del.id.k=which(Dnj==k)#remove k in the new row
  del.id.l=which(Dnj==l)#remove l in the new row
  del.non.activ.id=which(!is.element(Dnj,list_id))#remove sparse element in the new row
  Dnj=Dnj[-c(del.id.k,del.id.l,del.non.activ.id)]
  Dnx=Dnx[-c(del.id.k,del.id.l,del.non.activ.id)]
  Dni=rep(k,length(Dnj))
  find.switch=which(Dnj<k)
  if (length(find.switch)>0)
  {
    Dni[find.switch]=Dnj[find.switch]
    Dnj[find.switch]=k
  }
  if (length(Dn)==1)
  {
    Di=k
    Dj=k
    Dx=Dn
  }
  else
  {
    Di=c(Di,Dni)
    Dj=c(Dj,Dnj)
    Dx=c(Dx,Dnx)
  }
  c[kk]=c[kk]+c[ll]
  c=as.matrix(c[-ll],ncol=1)
  rownames(c)=list_id[-ll]
  return(list(Di=Di,Dj=Dj,Dx=Dx,c=c))
}

sp_ahc = function(Di,Dj,Dx,nb_ind,method)
{
  #Implement the AHC algorithm (tree building)
  #INPUT
  #D, mat (nxn): distances matrix
  #method, string: name of the method
  #OUTPUT
  #dend, mat ((n-1)x3): two clusters ids that are merged, height
  library("Matrix")
  list_id=0:(nb_ind-1)
  c=matrix(1,nrow=nb_ind,ncol=1)
  rownames(c)=list_id
  cpt=1
  while (cpt<nb_ind)
  {
    ind.d.min=which.min(Dx)
    height=Dx[ind.d.min]
    cl_to_merge=c(min(Di[ind.d.min],Dj[ind.d.min]),max(Di[ind.d.min],Dj[ind.d.min]))
    res=sp_lw_update_rule(Di,Dj,Dx,cl_to_merge[1],cl_to_merge[2],method,c)
    Di=res$Di;
    Dj=res$Dj;
    Dx=res$Dx;
    c=res$c;
    dend[cpt,]=c(cl_to_merge+1,height)
    cpt=cpt+1
  }
  return(cbind(dend,as.numeric(object.size(Di))+as.numeric(object.size(Dj))+as.numeric(object.size(Dx))))
}

#--------------------------------------------------------------

sp_kernel_lw_update_rule = function(Si,Sj,Sx,Sd,k,l,method,c)
{
  #Implement the LW updated formula
  #INPUT
  #Si, Sj , Sx, Sd; vecs (nx1): distance matrix in sparse triple format + diagonal
  #k, l, scalars: cluster ids to be merged
  #method, string: name of the method
  #c, vec (nx1): number of points in each cluster
  #OUTPUT
  #Si, Sj, Sx, Sd; vecs ((n-1)x1): updated distances matrix
  #c; vec ((n-1)x1): update number of points in each cluster
  list_id=rownames(c)
  id.ik=which(Si==k)#k is involved as a row
  id.jk=which(Sj==k)#k is involved as a column

  Sk=cbind(c(Sj[id.ik],Si[id.jk]),c(Sx[id.ik],Sx[id.jk]))#similarity row of k
  id.il=which(Si==l)#l is involved as a row
  id.jl=which(Sj==l)#l is involved as a column

  max.id=max(as.integer(list_id))
  Sl=cbind(c(Sj[id.il],Si[id.jl]),c(Sx[id.il],Sx[id.jl]))#similarity row of k
  Skl=Sx[intersect(id.ik,id.jl)]#Skl value to be used afterwards
  all.id=union(union(id.ik,id.jk),union(id.il,id.jl))#all id involving k and l
  Sk=spMatrix(nrow=1, ncol=max.id+1, i=rep(1,nrow(Sk)), j=as.integer(Sk[,1])+1, x=Sk[,2])
  Sl=spMatrix(nrow=1, ncol=max.id+1, i=rep(1,nrow(Sl)), j=as.integer(Sl[,1])+1, x=Sl[,2])
  kk=which(list_id==k)
  ll=which(list_id==l)
  ck=c[kk]
  cl=c[ll]
  Skk=Sd[kk]#Skk value to be used afterwards
  Sll=Sd[ll]#Sll value to be used afterwards
  #===
  switch(method,
         "single"={w=c(0.5,0.5,0.5,0.5,0.5,0)},
         "complete"={w=c(0.5,0.5,-0.5,0.5,0.5,0)},
         "average"={w=c(ck/(ck+cl),cl/(ck+cl),0,ck/(ck+cl),cl/(ck+cl),0)},
         "mcquitty"={w=c(0.5,0.5,0,0.5,0.5,0)},
         "median"={w=c(0.5,0.5,0,0.25,0.25,0.5)},
         "centroid"={w=c(ck/(ck+cl),cl/(ck+cl),0,(ck/(ck+cl))^2,(cl/(ck+cl))^2,2*(ck*cl)/(ck+cl)^2)},
         "ward"={w=c(ck/(ck+cl),cl/(ck+cl),0,(ck/(ck+cl))^2,(cl/(ck+cl))^2,2*(ck*cl)/(ck+cl)^2)}
  )
  Sn=w[1]*Sk+w[2]*Sl+w[3]*abs(Sk-Sl)
  Snn=w[4]*Skk+w[5]*Sll+w[6]*Skl
  #remove from previous values the obsolete row
  Si=Si[-all.id]
  Sj=Sj[-all.id]
  Sx=Sx[-all.id]
  #format the new row in sparse format
  Sn=as(Sn,"dgTMatrix")
  Snj=Sn@j
  Snx=Sn@x
  del.id.k=which(Snj==k)
  del.id.l=which(Snj==l)
  del.non.activ.id=which(!is.element(Snj,list_id))
  Snj=Snj[-c(del.id.k,del.id.l,del.non.activ.id)]
  Snx=Snx[-c(del.id.k,del.id.l,del.non.activ.id)]
  Sni=rep(k,length(Snj))
  find.switch=which(Snj<k)
  if (length(find.switch)>0)
  {
    Sni[find.switch]=Snj[find.switch]
    Snj[find.switch]=k
  }
  if (length(Sn)==1)
  {
    Si=k
    Sj=k
    Sx=Sn
    Sd=Snn
  }
  else
  {
    Si=c(Si,Sni)
    Sj=c(Sj,Snj)
    Sx=c(Sx,Snx)
    Sd[kk]=Snn
    Sd=Sd[-ll]
  }
  c[kk]=c[kk]+c[ll]
  c=as.matrix(c[-ll],ncol=1)
  rownames(c)=list_id[-ll]
  return(list(Si=Si,Sj=Sj,Sx=Sx,Sd=Sd,c=c))
}

sp_kernel_lw_merge_rule = function(Si,Sj,Sx,Sd,method,c)
{
  #Implement the AHC merging rule when using kernel matrix
  #INPUT
  #S, mat (nxn): similarities matrix
  #method, string: name of the method
  #c, vec (nx1): number of points in each cluster  
  #OUTPUT
  #cluster ids to merge
  #the distance measure (height in the dendrogram)
  list_id=rownames(c)#liste des id actifs
  max.id=max(as.numeric(list_id))
  switch(method,
         "single"={p=2;m=-1},
         "complete"={p=2;m=-1},
         "average"={p=2;m=-1},
         "mcquitty"={p=2;m=-1},
         "median"={p=2;m=-1},
         "centroid"={p=2;m=-1},
         "ward"={
           cc=spMatrix(nrow=1, ncol=max.id+1, i=rep(1,length(list_id)), j=as.integer(list_id)+1, x=c)
           p=2*cc[Si+1]*cc[Sj+1]/(cc[Si+1]+cc[Sj+1])
           m=-p/2
         }
  )
  SSd=spMatrix(nrow=1, ncol=max.id+1, i=rep(1,length(list_id)), j=as.integer(list_id)+1, x=Sd)
  S=p*Sx+m*(SSd[Si+1]+SSd[Sj+1])
  max.ind=which.max(S)
  m=S[max.ind]
  return(c(Si[max.ind],Sj[max.ind],m))  
}

sp_kernel_ahc = function(Si,Sj,Sx,Sd,method)
{
  #Implement the AHC algorithm (tree building) when using the kernel matrix
  #INPUT
  #S, mat (nxn): similarities matrix
  #method, string: name of the method
  #OUTPUT
  #dend, mat ((n-1)x4): two clusters ids that are merged, height, size of S
  library("Matrix")
  nb_ind=length(Sd)
  list_id=0:(nb_ind-1)
  dend=matrix(0,ncol=4,nrow=nb_ind-1)
  c=matrix(1,nrow=nb_ind,ncol=1)
  rownames(c)=list_id
  cpt=1
  while (cpt<(nb_ind) && length(Sx)>0)
  {
    cl_to_merge=sp_kernel_lw_merge_rule(Si,Sj,Sx,Sd,method,c)
    res=sp_kernel_lw_update_rule(Si,Sj,Sx,Sd,cl_to_merge[1],cl_to_merge[2],method,c)
    Si=res$Si;
    Sj=res$Sj;
    Sx=res$Sx;
    Sd=res$Sd;
    c=res$c;
    dend[cpt,]=c(cl_to_merge[1]+1,cl_to_merge[2]+1,cl_to_merge[3],as.numeric(object.size(Si))+as.numeric(object.size(Sj))+as.numeric(object.size(Sx))+as.numeric(object.size(Sd)))
    cpt=cpt+1
  }
  return(dend)
}

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

sp_prox_matrices = function(X,center,scale,normalize,kernel)
{
  #Center and/or scale and/or normalize and/or project (kernel) similarities and distances matrices 
  #INPUT
  #X, mat (nxp): a feature matrix
  #center (TRUE or FALSE): if TRUE center wrt the mean vector
  #scale (TRUE or FALSE): if TRUE scale each feature wrt standard deviation
  #normalize (TRUE or FALSE): if TRUE normalize the vectors lengths to 1
  #kernel (linear,quadratic,gaussian): 3 different types of kernel
  #OUTPUT
  #S, mat (nxn): dot products matrix
  #D, mat (nxn): euclidean distances matrix
  X=scale(X,center = center,scale = scale)
  library("kernlab")
  switch(kernel,
         "linear"={kernel=vanilladot()},
         "gaussian"={kernel=rbfdot(sigma=2)},
         "quadratic"={kernel=polydot(degree=2,scale=1,offset=1)},
         "sigmoid"={kernel=tanhdot(scale=1,offset=1)}
  )
  S=kernelMatrix(kernel,X)
  if (normalize)
  {
    S=S/sqrt((repvec(diag(S),1,nrow(S))*(repvec(diag(S),2,nrow(S)))))#cosine normalization
  }
  Sd=matrix(diag(S),nrow=nrow(S))  
  D=repvec(Sd,1,nrow(S))+repvec(Sd,2,nrow(S))-2*S
  library("Matrix")
  S=triu(as(S,"dgTMatrix"),k=1)
  D=triu(as(D,"dgTMatrix"),k=1)
  return(list(S=S,Sd=Sd,D=D))
}

#[150213] COPHENETIC MATRICES --------------------------------------------------------------------------------------------------------------------------------

cophenetic_matrices = function(D)
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


cut_dend = function(dend,k)
{
  #Cut the dendogram in order to have an assignment matrix that represents a flat partition
  #INPUT
  #dend, mat (nx3) : an ouput given by ahc or kernel_ahc representing a dendogram
  #k, sca : the number of clusters
  #OUTPUT
  #res, vec (nx1) : an assignment matrix
  nb_clus=max(dend[,2])
  res=seq(1,nb_clus)
  res[res==dend[1,2]]=dend[1,1]
  nb_clus=nb_clus-1
  i=2
  while (nb_clus>k)
  {
    res[res==dend[i,2]]=dend[i,1]
    nb_clus=nb_clus-1
    i=i+1
  }
  list_id_clus=sort(unique(res))
  for (i in 1:length(list_id_clus))
  {
    res[res==list_id_clus[i]]=i
  }
  return(res)
}

# added for using parVapply() ------------------------------

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

# parVapply() w/o thresholds ----------------
# ARI_SPkerAHC = function(m)
# { 
#   dend = sp_kernel_ahc(Si,Sj,Sx,Sd,m)
#   pre_tags = convert_dend(dend,list_id,15)
#   ARI = adjustedRandIndex(label, pre_tags)
#   return(ARI)
# }

# parVapply() with thresholds -------------

ARI_SPkerAHC = function(si,sj,sx,m)
{
  dend = sp_kernel_ahc(si,sj,sx,Sd,m)
  pre_tags = convert_dend(dend,list_id,2)
  ARI = adjustedRandIndex(c_y, pre_tags)
  return(ARI)
}

with_thre = function(s)
{
  filt.id = which(Sx>=s)
  si = Si[filt.id]
  sj = Sj[filt.id]
  sx = Sx[filt.id]
  get_ARI_SPkerAHC <- sapply(methL, function(m) ARI_SPkerAHC(si,sj,sx,m))
  return(get_ARI_SPkerAHC)
}

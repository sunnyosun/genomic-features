# Finding Transcription Termination Sites
# Unidirectional genes
# algorithm 2
# Xiaoji SUN
# 12/11/2013

# Read wiggle files ######################################################################################
t1 = read.table('mRNA_treplication.coverage.F.rpkm.wig', skip=1, fill=T, stringsAsFactors=F)
t2 = read.table('mRNA_treplication.coverage.Rpos.rpkm.wig', skip=1, fill=T, stringsAsFactors=F)

index1 = which(t1[,3]!="")
index2 = which(t2[,3]!="")

data1 = list()
for (i in 1:16)
{
  data1[[i]] = t1[(index1[i]+1):(index1[i+1]-1),1:2]
}

data2 = list()
for (i in 1:16)
{
  data2[[i]] = t2[(index2[i]+1):(index2[i+1]-1),1:2]
}

# read in Red1 signals
red1 = read.table('AH119B_afterfiting_norm_all.wig', skip=1, fill=T, stringsAsFactors=F)
index_red1 = which(red1[,3]!="")
red1_data = list()
for (i in 1:16)
{
  red1_data[[i]] = red1[(index_red1[i]+1):(index_red1[i+1]-1),1:2]
}

# read in MNase data
mnase = read.table('MNase_3h_sgrp_P10_treat_afterfiting_all.wig', skip=1, fill=T, stringsAsFactors=F)
index_mnase = which(mnase[,3]!="")
mnase_data = list()
for (i in 1:16)
{
  mnase_data[[i]] = mnase[(index_mnase[i]+1):(index_mnase[i+1]-1),1:2]
}

# Read gene gff files #####################################################################################
genes = read.table('sk1_sgrp_genes.gff',stringsAsFactors=F)
colnames(genes) = c('chr', 'start', 'end', 'strand', 'id')

# Watson genes
w = matrix(NA,ncol=5,nrow=1)
colnames(w) = c('chr', 'start', 'end', 'strand', 'id')
for (i in 1:nrow(genes))
{
  if (genes[i,'strand']=='+' & genes[i+1, 'strand']=='+' & genes[i,'chr']==genes[i+1,'chr'])
  {
    w = rbind(w, genes[i:(i+1),])
  }
}
w = w[2:nrow(w),]
rownames(w) = c(1:nrow(w))
write.table(w[,1:5], 'sk1_sgrp_w.txt', quote=F, sep='\t', col.names=F, row.names=F)

# Crick genes
c = matrix(NA,ncol=5,nrow=1)
colnames(c) = c('chr', 'start', 'end', 'strand', 'id')
for (i in 1:nrow(genes))
{
  if (genes[i,'strand']=='-' & genes[i+1, 'strand']=='-' & genes[i,'chr']==genes[i+1,'chr'])
  {
    c = rbind(c, genes[i:(i+1),])
  }
}
c = c[2:nrow(c),]
rownames(c) = c(1:nrow(c))
write.table(c[,1:5], 'sk1_sgrp_c.txt', quote=F, sep='\t', col.names=F, row.names=F)

# read in the w genes
w = read.table('sk1_sgrp_w.txt',stringsAsFactors=F)
colnames(w) = c('chr', 'start', 'end', 'strand', 'id')

# read in the c genes
c = read.table('sk1_sgrp_c.txt',stringsAsFactors=F)
colnames(c) = c('chr', 'start', 'end', 'strand', 'id')

# create gene pairs
# chr index
for (i in 1:16)
{
  if (i<=9)
  {
    chr = paste('chr0',i,sep='')
    w[which(w[,1]==chr),1]=i
  }
  if (i>9)
  {
    chr = paste('chr',i,sep='')
    w[which(w[,1]==chr),1]=i
  }
}

# chr index
for (i in 1:16)
{
  if (i<=9)
  {
    chr = paste('chr0',i,sep='')
    c[which(c[,1]==chr),1]=i
  }
  if (i>9)
  {
    chr = paste('chr',i,sep='')
    c[which(c[,1]==chr),1]=i
  }
}

# w
w = w[,1:5]
for (i in c(1:(nrow(w)/2))*2-1)
{
  print (i)
  c = as.numeric(w[i,'chr'])
  end1 = as.numeric(w[i,'end'])
  end2 = as.numeric(w[i+1, 'start'])  
  chr_data = data1[[c]]
  conv_data = chr_data[which(as.numeric(chr_data[,1])>=end1 & as.numeric(chr_data[,1])<end2),]
  if (nrow(conv_data)>1)
  {    
    for (j in 1:(nrow(conv_data)-1))
    {    
      if (as.numeric(conv_data[nrow(conv_data),1])>2)
      {
        w[i,6] = conv_data[nrow(conv_data),1]
      }
      if (as.numeric(conv_data[j+1,2])<2 & as.numeric(conv_data[j,2])>2)
      {
        w[i,6] = conv_data[j,1]
      }
      if (as.numeric(conv_data[j+1,1])-as.numeric(conv_data[j,1])>10)
      {
        w[i,6] = conv_data[j,1]
        break
      }
    }
  }
}

# c
c=c[,1:5]
for (i in c(1:(nrow(c)/2))*2-1)
{
  print (i)
  r = as.numeric(c[i,'chr'])
  end1 = as.numeric(c[i,'end'])
  end2 = as.numeric(c[i+1, 'start'])  
  chr_data = data2[[r]]
  conv_data = chr_data[which(as.numeric(chr_data[,1])<=end2 & as.numeric(chr_data[,1])>end1),]
  if (nrow(conv_data)>1)
  {
    for (j in (nrow(conv_data)):2)
    {    
      if (as.numeric(conv_data[1,1])>2)
      {
        c[i,6] = conv_data[1,1]
      }
      if (as.numeric(conv_data[j-1,2])<2 & as.numeric(conv_data[j,2])>2)
      {
        c[i,6] = conv_data[j,1]
      }
      if (as.numeric(conv_data[j,1])-as.numeric(conv_data[j-1,1])>10)
      {
        c[i,6] = conv_data[j,1]
        break
      }           
    }
  }
}

# create pairs
w_p = data.frame(matrix(NA, nrow=nrow(w)/2, ncol=10))
for (i in seq(1,nrow(w),by=2))
{
  w_p[(i+1)/2,1:6] = w[i,1:6]
  w_p[(i+1)/2,7:10] = w[i+1,2:5]
}
colnames(w_p)=c('chr','start1','end1','strand1','id1','TransEnd','start2','end2','strand2','id2')

c_p = data.frame(matrix(NA, nrow=nrow(c)/2, ncol=10))
for (i in seq(1,nrow(c),by=2))
{
  c_p[(i+1)/2,1:6] = c[i,1:6]
  c_p[(i+1)/2,7:10] = c[i+1,2:5]
}
colnames(c_p)=c('chr','start1','end1','strand1','id1','TransEnd','start2','end2','strand2','id2')

# remove NAs
w2 = w_p[which(!is.na(w_p[,6])),]
write.table(w2,'TransEnds_w_filtered.txt',quote=F, sep='\t', col.names=T,row.names=F)
c2 = c_p[which(!is.na(c_p[,6])),]
write.table(c2,'TransEnds_c_filtered.txt',quote=F, sep='\t', col.names=T,row.names=F)

###
out=w2[,c(1,3,6)]
for (i in 1:nrow(out))
{
  if (as.numeric(out[i,1])<=9)
  {
    out[i,1]=paste('chr0',out[i,1],sep='')
  }
  else
  {
    out[i,1]=paste('chr',out[i,1],sep='')
  }
}
write.table(out,'trans_end_ww.bed',sep='\t',quote=F,col.names=F,row.names=F)

out=c2[,c(1,6,7)]
for (i in 1:nrow(out))
{
  if (as.numeric(out[i,1])<=9)
  {
    out[i,1]=paste('chr0',out[i,1],sep='')
  }
  else
  {
    out[i,1]=paste('chr',out[i,1],sep='')
  }
}
write.table(out,'trans_end_cc.bed',sep='\t',quote=F,col.names=F,row.names=F)

#######################
# 500bp either side of transEnds
w2 = read.table('TransEnds_w_filtered.txt', header=T, stringsAsFactors=F)
c2 = read.table('TransEnds_c_filtered.txt', header=T, stringsAsFactors=F)

n = 500
w_matrix=data.frame(matrix(0,nrow=nrow(w2),ncol=2*n))
for (i in 1:nrow(w2))
{
  print (i)
  r = as.numeric(w2[i,'chr'])
  #chr_data = red1_data[[r]]
  chr_data = mnase_data[[r]]
  start = as.numeric(w2[i,'TransEnd'])-n
  end = as.numeric(w2[i,'TransEnd'])+n
  tmp = chr_data[which(as.numeric(chr_data[,1])>=start & as.numeric(chr_data[,1])<end),]
  if (dim(tmp)[1]==2*n)
  {
    w_matrix[i,]=as.numeric(tmp[,2])
  }
}
#write.table(w_matrix,'rpkm_w_500.txt', sep='\t', quote=F, col.names=F, row.names=F)
write.table(w_matrix/46.5,'rpkm_w_mnase_500.txt', sep='\t', quote=F, col.names=F, row.names=F)

n = 500
c_matrix=data.frame(matrix(0,nrow=nrow(c2),ncol=2*n))
for (i in 1:nrow(c2))
{
  print (i)
  r = as.numeric(c2[i,'chr'])
  #chr_data = red1_data[[r]]
  chr_data = mnase_data[[r]]
  start = as.numeric(c2[i,'TransEnd'])-n
  end = as.numeric(c2[i,'TransEnd'])+n
  tmp = chr_data[which(as.numeric(chr_data[,1])>=start & as.numeric(chr_data[,1])<end),]
  if (dim(tmp)[1]==2*n)
  {
    c_matrix[i,]=as.numeric(tmp[,2])
  }
}
#write.table(c_matrix, 'rpkm_c_500.txt', sep='\t', quote=F, col.names=F, row.names=F)
write.table(c_matrix/46.5, 'rpkm_c_mnase_500.txt', sep='\t', quote=F, col.names=F, row.names=F)

# normalize
#w_matrix=read.table('rpkm_w_500.txt')
w_matrix=read.table('rpkm_w_mnase_500.txt')
w_matrix_norm = w_matrix
for (i in 1:nrow(w_matrix))
{
  if (w_matrix[i,1]!=0)
  {
    x=as.vector(t(w_matrix[i,]))
    y=(x-mean(x))/sd(x)
    w_matrix_norm[i,]=y
  }
}
#write.table(w_matrix_norm,'w_matrix_500_norm.txt',col.names=F, row.names=F, sep='\t', quote=F)
write.table(w_matrix_norm,'w_matrix_mnase_500_norm.txt',col.names=F, row.names=F, sep='\t', quote=F)

#c_matrix=read.table('rpkm_c_500.txt')
c_matrix=read.table('rpkm_c_mnase_500.txt')
c_matrix_norm = c_matrix
for (i in 1:nrow(c_matrix))
{
  if (c_matrix[i,1]!=0)
  {
    x=as.vector(t(c_matrix[i,]))
    y=(x-mean(x))/sd(x)
    c_matrix_norm[i,]=y
  }
}
#write.table(c_matrix_norm,'c_matrix_500_norm.txt',col.names=F, row.names=F, sep='\t', quote=F)
write.table(c_matrix_norm,'c_matrix_mnase_500_norm.txt',col.names=F, row.names=F, sep='\t', quote=F)



# color plot
# sort by coverage
w_matrix_norm=read.table('w_matrix_500_norm.txt')
#w_matrix_norm=read.table('rpkm_w_500.txt')
range(w_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- cbind(w_matrix_norm, rowSums(w_matrix_norm))
o1<- rev(order(m.row.sum[,1001]))
m.row.sum<- m.row.sum[o1,]
bk = unique(c(seq(0,1, length=100),seq(1,3,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:1000], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=T, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

# sort by midpoint value
c_matrix_norm=read.table('c_matrix_500_norm.txt')
#conv_matrix=read.table('nucleosome_conv_500.txt')
range(c_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- cbind(c_matrix_norm, rowSums(c_matrix_norm))
o1<- rev(order(m.row.sum[,1001]))
m.row.sum<- m.row.sum[o1,]
bk = unique(c(seq(0,1, length=100),seq(1,3,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:1000], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=T, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

# sort by midpoint value
w_matrix_norm=read.table('w_matrix_500_norm.txt')
range(w_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- w_matrix_norm
o1<- rev(order(m.row.sum[,500]))
# index_mid_w=o1
m.row.sum<- m.row.sum[o1,]
#bk = unique(c(seq(-0.1,3, length=100),seq(3,7,length=100)))
bk = unique(c(seq(0,1, length=100),seq(1,3,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:1000], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

c_matrix_norm=read.table('c_matrix_500_norm.txt')
range(c_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- c_matrix_norm
o1<- rev(order(m.row.sum[,500]))
# index_mid_c=o1
m.row.sum<- m.row.sum[o1,]
#bk = unique(c(seq(-0.1,3, length=100),seq(3,7,length=100)))
bk = unique(c(seq(0,1, length=100),seq(1,3,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:1000], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

# mnase
# normalize

# color plot
# sort by coverage
w_matrix_norm=read.table('w_matrix_mnase_500_norm.txt')
#w_matrix_norm=read.table('rpkm_w_500.txt')
range(w_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- cbind(w_matrix_norm, rowSums(w_matrix_norm))
o1<- rev(order(m.row.sum[,1001]))
m.row.sum<- m.row.sum[o1,]
bk = unique(c(seq(0,1, length=100),seq(1,3,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:1000], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=T, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

# sort by coverage
c_matrix_norm=read.table('c_matrix_mnase_500_norm.txt')
#conv_matrix=read.table('nucleosome_conv_500.txt')
range(c_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- cbind(c_matrix_norm, rowSums(c_matrix_norm))
o1<- rev(order(m.row.sum[,1001]))
m.row.sum<- m.row.sum[o1,]
bk = unique(c(seq(0,1, length=100),seq(1,3,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:1000], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=T, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

# sort by midpoint value
w_matrix_norm=read.table('w_matrix_mnase_500_norm.txt')
range(w_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- w_matrix_norm
#o1<- rev(order(m.row.sum[,500]))
m.row.sum<- m.row.sum[index_mid_w,]
#bk = unique(c(seq(-0.1,3, length=100),seq(3,7,length=100)))
bk = unique(c(seq(0,1, length=100),seq(1,3,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:1000], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

c_matrix_norm=read.table('c_matrix_mnase_500_norm.txt')
range(c_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- c_matrix_norm
#o1<- rev(order(m.row.sum[,500]))
m.row.sum<- m.row.sum[index_mid_c,]
#bk = unique(c(seq(-0.1,3, length=100),seq(3,7,length=100)))
bk = unique(c(seq(0,1, length=100),seq(1,3,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:1000], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off()


# ave
w_matrix=read.table('rpkm_w_500.txt')
c_matrix=read.table('rpkm_c_500.txt')
red1_ave_vector_w=colMeans(w_matrix)
red1_ave_vector_c=colMeans(c_matrix)
png("test.png", width=700, height = 500)
plot(c(1:1000)-500, red1_ave_vector_w, pch=16, col='blue', xlab='Mid Point of Convergent Transcript Ends', ylab='Red1', frame.plot=F)
plot(c(1:1000)-500, red1_ave_vector_c, pch=16, col='red', xlab='Mid Point of Convergent Transcript Ends', ylab='Red1', frame.plot=F)
dev.off()

w_matrix=read.table('rpkm_w_mnase_500.txt')
c_matrix=read.table('rpkm_c_mnase_500.txt')
red1_ave_vector_w=colMeans(w_matrix)
red1_ave_vector_c=colMeans(c_matrix)
png("test.png", width=700, height = 500)
plot(c(1:1000)-500, red1_ave_vector_w, pch=16, col='blue', xlab='Mid Point of Convergent Transcript Ends', ylab='MNase', frame.plot=F)
plot(c(1:1000)-500, red1_ave_vector_c, pch=16, col='red', xlab='Mid Point of Convergent Transcript Ends', ylab='MNase', frame.plot=F)
dev.off()

w_matrix=read.table('w_matrix_500_norm.txt')
c_matrix=read.table('c_matrix_500_norm.txt')
red1_ave_vector_w=colMeans(w_matrix)
red1_ave_vector_c=colMeans(c_matrix)
w_mnase=read.table('w_matrix_mnase_500_norm.txt')
c_mnase=read.table('c_matrix_mnase_500_norm.txt')
mnase_ave_vector_w=colMeans(w_mnase)
mnase_ave_vector_c=colMeans(c_mnase)
png("test.png", width=700, height = 500)
plot(c(1:1000)-500, red1_ave_vector_c, pch=16, col='red', xlab='Mid Point of Convergent Transcript Ends', ylab='Red1', frame.plot=F)
points(c(1:1000)-500, mnase_ave_vector_c, pch=16, col='blue')
dev.off()

#####
a=read.table('TransEnds_conv_filtered.txt',header=T)
b=read.table('TransEnds_w_filtered.txt', header=T)
c=read.table('TransEnds_c_filtered.txt', header=T)
d=a
colnames(d)=colnames(b)[1:6]
d=rbind(d,b[,1:6])
colnames(d)=colnames(c)[c(1,7,8,9,10,6)]
d=rbind(d,c[,c(1,7,8,9,10,6)])
write.table(d,'TransEnds_genes.txt', quote=F, sep='\t', col.names=T, row.names=F)

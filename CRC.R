setwd("C:/project/neoantigen/CRC")

###mutation
inserm<-read.delim("./raw_data/coadread_dfci_2016/data_mutations.txt")
msk<-read.delim("./raw_data/coadread_genentech/data_mutations.txt",skip = 1)
meric<-read.delim("./raw_data/coadread_mskresistance_2022/data_mutations.txt")
mskimpact<-read.delim("./raw_data/coadread_tcga/data_mutations.txt")
prv<-read.delim("./raw_data/coadread_tcga_pan_can_atlas_2018/data_mutations.txt")
tcga<-read.delim("./raw_data/coadread_tcga_pub/data_mutations.txt",skip=2)
atlas<-read.delim("./raw_data/crc_apc_impact_2020/data_mutations.txt",skip=2)
dd<-read.delim("./raw_data/crc_dd_2022/data_mutations.txt",skip=2)


used_col=c("Hugo_Symbol","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","HGVSp_Short","Tumor_Sample_Barcode","Variant_Classification")
inserm=inserm[,which(colnames(inserm) %in% used_col)]
msk=msk[,which(colnames(msk) %in% used_col)]
meric=meric[,which(colnames(meric) %in% used_col)]
mskimpact=mskimpact[,which(colnames(mskimpact) %in% used_col)]
prv=prv[,which(colnames(prv) %in% used_col)]
tcga=tcga[,which(colnames(tcga) %in% used_col)]
atlas=atlas[,which(colnames(atlas) %in% used_col)]
dd=dd[,which(colnames(dd) %in% used_col)]
mut=rbind(inserm,msk,meric,mskimpact,prv,tcga,atlas,dd)
mut=unique(mut)
idx=grep("chr",mut$Chromosome)
mut$Chromosome[-idx]=paste0("chr",mut$Chromosome[-idx])
mut_delSilent=mut[-which(mut$Variant_Classification == "Silent"),]
mut_Silent=mut[which(mut$Variant_Classification == "Silent"),]
write.csv(mut_delSilent,"./result/mut_delSilent_info.csv")
mut_judge=paste0(mut_delSilent$Hugo_Symbol,":",mut_delSilent$Chromosome,":",mut_delSilent$Start_Position,":",mut_delSilent$HGVSp_Short)

########################################################################
#------------——————-计算所有位点的人群频率-------------------------------
########################################################################
mut_delSilent=read.csv("./result/mut_delSilent_info.csv")
length(unique(mut_delSilent$Tumor_Sample_Barcode)) #1815
mut_judge=paste0(mut_delSilent$Hugo_Symbol,":",mut_delSilent$Chromosome,":",mut_delSilent$Start_Position,":",mut_delSilent$HGVSp_Short)
a=sort(table(mut_judge),decreasing=T)
a=data.frame("HGVSp"=names(a),"sample_freq"=as.numeric(a))
mut_freq=a
write.csv(mut_freq,"./result/summutloc_freq1_delSilent.csv")

library(tidyr)
mut_freq=read.csv("./result/summutloc_freq1_delSilent.csv")
mut_freq2=mut_freq[order(mut_freq$sample_freq,decreasing=T),]
mut_freq2=separate(mut_freq2,HGVSp,c("Hugo_Symbol","Chromosome","Start_Position","HGVSp_Short"),":")
mut_freq2$HGVSp=paste0(mut_freq2$Hugo_Symbol,":",mut_freq2$Chromosome,":",mut_freq2$Start_Position,":",mut_freq2$HGVSp_Short)
mut_freq2=mut_freq2[which(mut_freq2$sample_freq > 1),]##挑选人群频率大于2的突变位点
total_sample=c()
total_sample_len=c()
for(i in 1:nrow(mut_freq2)){
  sample_idx=which(mut_judge == mut_freq2$HGVSp[i])
  
  total_sample=c(total_sample,unique(mut_delSilent$Tumor_Sample_Barcode[sample_idx]))
  total_sample_len=c(total_sample_len,length(unique(total_sample)))
}
mut_freq2$total_sample_len=total_sample_len
mut_freq2$accu_percent=round(total_sample_len/length(unique(mut_delSilent$Tumor_Sample_Barcode)),digits = 2)
write.csv(mut_freq2,"./result/summutloc_freq2_delSilent2.csv")



############################################################################################
#---------------------------------netmhcpan输入数据的构建-----------------------------------
############################################################################################
mut=read.csv("./result/mut_delSilent_info.csv")
mut_judge=paste0(mut$Hugo_Symbol,":",mut$Chromosome,":",mut$Start_Position,":",mut$HGVSp_Short)
mut_freq2=read.csv("./result/summutloc_freq2_delSilent2.csv")

mut_freq2$index <- 1:nrow(mut_freq2)
match_idx=match(mut_judge,mut_freq2$HGVSp)
mut$freq=mut_freq2$sample_freq[match_idx]
mut$index=mut_freq2$index[match_idx]
mut2=mut[-which(is.na(mut$freq)),-1]


mut <- na.omit(mut, select = index)
index <- min(which(mut_freq2$accu_percent >= 0.95)) #人群覆盖度达到95%
mut2 <- mut[mut$index <= index, ]

mut2=unique(mut2)
dim(mut2)
CHROM	= mut2$Chromosome
POS	= mut2$Start_Position
ID = rep(".",nrow(mut2))	
REF	= mut2$Reference_Allele
ALT	= mut2$Tumor_Seq_Allele2
QUAL = rep(100,nrow(mut2))
FILTER = rep("PASS",nrow(mut2))
INFO = paste0(mut2$Hugo_Symbol,":",mut2$HGVSp_Short)	
FORMAT = rep("GT:DP:VD:AF",nrow(mut2))
TCGA = rep("0/1:1000:100:0.1",nrow(mut2))
mut_vcf=data.frame(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,TCGA)
mut_vcf=unique(mut_vcf) 
table(nchar(mut_vcf$REF))
table(nchar(mut_vcf$ALT))
dim(mut_vcf)
length(which(nchar(mut_vcf$REF) == 1))
snv_idx=intersect(which(nchar(mut_vcf$REF) == 1),which(nchar(mut_vcf$ALT) == 1))
mut_vcf=mut_vcf[snv_idx,]
write.table(mut_vcf,"./result/summutloc_delSilent2.vcf",row.names = F,quote =F,sep = "\t")


########################################################################
#----------------/home/mall/neopep/Neopredpipe/-------------------------
#------------进行netmhcpan，预测peptide与MHC的结合亲性-----------------
#----------------------shell_netmhcpan.R---------------------------------
# ########################################################################
##shell
# source ~/.bashrc
# conda activate py27
# 查看配置文件usr_paths
# cd /home/mall/neopep/Neopredpipe/result
# 首先运行/home/mall/neopep/Neopredpipe/result/run.sh
# mkdir result
# 之后运行perl getResult.pl
# ####
# cd /home/mall/neopep/Neopredpipe/result/result
# for line in $(ls /home/mall/neopep/Neopredpipe/result/result)  #遍历文件夹下文件
# do
# sed -i '/^#/d' $line ###删除开头为#的行
# sed -i '/^-/d' $line  ###删除开头为-的行
# sed -i '/^$/d' $line  ###删除空行
# sed -i '/^[a-zA-Z]/d' $line ##删除以字母开头的所有行
# awk '!a[$0]++' $line  > $"{line}.tmp" && mv -f $"{line}.tmp" $line  ###删除重复行
# sed 's/^ *//' $line > $"{line}.tmp" && mv -f $"{line}.tmp" $line ##删除每行的空格
# done
########################################################################
#----------------------------所有预测多肽-------------------------------
########################################################################
sample=dir("./test")
library(tidyr)
sample_2=data.frame("x"=sample)
sample_2=separate(sample_2,x,c("x","sample","pep_len"),"_")
sample_3=unique(paste0(sample_2$x,"_",sample_2$sample))


for(i in 1:length(sample_3)){
  pep_8=read.delim(paste0("./test/",sample_3[i],"_8.txt"),header=F)
  pep_9=read.delim(paste0("./test/",sample_3[i],"_9.txt"),header=F)
  pep_10=read.delim(paste0("./test/",sample_3[i],"_10.txt"),header=F)
  pep_11=read.delim(paste0("./test/",sample_3[i],"_11.txt"),header=F)
  pep_12=read.delim(paste0("./test/",sample_3[i],"_12.txt"),header=F)
  pep_13=read.delim(paste0("./test/",sample_3[i],"_13.txt"),header=F)
  
  hla=paste0("HLA-A",unique(sample_2$sample)[i])
  
  
  pep_8$hla=rep(hla,nrow(pep_8))
  pep_9$hla=rep(hla,nrow(pep_9))
  pep_10$hla=rep(hla,nrow(pep_10))
  pep_11$hla=rep(hla,nrow(pep_11))
  pep_12$hla=rep(hla,nrow(pep_12))
  pep_13$hla=rep(hla,nrow(pep_13))
  
  
  if(i == 1){
    file=pep_8
    file=rbind(file,pep_9,pep_10,pep_11,pep_12,pep_13)
  }else{
    file=rbind(file,pep_8,pep_9,pep_10,pep_11,pep_12,pep_13)
  }
  
}
colnames(file)=c("Pos","Peptide","Core","Mut","Score","Aff(nM)","%Rank","BindLevel","HLA")
library(tidyr)
file_mut=data.frame("x"=file$Mut)
file_mut=separate(file_mut,x,c("gene","x"),":")
file$gene=file_mut$gene
file$loc=file_mut$x
write.csv(file,"./test/new_netmhcpan_result_delSilent.csv")


######################################################################################
######################################################################################
accu_sample_freq<-function(mut_delSilent,mut_freq){
  #mut_freq=mut_freq2
  total_sample=c()
  total_sample_len=c()
  #mut_judge=paste0(mut_delSilent$Hugo_Symbol,":",mut_delSilent$Chromosome,":",mut_delSilent$Start_Position,":",mut_delSilent$HGVSp_Short)
  mut_judge=paste0(mut_delSilent$Hugo_Symbol,":",mut_delSilent$HGVSp_Short)
  for(i in 1:nrow(mut_freq)){
    #sample_idx=which(mut_judge == mut_freq$HGVSp[i])
    sample_idx=which(mut_judge == mut_freq$Mut[i])
    total_sample=c(total_sample,unique(mut_delSilent$Tumor_Sample_Barcode[sample_idx]))
    total_sample_len=c(total_sample_len,length(unique(total_sample)))
  }
  mut_freq$total_sample_len=total_sample_len
  mut_freq$accu_percent=round(total_sample_len/length(unique(mut_delSilent$Tumor_Sample_Barcode)),digits = 2)
  
  return(mut_freq)
}

mut_delSilent=read.csv("./result/mut_delSilent_info.csv")
#mut_judge=paste0(mut_delSilent$Hugo_Symbol,":",mut_delSilent$Chromosome,":",mut_delSilent$Start_Position,":",mut_delSilent$HGVSp_Short)
mut_judge=paste0(mut_delSilent$Hugo_Symbol,":",mut_delSilent$HGVSp_Short)


file=read.csv("./test/new_netmhcpan_result_delSilent.csv")
file2=file[which(file$X.Rank <= 2),]
file_SB=file[which(file$X.Rank <= 0.5),]
file2_A0201=file2[which(file2$HLA == "HLA-A0201"),]
length(unique(file$gene) ) #9624
length(unique(file$Mut))   #21834
length(unique(paste0(file$Peptide,file$Mut))) #1374055
dim(file)                  #5441754      12
length(unique(file2$gene) ) #8229
length(unique(file2$Mut))   #16202
length(unique(file2$Peptide)) #75807
length(unique(paste0(file2$Peptide,file2$Mut)))
length(unique(paste0(file_SB$Peptide,file_SB$Mut)))
dim(file2)                  #95516    12
table(file$Mut)
table(table(file$Mut))     

mut_freq=read.csv("D:/CRC/result/big_library/summutloc_freq2_delSilent2.csv")  ##人群频率
mut_freq_judge=paste0(mut_freq$Hugo_Symbol,":",mut_freq$HGVSp_Short)
##提取高频突变位点
mut_freq_gene=mut_freq[1:max(which(mut_freq$accu_percent == 0.8)),]
mut_freq_gene_judge=paste0(mut_freq_gene$Hugo_Symbol,":",mut_freq_gene$HGVSp_Short)

##mut_freq2中仅保存netmhcpan出现的突变位点
mut_freq2=mut_freq[which(mut_freq_judge %in% file$Mut),]
mut_freq2_judge=paste0(mut_freq2$Hugo_Symbol,":",mut_freq2$HGVSp_Short)

mut_freq2$Mut=mut_freq2_judge
temp=mut_freq2[which(mut_freq2_judge %in% file$Mut),]
netmhcpan_samplefreq=accu_sample_freq(mut_delSilent,mut_freq2)  ##人群频率

pep_data_freq<-function(mut_delSilent,file2,mut_freq3){
  #mut_judge=paste0(mut_delSilent$Hugo_Symbol,":",mut_delSilent$HGVSp_Short)
  mut_freq3_judge=paste0(mut_freq3$Hugo_Symbol,":",mut_freq3$HGVSp_Short)
  total_sample=c()
  pep_data=c()
  for(i in 1:nrow(mut_freq3)){
    ##确定样本名
    # sample_idx=which(mut_judge == mut_freq3_judge[i])
    # sample_temp=unique(mut_delSilent$Tumor_Sample_Barcode[sample_idx])
    ##筛选netmhcpan的结果
    pep=file2[which(file2$Mut %in% mut_freq3_judge[i]),]
    pep=pep[order(pep$X.Rank),]
    
    
    if(nrow(pep)==1){
      pep_data=rbind(pep_data,pep[1,])
    }else if(nrow(pep)>1){
      pep_data=rbind(pep_data,pep[1:2,])
    }
    
    #total_sample=c(total_sample,sample_temp)
  }
  match_idx=match(pep_data$Mut,mut_freq3_judge)
  pep_data$HGVSp=mut_freq3$HGVSp[match_idx]
  return(pep_data)
}


mut_freq3=mut_freq2[which(mut_freq2_judge %in% file2$Mut),]
mut_freq3_judge=paste0(mut_freq3$Hugo_Symbol,":",mut_freq3$HGVSp_Short)
pep_data=pep_data_freq(mut_delSilent,file2,mut_freq3)
length(unique(mut_freq3$Hugo_Symbol) ) #8229
length(unique(mut_freq3_judge))        #16202
dim(mut_freq3)                         #16309


length(unique(pep_data$gene) ) #8229
length(unique(pep_data$Mut))   #16202
dim(pep_data)                  #29897   12

pep_datafreq=accu_sample_freq(mut_delSilent,pep_data)  
write.csv(pep_datafreq,"D:/CRC/result/big_library/pep_datafreq.csv")


min(which(pep_datafreq$accu_percent == 0.8))
pep_datafreq2=pep_datafreq[1:1111,]
A1101=pep_datafreq2[which(pep_datafreq2$HLA == "HLA-A1101"),]
A0201=pep_datafreq2[which(pep_datafreq2$HLA == "HLA-A0201"),]
A0301=pep_datafreq2[which(pep_datafreq2$HLA == "HLA-A0301"),]
A2402=pep_datafreq2[which(pep_datafreq2$HLA == "HLA-A2402"),]
A0201=accu_sample_freq(mut_delSilent,A0201) 
A0301=accu_sample_freq(mut_delSilent,A0301) 
A1101=accu_sample_freq(mut_delSilent,A1101) 
A2402=accu_sample_freq(mut_delSilent,A2402) 

A0201$accu_percent[nrow(A0201)]
A0301$accu_percent[nrow(A0301)]
A1101$accu_percent[nrow(A1101)]
A2402$accu_percent[nrow(A2402)]


##add high-frequency mut site
loc=mut_freq_gene_judge
HLA=c("HLA-A0201","HLA-A0301","HLA-A1101","HLA-A2402")
pep_data=c()
for(i in 1:length(loc)){
  file_top=file[file$Mut == loc[i],]
  for(j in 1:length(HLA)){
    pep=file_top[which(file_top$HLA==HLA[j]),]
    pep=pep[order(pep$X.Rank),]
    if(nrow(pep)==1){
      pep_data=rbind(pep_data,pep[1,])
    }else if(nrow(pep)>1){
      pep_data=rbind(pep_data,pep[1:2,])
    }
  }
}


a1=pep_data[which(pep_data$HLA == "HLA-A0201"),]
a2=pep_data[which(pep_data$HLA == "HLA-A0301"),]
a3=pep_data[which(pep_data$HLA == "HLA-A1101"),]
a4=pep_data[which(pep_data$HLA == "HLA-A2402"),]

A0201_2=unique(rbind(A0201[,2:12],pep_data[which(pep_data$HLA == "HLA-A0201"),-1]))
idx1=match(A0201_2$Mut,mut_freq2$Mut)
A0201_2$sample_freq=mut_freq2$sample_freq[idx1]
A0201_2=A0201_2[order(A0201_2$sample_freq,decreasing = T),]
A0201_2=accu_sample_freq(mut_delSilent,A0201_2)


A0301_2=unique(rbind(A0301[,2:12],pep_data[which(pep_data$HLA == "HLA-A0301"),-1]))
idx1=match(A0301_2$Mut,mut_freq2$Mut)
A0301_2$sample_freq=mut_freq2$sample_freq[idx1]
A0301_2=A0301_2[order(A0301_2$sample_freq,decreasing = T),]
A0301_2=accu_sample_freq(mut_delSilent,A0301_2)


A1101_2=unique(rbind(A1101[,2:12],pep_data[which(pep_data$HLA == "HLA-A1101"),-1]))
idx1=match(A1101_2$Mut,mut_freq2$Mut)
A1101_2$sample_freq=mut_freq2$sample_freq[idx1]
A1101_2=A1101_2[order(A1101_2$sample_freq,decreasing = T),]
A1101_2=accu_sample_freq(mut_delSilent,A1101_2)


A2402_2=unique(rbind(A2402[,2:12],pep_data[which(pep_data$HLA == "HLA-A2402"),-1]))
idx1=match(A2402_2$Mut,mut_freq2$Mut)
A2402_2$sample_freq=mut_freq2$sample_freq[idx1]
A2402_2=A2402_2[order(A2402_2$sample_freq,decreasing = T),]
A2402_2=accu_sample_freq(mut_delSilent,A2402_2)

A0201_2$accu_percent[nrow(A0201_2)]
A0301_2$accu_percent[nrow(A0301_2)]
A1101_2$accu_percent[nrow(A1101_2)]
A2402_2$accu_percent[nrow(A2402_2)]

data=rbind(A0201_2,A0301_2,A1101_2,A2402_2)
write.csv(data,"D:/CRC/result/big_library/pep_V1.csv")

data=read.csv("D:/CRC/result/big_library/pep_V1.csv")
length(unique(data$Peptide))
library(VennDiagram)
venn<-venn.diagram(x=list(A0201$Peptide,A0301$Peptide,A1101$Peptide,A2402$Peptide),
                   scaled = T, # 根据比例显示大小
                   alpha= 0.5, #透明度
                   lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF','#FFCCCC','#CCCCFF'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
                   label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
                   cex = 1, # 数字大小
                   fontface = "bold",  # 字体粗细；加粗bold
                   fill=c('#FFFFCC','#CCFFFF','#FFCCCC','#CCCCFF'), # 填充色 配色https://www.58pic.com/
                   category.names = c("HLA0201","HLA0301","HLA1101", "HLA2402") , #标签名
                   cat.dist = 0.02, # 标签距离圆圈的远近
                   cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
                   cat.cex = 1.5, #标签字体大小
                   cat.fontface = "bold",  # 标签字体加粗
                   cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
                   cat.default.pos = "outer",  # 标签位置, outer内;text 外
                   filename=NULL)
pdf('D:/CRC/result/big_library/fig/HLA_peptide.pdf')
grid.draw(venn)
dev.off()


venn<-venn.diagram(x=list(A0201_2$Peptide,A0301_2$Peptide,A1101_2$Peptide,A2402_2$Peptide),
                   scaled = T, # 根据比例显示大小
                   alpha= 0.5, #透明度
                   lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF','#FFCCCC','#CCCCFF'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
                   label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
                   cex = 1, # 数字大小
                   fontface = "bold",  # 字体粗细；加粗bold
                   fill=c('#FFFFCC','#CCFFFF','#FFCCCC','#CCCCFF'), # 填充色 配色https://www.58pic.com/
                   category.names = c("HLA0201","HLA0301","HLA1101", "HLA2402") , #标签名
                   cat.dist = 0.02, # 标签距离圆圈的远近
                   cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
                   cat.cex = 1.5, #标签字体大小
                   cat.fontface = "bold",  # 标签字体加粗
                   cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
                   cat.default.pos = "outer",  # 标签位置, outer内;text 外
                   filename=NULL)
pdf('D:/CRC/result/big_library/fig/new_HLA_peptide.pdf')
grid.draw(venn)
dev.off()


library(wordcloud2)
A0201_missingene=mut_freq_gene[which(mut_freq_gene_judge %in% setdiff(mut_freq_gene_judge,A0201$Mut)),] 
A0201_missingene$loc=paste0(A0201_missingene$Hugo_Symbol,":",A0201_missingene$HGVSp_Short)
A0201_missingene=A0201_missingene[which(A0201_missingene$loc %in% file$Mut),]
A0201_missingene2=data.frame("word"=A0201_missingene$loc,"freq"=A0201_missingene$sample_freq)
wordcloud2(A0201_missingene2,size = 0.5)

A0301_missingene=mut_freq_gene[which(mut_freq_gene_judge %in% setdiff(mut_freq_gene_judge,A0301$Mut)),] 
A0301_missingene$loc=paste0(A0301_missingene$Hugo_Symbol,":",A0301_missingene$HGVSp_Short)
A0301_missingene=A0301_missingene[which(A0301_missingene$loc %in% file$Mut),]
A0301_missingene2=data.frame("word"=A0301_missingene$loc,"freq"=A0301_missingene$sample_freq)
wordcloud2(A0301_missingene2,size = 0.5)


A1101_missingene=mut_freq_gene[which(mut_freq_gene_judge %in% setdiff(mut_freq_gene_judge,A1101$Mut)),] 
A1101_missingene$loc=paste0(A1101_missingene$Hugo_Symbol,":",A1101_missingene$HGVSp_Short)
A1101_missingene=A1101_missingene[which(A1101_missingene$loc %in% file$Mut),]
A1101_missingene2=data.frame("word"=A1101_missingene$loc,"freq"=A1101_missingene$sample_freq)
wordcloud2(A1101_missingene2,size = 0.5)



A2402_missingene=mut_freq_gene[which(mut_freq_gene_judge %in% setdiff(mut_freq_gene_judge,A2402$Mut)),] 
A2402_missingene$loc=paste0(A2402_missingene$Hugo_Symbol,":",A2402_missingene$HGVSp_Short)
A2402_missingene=A2402_missingene[which(A2402_missingene$loc %in% file$Mut),]
A2402_missingene2=data.frame("word"=A2402_missingene$loc,"freq"=A2402_missingene$sample_freq)
wordcloud2(A2402_missingene2,size = 0.5)




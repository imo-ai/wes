setwd("C:/project/neoantigen/pancreas_dna/")

library(tidyr)

###mutation
paad_cptac_2021<-read.delim("./raw_data/cbioportal/paad_cptac_2021/data_mutations.txt")
paad_icgc<-read.delim("./raw_data/cbioportal/paad_icgc/data_mutations.txt")
paad_qcmg_uq_2016<-read.delim("./raw_data/cbioportal/paad_qcmg_uq_2016/data_mutations_mskcc.txt")
paad_tcga<-read.delim("./raw_data/cbioportal/paad_tcga/data_mutations.txt")
#paad_tcga<-read.delim("./raw_data/cbioportal/paad_tcga/data_mutations2.xlsx.txt")
paad_tcga_pan_can_atlas_2018<-read.delim("./raw_data/cbioportal/paad_tcga_pan_can_atlas_2018/data_mutations.txt")
paad_utsw_2015<-read.delim("./raw_data/cbioportal/paad_utsw_2015/data_mutations_mskcc.txt")
pact_jhu<-read.delim("./raw_data/cbioportal/pact_jhu_2011/data_mutations_mskcc.txt")
paac_jhu<-read.delim("./raw_data/cbioportal/paac_jhu_2014/data_mutations_mskcc.txt")
panet_jhu<-read.delim("./raw_data/cbioportal/panet_jhu_2011/data_mutations_mskcc.txt")
panet<-read.delim("./raw_data/cbioportal/panet_shanghai_2013/data_mutations_mskcc.txt",skip=1)
panet_arcnet<-read.delim("./raw_data/cbioportal/panet_arcnet_2017/data_mutations_mskcc.txt",skip=1)


###
used_col=c("Hugo_Symbol","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","HGVSp_Short","Tumor_Sample_Barcode","Variant_Classification")
paad_cptac_2021=paad_cptac_2021[,which(colnames(paad_cptac_2021) %in% used_col)]
paad_icgc=paad_icgc[,which(colnames(paad_icgc) %in% used_col)]
paad_qcmg_uq_2016=paad_qcmg_uq_2016[,which(colnames(paad_qcmg_uq_2016) %in% used_col)]
paad_tcga=paad_tcga[,which(colnames(paad_tcga) %in% used_col)]
paad_tcga_pan_can_atlas_2018=paad_tcga_pan_can_atlas_2018[,which(colnames(paad_tcga_pan_can_atlas_2018) %in% used_col)]
paad_utsw_2015=paad_utsw_2015[,which(colnames(paad_utsw_2015) %in% used_col)]
pact_jhu=pact_jhu[,which(colnames(pact_jhu) %in% used_col)]
paac_jhu=paac_jhu[,which(colnames(paac_jhu) %in% used_col)]
panet_jhu=panet_jhu[,which(colnames(panet_jhu) %in% used_col)]
panet=panet[,which(colnames(panet) %in% used_col)]
panet_arcnet=panet_arcnet[,which(colnames(panet_arcnet) %in% used_col)]
mut=rbind(paad_cptac_2021,paad_icgc,paad_qcmg_uq_2016,paad_tcga,paad_tcga_pan_can_atlas_2018,
          paad_utsw_2015,pact_jhu,paac_jhu,panet_jhu,panet,panet_arcnet)
# length(unique(mut$Hugo_Symbol)) #paad癌肿共17538个基因，胰腺癌所有癌肿涉及17946个基因
# length(unique(paste0(mut$Hugo_Symbol,":",mut$HGVSp_Short))) #[1] paad癌肿共75451个位点，胰腺癌所有癌肿涉及81726个位点

# 数据清洗
mut=unique(mut)
#mut=mut[-which(mut$HGVSp_Short == ""),]
idx=grep("chr",mut$Chromosome)
mut$Chromosome[-idx]=paste0("chr",mut$Chromosome[-idx])
mut_delSilent=mut[-which(mut$Variant_Classification == "Silent"),]
mut_Silent=mut[which(mut$Variant_Classification == "Silent"),]
write.csv(mut,"./result/mut_info2.csv")
write.csv(mut_delSilent,"./result/mut_delSilent_info.csv")
write.csv(mut_Silent,"./result/mut_Silent_info.csv")




############################################################################################
#------------------------------计算突变位点的人群累计人群频率-------------------------------
############################################################################################

mut_delSilent <- read.csv("./result/mut_delSilent_info.csv")
length(unique(mut_delSilent$Tumor_Sample_Barcode))  ##1079
mut_judge=paste0(mut_delSilent$Hugo_Symbol,":",mut_delSilent$Chromosome,":",mut_delSilent$Start_Position,":",mut_delSilent$HGVSp_Short)


# 初始化一个空列表来存储每个突变位点的样本ID
sample_ids_list <- list()

# 计算突变在样本中发生频率
sample_freq=c()
for(i in 1:length(unique(mut_judge))){
  sample_idx=which(mut_judge == unique(mut_judge)[i])
  sample_freq=c(sample_freq,length(unique(mut_delSilent$Tumor_Sample_Barcode[sample_idx])))        ##根据突变id提取对应的sample_barcode，并去重，计算不同样本的数量
  unique_samples <- unique(mut_delSilent$Tumor_Sample_Barcode[sample_idx])
  sample_ids_list[[unique(mut_judge)[i]]] <- unique_samples # 存储每个突变位点的样本ID
}
mut_freq=data.frame("HGVSp"=unique(mut_judge),sample_freq)
mut_freq2=mut_freq[order(mut_freq$sample_freq,decreasing=T),]     # 排序



# 计算每个突变的累积样本数量
#total_sample=c()
#total_sample_len=c()
#for(i in 1:nrow(mut_freq2)){
#  sample_idx=which(mut_judge == mut_freq2$HGVSp[i])
#  total_sample=c(total_sample,unique(mut_delSilent$Tumor_Sample_Barcode[sample_idx]))
#  total_sample_len=c(total_sample_len,length(unique(total_sample)))
#}


total_sample <- c()  # 初始化累积样本向量
total_sample_len <- c()  # 初始化累积样本长度向量

# 遍历mut_freq2中的每个突变位点
for (i in 1:nrow(mut_freq2)) {
  sample_idx <- which(mut_judge == mut_freq2$HGVSp[i])  # 找到当前突变位点对应的样本索引
  
  # 提取当前突变位点的唯一样本ID
  current_samples <- unique(mut_delSilent$Tumor_Sample_Barcode[sample_idx])
  
  # 如果是第一个突变位点，初始化total_sample，否则将当前样本ID合并到total_sample中
  if (i == 1) {
    total_sample <- current_samples
  } else {
    # 使用union函数合并两个向量，并保留唯一值
    total_sample <- union(total_sample, current_samples)
  }

  # 计算累积的唯一样本数量
  total_sample_len <- c(total_sample_len, length(unique(total_sample)))
}

# 将累积样本长度添加到mut_freq2数据框中
mut_freq2$total_sample_len <- total_sample_len

# 打印结果
head(mut_freq2)



mut_freq2$accu_percent=round(total_sample_len/1079,digits = 2)
write.csv(mut_freq2,"./result/summutloc_freq_delSilent.csv")







############################################################################################
#---------------------------------netmhcpan输入数据的构建-----------------------------------
############################################################################################
mut=read.csv("./result/mut_delSilent_info.csv")
mut_judge=paste0(mut$Hugo_Symbol,":",mut$Chromosome,":",mut$Start_Position,":",mut$HGVSp_Short)



#------------->95%--------------#
mut_freq2 <- read.csv("./result/summutloc_freq_delSilent.csv")
dim(mut_freq2)
mut_freq2$index <- 1:nrow(mut_freq2)
match_idx=match(mut_judge,mut_freq2$HGVSp)
mut$freq=mut_freq2$sample_freq[match_idx]
mut$index=mut_freq2$index[match_idx]
head(mut)

index <- min(which(mut_freq2$accu_percent >= 0.8)) #人群覆盖度达到95%
mut2 <- mut[mut$index <= 59176, ]
#--------------------------------#



#-------------freq>2-------------#
match_idx=match(mut_judge,mut_freq2$HGVSp)
mut$freq=mut_freq2$sample_freq[match_idx]
#mut2=mut[which(mut$freq > 1),]
length(which(mut_freq2$sample_freq >2)) #201个突变位点  184
length(which(mut_freq2$sample_freq >3)) #93个突变位点   89
mut2=mut[which(mut$freq > 2),-1]
#--------------------------------#



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
mut_vcf=unique(mut_vcf)  #213 与201数值不一致的原因是，ALT中碱基不一致，但是HGVSp相同
dim(mut_vcf)
# table(mut_vcf$INFO)
# mut_vcf[which(mut_vcf$INFO == "TP53:p.X126_splice"),]
# mut_vcf[which(mut_vcf$INFO == "TP53:p.X261_splice"),]
write.table(mut_vcf,"./result/summutloc_delSilent.vcf",row.names = F,quote =F,sep = "\t")



########################################################################
#----------------/home/mall/neopep/Neopredpipe/-------------------------
#------------进行netmhcpan，预测peptide与MHC的结合亲性-----------------
#----------------------shell_netmhcpan.R---------------------------------
########################################################################
##shell
# source ~/.bashrc
# conda activate py27
# 查看配置文件usr_paths
# cd /home/mall/neopep/Neopredpipe/result
# 首先运行/home/mall/neopep/Neopredpipe/result/run.sh
## nohup bash run.sh &

###---------------run.sh----------------###
###########################################
##shell

#nohup python /home/mall/neopep/Neopredpipe/NeoPredPipe.py -I vcf/ -H vcf/HLA_0101.txt -o result_0101/ -n tcga -c 0 -d -l -E 8 9 10 11 12 13 &
#nohup python /home/mall/neopep/Neopredpipe/NeoPredPipe.py -I vcf/ -H vcf/HLA_0201.txt -o result_0201/ -n tcga -c 0 -d -l -E 8 9 10 11 12 13 &
#nohup python /home/mall/neopep/Neopredpipe/NeoPredPipe.py -I vcf/ -H vcf/HLA_0203.txt -o result_0203/ -n tcga -c 0 -d -l -E 8 9 10 11 12 13 &

## 必须输入
# -I 输入VCF文件目录位置
# -H VCF患者样本的HLA文件，或运行POLYSOLVER后包含患者特定目录的文件
# -o 输出结果的pathway
# -n 输出结果的名字

## 选择输入
#-l：指定是否删除ANNOVAR日志文件。默认为True。注意：用于调试。
#-d：指定是否删除程序创建的中间文件。默认为True。注意：设置此标志以恢复作业。
#-E EPITOPES [EPITOPES ...]， --epitopes EPITOPES [EPITOPES ...]：预测表位的长度。默认为8、9、10。
#-c COLREGIONS [COLREGIONS ...]：在多个区域VCF文件中，格式字段之后的列中不正常区域的列。例如：0是测试样本中的正常，肿瘤是其他列。程序可以处理每个VCF文件的不同数量的区域。

###########################################
###---------------run.sh----------------###





#mkdir result
# 之后运行perl getResult.pl

result_HLAxxxx_length.txt

#####
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
setwd("C:/project/neoantigen/pancreas_dna")
sample=dir("./test/")
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
write.csv(file,"./test/netmhcpan_result_delSilent.csv")


##筛选每个突变位点结合亲性最高的5条多肽
mut_freq2=read.csv("./result/summutloc_freq_delSilent.csv")
library(tidyr)
mut_freq3=separate(mut_freq2,HGVSp,c("gene","chr","start","Hgvsp"),":")
mut_freq3_judge=paste0(mut_freq3$gene,":",mut_freq3$Hgvsp)
mut_freq3_judge[1:10]
file=read.csv("./test/netmhcpan_result_delSilent.csv",row.names = 1)
file_Mut=unique(file$Mut)
result=c()
pop_num=c()
pop_num_idx=c()
for(i in 1:length(file_Mut)){
  mut_idx=which(file$Mut == file_Mut[i])
  temp=file[mut_idx,]
  temp=temp[order(temp$X.Rank,decreasing = F),]
  temp=temp[1:5,]
  result=rbind(result,temp)
  freq_idx=which(mut_freq3_judge == file_Mut[i])
  pop_num=c(pop_num,mut_freq3$sample_freq[min(freq_idx)])
  pop_num_idx=c(pop_num_idx,min(freq_idx))
}
result$pop_num=rep(pop_num,each=5)
result2=result[order(result$pop_num,decreasing = T),]


mut_delSilent=read.csv("./result/mut_delSilent_info.csv")
HGVSp=mut_freq2$HGVSp[sort(pop_num_idx)]
mut_judge=paste0(mut_delSilent$Hugo_Symbol,":",mut_delSilent$Chromosome,":",mut_delSilent$Start_Position,":",mut_delSilent$HGVSp_Short)
sample_freq=c()
total_sample=c()
total_sample_len=c()
for(i in 1:length(HGVSp)){
  sample_idx=which(mut_judge == HGVSp[i])
  sample_freq=c(sample_freq,length(unique(mut_delSilent$Tumor_Sample_Barcode[sample_idx])))
  
  total_sample=c(total_sample,unique(mut_delSilent$Tumor_Sample_Barcode[sample_idx]))
  total_sample_len=c(total_sample_len,length(unique(total_sample)))
}
freq=data.frame("HGVSp"=HGVSp,sample_freq,total_sample_len)
freq$accu_percent=round(total_sample_len/1082,digits = 2)
freq2=freq[order(freq$sample_freq,decreasing=T),]
result2$sample_freq=rep(freq2$sample_freq,each=5)
result2$total_sample_len=rep(freq2$total_sample_len,each=5)
result2$accu_percent=rep(freq2$accu_percent,each=5)
write.csv(result2[,-12],"./test/peptide.csv")

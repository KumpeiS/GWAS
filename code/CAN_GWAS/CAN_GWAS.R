# 設定
n <- 60 #何回繰り返し共変数を加えGWASするのか

# 必要パッケージのロード
library(vcfR)
library(rrBLUP)
library(qqman)
library(tidyverse)


# vcf形式の遺伝子型データの読み込み
vcf <- read.vcfR("360.vcf")

# 遺伝子型データの抽出
gt <- extract.gt(vcf)
# マーカー情報：染色体番号・位置データの抽出
chrom <- getCHROM(vcf)
pos <- getPOS(vcf)


# 遺伝子型を文字型（塩基）から数字に変換
# 変換後の遺伝子型の表の入れ物作成
gt.score <- matrix(NA, nrow(gt), ncol(gt))
# "1", "-1"はそれぞれ1 −1へ変換
gt.score[gt == "1"] <- 1
gt.score[gt == "-1"] <- -1
# 行列名を表に追加
rownames(gt.score) <- rownames(gt)
colnames(gt.score) <- colnames(gt)
# 行列を転座
gt.score <- t(gt.score)

# 表現型データの読み込み
pheno <- read.csv("360_pheno.csv", row.names = 1, check.names=FALSE)

# 共変数なしのGWAS
Pheno<-na.omit(pheno) # 形質ごとに欠測を除く
line <- intersect(rownames(Pheno), colnames(gt))
y <- Pheno[line,]   # データを読み込み

# 遺伝子型の読み込み
Geno <- gt.score[line,]
x <- as.matrix(Geno)

# 血縁関係も含めた解析
amat <- A.mat(x, shrink = T)        # 血縁行列の計算

# rrBLUPの準備
g <- data.frame(colnames(x), as.numeric(chrom), pos, t(x))
rownames(g) <- 1:nrow(g)
colnames(g) <- c("marker", "chrom", "pos", rownames(x))
p <- data.frame(rownames(x), y)
colnames(p)
colnames(p) <- c("gid","y")
colnames(amat) <- rownames(amat) <- rownames(x)

# rrBLUPパッケージのGWAS関数を使用してGAWAS
gwas <- GWAS(p, g, n.PC = 4, K = amat, plot = F, min.MAF=0.02)

# 結果を入れておく準備
result<- matrix(NA, nrow = ncol(Geno), ncol = 1)
rownames(result)<-colnames(Geno)
colnames(result)<-colnames(Pheno)[2]

# 結果の保存
result[,1]<-gwas[,4]

#make folder
dir.create("./")

# 結果をcsvに保存
write.csv(result,paste0("./0/QK-0qtn.csv"),quote=F)

# 結果をマンハッタンプロットとして散布
mht <- data.frame(SNP = gwas$marker, CHR = gwas$chrom, BP = gwas$pos, P = 10^(-gwas$y))
mht <- na.omit(mht)
pdf(paste0("../0/QK-0qtn.pdf"), width = 10, height = 5)
manhattan(mht)
dev.off()

# qqplot
pdf(paste0("./0/QK-0qtn_qqplot.pdf"), width = 5, height = 5)
qq(mht$P)
dev.off()

# データの読み込み
data_id <- read.csv("all_flowering_SNPs.csv")
colnames(data_id)

data_gwas <- read.csv(paste0("./0/QK-0qtn.csv"))
colnames(data_gwas) <- c("ID", "LOD")
colnames(data_gwas)

# IDで情報を結合する
left_joined_data <- left_join(data_id, data_gwas, by = "ID")

write.csv(left_joined_data,paste0("./0/flowering_genes_0cofactors.csv"))


############################################################
############################################################
# 共変数ありのGWAS
for(k in 1:n){
  #k<-32
  
  #highest SNP in known flowering genes
  SNP_high <- read.csv(paste0("../",k-1,"/flowering_genes_",k-1,"cofactors.csv"))
  SNP_high_ID <- SNP_high %>% filter(LOD == max(LOD)) %>% select(ID)
  # make pheno data
  pheno_add <- read.csv("Geno_72HD.csv")
  pheno_add2 <- pheno_add %>% filter(ID==c(SNP_high_ID$ID[1]))
  pheno_add3 <- t(pheno_add2)
  pheno_add3 <- pheno_add3[-1,]
  pheno_add3 <- data.frame(pheno_add3)
  colnames(pheno_add3) <- c(SNP_high_ID$ID[1])
  # joint pheno and new fix
  pheno <- cbind(pheno,pheno_add3)
  
  Pheno<-na.omit(pheno) # 形質ごとに欠測を除く
  line <- intersect(rownames(Pheno), colnames(gt))
  y <- Pheno[line,]   # データを読み込み
  
  # 遺伝子型の読み込み
  Geno <- gt.score[line,]
  x <- as.matrix(Geno)
  
  # 血縁関係も含めた解析
  amat <- A.mat(x, shrink = T)        # 血縁行列の計算
  
  # rrBLUPの準備
  g <- data.frame(colnames(x), as.numeric(chrom), pos, t(x))
  rownames(g) <- 1:nrow(g)
  colnames(g) <- c("marker", "chrom", "pos", rownames(x))
  p <- data.frame(rownames(x), y)
  colnames(p)
  colnames(p) <- c("gid","y",as.character(seq(1,k)))
  colnames(amat) <- rownames(amat) <- rownames(x)
  
  # rrBLUPパッケージのGWAS関数を使用してGAWAS
  Fix <- as.character(seq(1,k))
  gwas <- GWAS(p, g, n.PC = 4, K = amat, plot = F, fixed=Fix, min.MAF=0.01)
  head(gwas)
  
  # 結果を入れておく準備
  result<- matrix(NA, nrow = ncol(Geno), ncol = 1)
  rownames(result)<-colnames(Geno)
  colnames(result)<-colnames(Pheno)[2]
  
  # 結果の保存
  result[,1]<-gwas[,4]
  
  #make folder
  dir.create(paste0("./",k))
  
  # 結果をcsvに保存
  write.csv(result,paste0("./",k,"/QK-",k,"qtn.csv"),quote=F)
  
  # 結果をマンハッタンプロットとして散布
  mht <- data.frame(SNP = gwas$marker, CHR = gwas$chrom, BP = gwas$pos, P = 10^(-gwas$y))
  mht <- na.omit(mht)
  pdf(paste0("./",k,"/QK-",k,"qtn.pdf"), width = 10, height = 5)
  manhattan(mht)
  dev.off()
  
  # qqplot
  pdf(paste0("./",k,"/",k,"qqplot.pdf"), width = 5, height = 5)
  qq(mht$P)
  dev.off()
  
  ############################################################
  ############################################################
  # データの読み込み
  data_id <- read.csv("all_flowering_SNPs.csv")
  colnames(data_id)
  
  data_gwas <- read.csv(paste0("./",k,"/QK-",k,"qtn.csv"))
  colnames(data_gwas) <- c("ID", "LOD")
  colnames(data_gwas)
  
  # IDで情報を結合する
  left_joined_data <- left_join(data_id, data_gwas, by = "ID")
  
  write.csv(left_joined_data,paste0("./",k,"/flowering_genes_",k,"cofactors.csv"))

}


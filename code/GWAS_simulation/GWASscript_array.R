# gwas_array_jpb

# 必要パッケージのロード
require(vcfR)
require(rrBLUP)

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
pheno <- read.csv("Pheno.csv", row.names = 1, check.names=FALSE)

# 遺伝子型の読み込み
Geno <- gt.score
x <- as.matrix(Geno)

# 血縁関係も含めた解析
amat <- A.mat(x, shrink = T)        # 血縁行列の計算

# phenotypeの読み込み
args <- commandArgs(trailingOnly = TRUE)   # 引数を取得
k <- as.numeric(args[1])                 # kに引1を設定
Pheno<-pheno[,k,drop=F]  #  表現型のうちk列目を取り出す
Pheno<-na.omit(Pheno) # 形質ごとに欠測を除く
line <- intersect(rownames(Pheno), colnames(gt))
y <- Pheno[line,]   # データを読み込み

# rrBLUPの準備
g <- data.frame(colnames(x), as.numeric(chrom), pos, t(x))
rownames(g) <- 1:nrow(g)
colnames(g) <- c("marker", "chrom", "pos", rownames(x))
p <- data.frame(rownames(x), y)
colnames(p) <- c("gid", "y")
colnames(amat) <- rownames(amat) <- rownames(x)

# rrBLUPパッケージのGWAS関数を使用してGAWAS
gwas <- GWAS(p, g, n.PC = 4, K = amat, plot = F, min.MAF = 0.02)
head(gwas)

# 結果を入れておく準備
result<- matrix(NA, nrow = ncol(gt.score), ncol = 1)
rownames(result)<-colnames(gt.score)
colnames(result)<-colnames(Pheno)

# 結果の保存
result[,1]<-gwas[,4]

# 結果をcsvに保存
write.csv(result,paste0(k,"_qk.csv"),quote=F)

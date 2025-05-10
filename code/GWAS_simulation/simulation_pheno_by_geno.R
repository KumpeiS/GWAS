# 遺伝子型から仮想形質データ算出スクリプト
# 設定
QTNs <- 20     # how many qtns are set as causal genes
kmax <- 1000     # how many patterns of QTNs
#############################################

for (k in 1:kmax){
  geno <- read.csv(paste0("./QTN/", k, "_snp", QTNs,"info.csv"), row.names = 1, header = T)
  dim(geno)

  # k=1のときに結果の入れ物を用意する
  if(k == 1){
    result <- matrix(NA, nrow = ncol(geno), ncol = kmax)
    rownames(result) <- colnames(geno)
  }
  
  x <- as.matrix(geno)
  # NAを0に置き換える
  x <- ifelse(is.na(x), 0, x)
  # 遺伝子型ごとに表現型データ作成
  x[x == "0"] <- 1
  x[x == "-1"] <- 0
  x[x == "1"] <- 2
  # 出穂期は50日を基準とする
  day <- matrix(50, 1, ncol(x)) 
  # QTN効果の加算
  for (i in 1:QTNs){
    day <- day + x[i,] + runif(1*ncol(x), min=-1, max=1)
  }
  
  y1 <- day 
  y1 <- t(y1)
  result[,k] <- y1
  #colnames(result[,k]) <- "test"#paste0(k,"_noise1")
}

# 結果の出力
write.csv(result, "Pheno.csv", quote=F)
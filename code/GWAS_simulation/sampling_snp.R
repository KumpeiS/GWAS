##################################
#同時に抽出するSNP数
n.qtn <- 20
#データセット数
n.dataset <- 1000
##################################

# SNPデータの読み込み
# QTN数ループ
#for (k in 1:n.qtn){
  data <- read.csv("360linesMAF0.05NA0.1_pruning.csv"), sep = ",", header = T, row.names = 1, check.names=FALSE)
#}


# データセット数ループ
for (i in 1:n.dataset) {
#i=1
List <- NULL
  #サンプル抽出
  # QTN数ループ
  for (j in 1:n.qtn){
    
    #data.select <- get(paste0("snp_data_", j))
    rows.sampled <- assign(paste0(j,".rows.sampled"), sample(nrow(data), 1))
    snp.sampled <- data[rows.sampled, , drop=F]
    List <- rbind(List, snp.sampled)
  }
  
  #出力
	write.csv(List, paste0("./QTN/",i,"_snp", n.qtn, "info.csv", sep=""))
}
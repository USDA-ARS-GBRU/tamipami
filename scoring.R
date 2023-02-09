library(rjson)

result <- fromJSON(file = "~/Documents/HT-TAMDA/output.txt")
df <- as.data.frame(result)

df$diff = df$cont_clr -df$exp_clr
df$zscore <- (df$diff - mean(df$diff)) / sd(df$diff)
df$pval <- pnorm(df$zscore, lower.tail=FALSE)
df$p.adjust <- p.adjust(df$pval, method="BH")

df[df$p.adjust < 0.05,]

df <- df[order(df$p.adjust),]
write.csv(df, "4mer_tam.csv")
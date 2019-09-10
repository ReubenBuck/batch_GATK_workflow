args = commandArgs(trailingOnly=TRUE)

ref= args[1]

amb <- read.table(paste(ref,"amb", sep = "."))
fai <- read.table(paste(ref,"fai", sep = "."))

rownames(fai) <- fai$V1
fai$cumsum <- c(0,cumsum(as.numeric(fai[1:(nrow(fai)-1),2])))
gapChr <- cut(amb$V1,breaks = c(fai$cumsum,amb[1,1]), labels = fai$V1, right = TRUE)
gapStart <- (amb$V1 - fai[as.character(gapChr),"cumsum"]) + 1
gapEnd <- gapStart + amb$V2

gaps <- data.frame(chr = as.character(gapChr), start = gapStart - 1, end = gapEnd - 1)[2:nrow(amb),]

write.table(x = gaps, file = stdout(), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)




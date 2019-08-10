library (SNPediaR)

res <- getCategoryElements(category = "In_dbSNP")

# Select 5000 SNPs randomly
randomIndexes <- sample.int(107085, 5000) 
selection <- res[randomIndexes]

# Extract chromosome and position info
dataList <- getPages(titles = selection, wikiParseFunction = extractSnpTags, tags = c("Chromosome", "position"))

# Put information into a data frame 
snpInfo <- data.frame(matrix(unlist(dataList), nrow=length(dataList), byrow=T))
colnames(snpInfo) <- c("chromosome","position")
rownames(snpInfo) <- selection

write.csv(snpInfo, file="snpList copy.csv")
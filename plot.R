library(tidyverse)
library(cowplot)
anot <- read.csv2("~/projects/Marine_Lab/data/ref/annotation.bed", header=FALSE, stringsAsFactors=FALSE, sep="\t")
anot$V3 <- as.character(anot$V3)
anot <- anot[, c("V1","V2", "V3", "V7", "V9")]

data <- read.csv2("TresCombo_toR.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
data <- data %>%
	mutate(POS = paste(X.1.CHROM, X.2.POS, sep="_")) %>%
	select(POS, ends_with("AD"))

data <- split(data, data$POS)
data <- lapply(data, function(x) as.matrix(x[-1]))

dataNames <- lapply(data, function(x) apply(x, 1, function(y) names(subset(y, y>0))))
dataNames <- lapply(dataNames, function(x) if(class(x)=="character"){return(as.list(x))}else{return(x)}) 
dataNames <- lapply(dataNames, function(x) list(REF=unname(unlist(x[1])), ALT=unique(Reduce(c,x[-1]))))
dataNames <- lapply(dataNames, function(x) {x$REF <- x$REF[!x$REF %in% x$ALT]; x})
dataNames <- lapply(dataNames, function(x) lapply(x, length))

procData <- as.data.frame(unlist(dataNames))
procData$ID <- rownames(procData)
rownames(procData) <- c()
colnames(procData)  <-  c("CELLS", "ID")
P1 <- procData %>%
	separate(ID, c("CHROM", "LOC"), "_" ) %>%
	separate(LOC, c("LOC", "ALLELE")) %>%
	group_by(CHROM, LOC) %>%
	mutate(ALL=sum(CELLS)) %>%
	mutate(PERC=CELLS/ALL*100) %>%
	left_join(., anot, by=c("CHROM"= "V1", "LOC" = "V3")) %>%
	select(-V2) %>%
	mutate(ID = paste(V7, LOC, sep="_")) %>% 
	filter(ALL>10) %>%
	group_by(CHROM, LOC) %>%
	mutate(diff = PERC[1]-PERC[2]) %>%
	filter(diff<60)

P1 %>%ggplot(aes(ID, CELLS, fill=ALLELE)) + 
geom_bar(stat="identity") + 
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
scale_y_continuous(expand = c(0,0)) 


#write.csv(procData, file = "procData_Trescombo")
#write.csv(P1$ID, file = "ID_List_Trescombo")


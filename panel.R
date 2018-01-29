library(tidyverse)
library(cowplot)


data_T0 <- read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/combo/procData_T0combo", colClasses=c("NULL", NA, NA))
data_T4 <- read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/combo/procData_T4combo", colClasses=c("NULL", NA, NA))
data_T28 <- read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/combo/procData_T28combo", colClasses=c("NULL", NA, NA))
data_Tres <- read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/combo/procData_Trescombo", colClasses=c("NULL", NA, NA))

data_T0$ORIGIN <- "T0"
data_T4$ORIGIN <- "T4"
data_T28$ORIGIN <- "T28"
data_Tres$ORIGIN <- "Tres"

data <- rbind(data_T0,
	      data_T4,
	      data_T28,
	      data_Tres)

anot <- read.csv2("~/projects/Marine_Lab/data/ref/annotation.bed", header=FALSE, stringsAsFactors = FALSE, sep = "\t")
anot$V3 <- as.character(anot$V3)
anot <- anot[, c("V1", "V2", "V3", "V7", "V9")]


#ID_T0 <- unique(as.character(read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/T0/ID_List_T0")$x))
#ID_T4 <- unique(as.character(read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/T4/ID_List_T4")$x))
#ID_T28 <- unique(as.character(read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/T28/ID_List_T28")$x))
#ID_Tres <- unique(as.character(read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/Tres/ID_List_Tres")$x))
#ID_List <- unique(c(ID_T0, ID_T4, ID_T28, ID_Tres))

ID_T0 <- unique(as.character(read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/combo/ID_List_T0combo")$x))
ID_T4 <- unique(as.character(read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/combo/ID_List_T4combo")$x))
ID_T28 <- unique(as.character(read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/combo/ID_List_T28combo")$x))
ID_Tres <- unique(as.character(read.csv("~/projects/Marine_Lab/data/smartseq2/noMerge/combo/ID_List_Trescombo")$x))
ID_List <- unique(c(ID_T0, ID_T4, ID_T28, ID_Tres))

pdf("panel_combo60_comboselCELLS.pdf")
data %>%
	separate(ID, c("CHROM", "LOC"), "_") %>%
	separate(LOC, c("LOC", "ALLELE")) %>%
	mutate(ORIGIN = factor(ORIGIN, levels = c("T0", "T4", "T28", "Tres"))) %>%
	group_by(ORIGIN, CHROM, LOC) %>%
	mutate(ALL=sum(CELLS)) %>%
	mutate(PERC=CELLS/ALL*100) %>%
	left_join(., anot, by=c("CHROM" = "V1", "LOC" = "V3")) %>%
	select(-V2) %>%
	mutate(ID = paste(V7, LOC, sep="_")) %>%
	filter(ID %in% ID_List) %>%
	ggplot(aes(ID, CELLS, fill = ALLELE)) +
	geom_bar(stat = "identity") +
	theme(axis.text.x = element_text(angle=90, hjust=1, size = 8)) +
	scale_y_continuous(expand=c(0,0))+
	facet_wrap(~ORIGIN, ncol=1)
dev.off()

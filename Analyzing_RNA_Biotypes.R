library(dplyr)
library(ggplot)

MI_Cancer <- read.table("/Users/kenzuke/Documents/Paper_1/MM_MI_100k.sif", header = TRUE)
Annot <- read.csv("/Users/kenzuke/Documents/R_projects/Intersections_2/mart_export.txt")
Annot$name <- gsub("-", ".", Annot$name)
Annot$name <- gsub("_", ".", Annot$name)
table(as.factor(Annot$Gene.type))
Annot$Gene.type <- gsub(".*pseudogene", "pseudogene", Annot$Gene.type)
Annot$Gene.type <- gsub("IG.*", "IG", Annot$Gene.type)
Annot$Gene.type <- gsub("TR.*", "TR", Annot$Gene.type)
table(as.factor(Annot$Gene.type))


# Hay que revisar Annot para editar con Annot$Gene.type <- gsub("TR.*", "TR", Annot$Gene.type) #

Sources <- merge(MI_Cancer[1:100000,], Annot, by.x = "source", by.y = "name")
Sources <- Sources %>% dplyr::select("source", "mi", "target", "Gene.type")
colnames(Sources) <- c("source", "mi", "target", "Gene.type.source")
Targets <- merge(Sources, Annot, by.x = "target", by.y = "name")
Targets <- Targets %>% dplyr::select("source", "mi", "target", "Gene.type.source", "Gene.type")
colnames(Targets) <- c("source", "mi", "target", "Gene.type.source", "Gene.type.target")

for(i in 1:nrow(Targets)) {
  
  if(Targets$Gene.type.source[i] == Targets$Gene.type.target[i]) {
    Targets$Homo_Hetero[i] <- "Homo"
    Targets$Genes.linked[i] <- paste(Targets$Gene.type.source[i])} else {
      Targets$Homo_Hetero[i] <- "Hetero"
      Targets$Genes.linked[i] <- paste(pmin(Targets$Gene.type.source[i],Targets$Gene.type.target[i]), 
                                       pmax(Targets$Gene.type.source[i],Targets$Gene.type.target[i]),sep="-")
      
      Targets[order(abs(Targets$mi), decreasing = TRUE),]
    }
  
}


AML_100k_EdgesTypes <- as.data.frame((table(as.factor(Targets$Genes.linked))))
MM_100k_EdgesTypes <- as.data.frame((table(as.factor(Targets$Genes.linked))))
TALL_100k_EdgesTypes <- as.data.frame((table(as.factor(Targets$Genes.linked))))
BALL_100k_EdgesTypes <- as.data.frame((table(as.factor(Targets$Genes.linked))))

AML_100k_EdgesTypes <- AML_100k_EdgesTypes %>% mutate(fraction = Freq/100000, cancer = )
MM_100k_EdgesTypes <- MM_100k_EdgesTypes %>% mutate(fraction = Freq/100000)
TALL_100k_EdgesTypes <- TALL_100k_EdgesTypes %>% mutate(fraction = Freq/100000)
BALL_100k_EdgesTypes <- BALL_100k_EdgesTypes %>% mutate(fraction = Freq/100000)

for(i in 1:nrow(Targets)) {
  Hom_links <- data.frame(cancer="AML", link_type=paste(Targets$Genes.linked[i]), 
                          fraction=(nrow(Targets[Targets$Genes.linked == paste(Targets$Genes.linked[i]),])/nrow(Targets)))
  
  if(i==1){
    x <- Hom_links
  } else {
    x <- rbind(x, Hom_links)
  }
}
x <- x[!duplicated(x$link_type),]

x <- x[order(x$fraction, decreasing = TRUE),]
x <- x[1:5,]
x

Hom_fraction <- data.frame(cancer="NormalAML", link_type="Homo", fraction=(nrow(Targets[Targets$Homo_Hetero == "Homo",])/nrow(Targets)))
Het_fraction <- data.frame(cancer="NormalAML", link_type="Hetero", fraction=(nrow(Targets[Targets$Homo_Hetero == "Hetero",])/nrow(Targets)))
HH_fraction <- rbind(Hom_fraction, Het_fraction)

hh_fraction_plot <- ggplot(HH_fraction, aes(fill=link_type, x=cancer, y=fraction)) + geom_bar(position="stack", stat="identity")

link_types <- ggplot(binding, aes(fill=link_type, x=cancer, y=fraction)) + geom_bar(position="stack", stat="identity") + ylim(0,1)
link_types + facet_grid(rows = vars(V4))  + scale_fill_brewer(palette = "Paired") +
  theme(plot.title = element_text(face = "bold")) +
  ggtitle("Fractions of Top 5 of interactions \nbetween RNA biotype in the networks") +
  labs(y = "Fraction of initeractions by biotype", x = "Cancer and respective control")

binding$V4_f = factor(binding$V4, levels=c("Normal", "Cancer"))

AML_top5Links <- x
MM_top5Links <- x
BALL_top5Links <- x
TALL_top5Links <- x

NormalAML_top5Links <- x
NormalBALL_top5Links <- x
NormalTALL_top5Links <- x
NormalMM_top5Links <- x

binding <- rbind(AML_top5Links, NormalAML_top5Links, NormalBALL_top5Links, NormalTALL_top5Links, NormalMM_top5Links, MM_top5Links,BALL_top5Links,TALL_top5Links)

binding[6:25, 4] <- c(rep("Normal", 20))
binding[26:40, 4] <- c(rep("Cancer", 15))
binding[1:5, 4] <- c(rep("Cancer", 5))

binding[6:10, 1] <- c(rep("AML", 5))
binding[11:15, 1] <- c(rep("BALL", 5))
binding[16:20, 1] <- c(rep("TALL", 5))
binding[21:25, 1] <- c(rep("MM", 5))

binding_new <- binding
binding_new$V4_f <- factor(binding_new$V4, levels = c("Normal", "Cancer"))
binding_new$cancer <- factor(binding_new$cancer, levels = c("AML", "BALL", "TALL", "MM"))

link_types <- ggplot(binding_new, aes(fill=link_type, x=cancer, y=fraction)) + geom_bar(position="stack", stat="identity") + ylim(0,1)
BM_BiotypeLinked <- link_types + facet_grid(rows = vars(V4_f))  +  scale_fill_manual(values = cols_barplot) +theme(plot.title = element_text(face = "bold")) +
  ggtitle("Fractions of Top 5 of interactions \nbetween RNA biotype in the MI networks") +
  labs(y = "Fraction of initeractions by biotype", x = "Cancer and respective control")
BM_BiotypeLinked + scale_fill_manual(values = cols_barplot)



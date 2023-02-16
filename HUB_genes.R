# Set the working directory (this is where your source files are and where the output figures will be saved).
setwd(choose.dir(getwd()))

# Install these packages, if you don't have them installed.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("STRINGdb")
install.packages("tidyverse")
install.packages("igraph")

# Load libraries.
library(STRINGdb)
library(tidyverse)
library(igraph)

# Check the newest string version available and the code of the organism that you are working on: https://string-db.org
string_db <- STRINGdb$new( version="11", species=3702, 
                           score_threshold=10000, input_directory="")

# Download protein network data (full network, scored links between proteins) for your organism from "Download" tab on: https://string-db.org
# Here is the network for Arabidopsis thaliana
network <- read.delim("3702.protein.links.v11.5.txt", sep=" ")

# Load sample gene set
set1 <- read.csv("sample_set.csv", sep = ",")

# Map the genes and create network
mapped <- string_db$map( set1, "target", removeUnmappedRows = TRUE )
m <- merge(mapped, network, by.x = "STRING_id", by.y = "protein1")
m2 <- merge(mapped, m, by.x = "STRING_id", by.y = "protein2")
edge_list <- tibble(from = m2$target.x, to = m2$target.y)
routes_igraph <- graph_from_data_frame(d = edge_list, directed = TRUE)

# Now we will find genes with top 10 "degree" score. You can calculate hubs with other methods like "eccentricity".
deg <- as.data.frame(degree(routes_igraph))
deg2 <- deg[ order((deg$`degree(routes_igraph)`), row.names(deg),decreasing=T), ,drop=F]
deg2 <- cbind(rownames(deg2), data.frame(deg2, row.names=NULL))
deg3 <- as.data.frame(deg2[1:10,])
colnames(deg3) = c("target", "value")

# Map the genes and create new network for the identified hub genes.
mapped <- string_db$map( deg3, "target", removeUnmappedRows = TRUE )
m <- merge(mapped, network, by.x = "STRING_id", by.y = "protein1")
m2 <- merge(mapped, m, by.x = "STRING_id", by.y = "protein2")
edge_list <- tibble(from = m2$target.x, to = m2$target.y)
routes_igraph <- graph_from_data_frame(d = edge_list, directed = TRUE)

# You can save identified hub gene list.
name <- paste("degree.csv")
write.csv(deg3, name, quote = F)


# Plot the network. You can mark group of functionally related genes with different colors.
# Numbers used for the marking are the positions of selected genes on the list that you can check:
# V(routes_igraph)$name
plot(routes_igraph ,mark.groups=list(c(1,4,9),c(8,2,3,5,6,7)),
     mark.col=c("#F2F9EC","#E5F9FF"), mark.border=NA, vertex.color="grey", vertex.label.color="black",
     vertex.frame.color="grey",vertex.label.family="Arial", vertex.label.cex=1,edge.arrow.size 	= .2, 
     vertex.shape="none",vertex.label.font=2)

# Save the plot.
png("HUB genes degree test.png", 500, 500)
plot(routes_igraph ,mark.groups=list(c(1,4,9),c(8,2,3,5,6,7)),
     mark.col=c("#F2F9EC","#E5F9FF"), mark.border=NA, vertex.color="grey", vertex.label.color="black",
     vertex.frame.color="grey",vertex.label.family="Arial", vertex.label.cex=1,edge.arrow.size 	= .2, 
     vertex.shape="none",vertex.label.font=2)
dev.off()

# End of the code.
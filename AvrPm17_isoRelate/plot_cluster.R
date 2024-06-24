setwd("~/projects/project_avrpm17/isoRelate/")
library(isoRelate,lib = "~/data/R_lib")
#load("BgtE+r_avrpm17_geno.RData")
load("BgtE+_avrpm17_geno.RData")

library(ggplot2)
library(maps)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(patchwork)

sessionInfo()

## make clusters with isoRelate
my_i_clusters <- getIBDiclusters(ped.genotypes = my_genotypes, 
                                 ibd.segments = my_ibd, 
                                 interval = c("1", 4365017, 4370466), # From first nucleotide of CDS of first copy, to last nucleotide of CDS of second copy
                                 prop=1, 
                                 hi.clust = FALSE)

g<-my_i_clusters$i.network
Samples <- my_genotypes$pedigree$fid
samples <- paste(Samples, Samples, sep = "/")



my_groups <- my_genotypes[[1]][,1:3]
my_groups


# read variants (mature protein aa)
variants<-read.csv("~/projects/project_avrpm17/haplotypes/hap_class_mp_trns_no_miss.csv")
colnames(variants) <- c("variants","Sample.Name")
variants$Sample.Name <- gsub("_1_avrpm17_cds", "",variants$Sample.Name)
variants$Sample.Name <- gsub("_2_avrpm17_cds", "",variants$Sample.Name)

# read dna variants
dna_variants<-read.csv("~/projects/project_avrpm17/haplotypes/hap_class_cds_RC_no_miss.csv")
colnames(dna_variants) <- c("dna_variants","Sample.Name")
dna_variants$Sample.Name <- gsub("_1_avrpm17_cds", "",dna_variants$Sample.Name)
dna_variants$Sample.Name <- gsub("_2_avrpm17_cds", "",dna_variants$Sample.Name)


# read phenotype
phenotype<-read.csv("~/projects/project_avrpm17/Bgt_Amigo_pheno_2022_2023_Poland.csv")


# read meta, extract samples and reorder them based on graph vertex order

meta<-read.csv("~/projects/vcf_project_tritici/2022+before2022+2023+ncsu_metadata+fs+admxK9_03052024.csv")
meta1<-data.frame(meta$Sample.Name,meta$Country,meta$fs_level_4,meta$Longitude,meta$Latitude,meta$Year.of.Collection)
colnames(meta1) <- c("Sample.Name","Country", "Population","Longitude","Latitude","Year")

groups<- merge(my_groups,meta1,by.x="fid",by.y="Sample.Name")
groups
groups$merged <- paste(groups$fid, groups$iid, sep = "/")
colnames(groups) <- c("fid","iid","pid","Country","Population","Longitude","Latitude","Year","ID")



vertices<-V(g)$name
vertices

## identify samples not in the graph (not clustered) and add them
unique_elements <- setdiff(samples, vertices)
length(unique_elements)
g <- add_vertices(g, length(unique_elements), attr = list(name = unique_elements))
vertices<-V(g)$name



groups_s <- groups[match(vertices, groups$ID), ]

#### parse through variants and recode based on nomenclature Mueller et al.
var <- data.frame(Sample.Name=groups_s$fid)
var$var <-""
for (i in 1:nrow(var)){
  flag=0
  for (z in 1:nrow(variants)){
    if (var$Sample.Name[i] == variants$Sample.Name[z]){
      flag = flag+1
      if (flag > 1){
        var$var[i] = paste0(var$var[i],"/",variants$variants[z])
      }
      if (flag ==1){
        var$var[i] = variants$variants[z]
      }
    }
  }
}

var

for (i in 1:nrow(var)){
  if (var$var[i] == ""){var$var[i]  <- "NA"}
  }
var$var <- gsub("0/0", "B",var$var)
var$var <- gsub("1/1", "C",var$var)
var$var <- gsub("2/2", "A",var$var)
var$var <- gsub("0/1", "B/C",var$var)
var$var <- gsub("0/2", "A/B",var$var)
var$var <- gsub("1/2", "A/C",var$var)
var$var <- gsub("1/4", "C/I",var$var)
var$var <- gsub("2/6", "A/K",var$var)
var$var <- gsub("2/7", "A/F",var$var)
var$var <- gsub("1/8", "C/L",var$var)
var$var <- gsub("1/5", "C/J",var$var)
var$var <- gsub("0", "B",var$var)
var$var <- gsub("1", "C",var$var)
var$var <- gsub("2", "A",var$var)
var$var <- gsub("3", "H",var$var)



# var
# 0 = varB
# 1 = varC
# 2 = varA
# 3 = varH  Y31H + A53V compared to varA
# 4 = varI  Y31H + A53V + R80S compared to varA
# 5 = varJ  R80S
# 6 = varK E55R + G61A
# 7 = varF
# 8 = varL E55R + G61A + R80S



#### parse through dna_variants
dna_var <- data.frame(Sample.Name=groups_s$fid)
dna_var$dna_var <-""
for (i in 1:nrow(dna_var)){
  flag=0
  for (z in 1:nrow(dna_variants)){
    if (dna_var$Sample.Name[i] == dna_variants$Sample.Name[z]){
      flag = flag+1
      if (flag > 1){
        dna_var$dna_var[i] = paste0(dna_var$dna_var[i],"/",dna_variants$dna_variants[z])
      }
      if (flag ==1){
        dna_var$dna_var[i] = dna_variants$dna_variants[z]
      }
    }
  }
}


for (i in 1:nrow(dna_var)){
  if (dna_var$dna_var[i] == ""){dna_var$dna_var[i]  <- "NA"}
}

unique(dna_var$dna_var)


# merge everything together and reorder
groups_sor<- merge(groups_s,var,by.x="fid",by.y="Sample.Name")
groups_sort<- merge(groups_sor,dna_var,by.x="fid",by.y="Sample.Name")
groups_sort_p <- merge(groups_sort,phenotype,by.x="fid",by.y="Sample.Name",all.x=TRUE)

vertices<-V(g)$name
groups_sorted <- groups_sort_p[match(vertices, groups_sort_p$ID), ]



#plot distr of variants

groups_sort_clean <- groups_sort[!is.na(groups_sort_p$var), ]


color_mapping_var <- c("A" = "#76EEC6", "B" = "#556B2F","C" = "#8B4500", "H" = "#EE9A49", "A/B" = "#EE3A8C", "A/C" = "#008B8B","B/C" = "#458B00")#,"A/K" = "purple","A/F" = "violet","C/L" = "gold", "NA" = "gray", "C/I" = "aquamarine","C/J" = "antiquewhite")

ggplot(groups_sort_clean, aes(x = var, fill = var)) +
  geom_bar() +
  scale_fill_manual(values = color_mapping_var)+
    theme_minimal()+
  labs(x= "Variants", y= "Number of isolates")+
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(fill = NA))  # Keep the box around the plot


ggsave("Avrpm17_variants_stack_bplot.pdf",width=10,height=10)


proportions <- groups_sort_clean %>%
  count(var) %>%
  mutate(proportion = n / sum(n))

ggplot(proportions, aes(x = var, y = proportion, fill = var)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping_var)+
  theme_minimal()+
  labs(x= "Variants", y= "Proportion of isolates")+
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(fill = NA))  # Keep the box around the plot

ggsave("Avrpm17_variants_proportion_stack_bplot.pdf",width=10,height=10)

#for (i in 1:nrow(groups_sorted)){
#  if (groups_sorted$Year[i] < 2000){groups_sorted$period[i] = "before 2000"}
#  else {groups_sorted$period[i] = "after 2000"}
#}


######## plot all clusters
## attach attributes to graph
vertex_attr(g, "Population", index = V(g)) <- groups_sorted$Population
vertex_attr(g, "Variant", index = V(g)) <- groups_sorted$var
vertex_attr(g, "dna_Variant", index = V(g)) <- groups_sorted$dna_var

#vertex_attr(g, "Period", index = V(g)) <- groups_sorted$period

#Define a color mapping
#color_mapping_period <- c("before 2000" = "blue", "after 2000" = "red")
color_mapping_pop <- c("N_EUR" = "#377EB8", "S_EUR1" = "#EA9999","S_EUR2" = "#E41A1C", "TUR" = "#E5B110", "ME" = "#984EA3")
#Assign colors to the vertices based on their population
vertex_colors_pop <- color_mapping_pop[V(g)$Population]
vertex_colors_var <- color_mapping_var[V(g)$Variant]
#vertex_colors_period <- color_mapping_period[V(g)$Period]

pdf("Avrpm17_ibd_all_clusters_main.pdf")
par(oma = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

par(mfrow = c(1, 2))
layout= layout_nicely(g)

#plot(g,
#     layout=layout,
#     #vertex.shapes=
#     vertex.size=2.5,
#     vertex.label.cex=0.5,
#     vertex.label.dist=0.4,
#     vertex.label.color="black",
#     vertex.label=NA,
#     vertex.color=vertex_colors_period
#     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Variant"))
#)
#legend("bottomleft", legend = names(color_mapping_period), fill = color_mapping_period, cex = 0.5)

par(mar = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

plot(g,
     layout=layout,
     #vertex.shapes=
     vertex.size=2.5,
     vertex.label.cex=0.5,
     vertex.label.dist=0.4,
     vertex.label.color="black",
     vertex.label=NA,
     vertex.frame.width=0.5,
     edge.width=0.3,
     vertex.color=vertex_colors_pop
     #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Variant"))
)
legend(x=-1,y=-1, legend = names(color_mapping_pop), fill = color_mapping_pop, cex = 0.6, bty="n",horiz= "True")

par(mar = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

plot(g,
     layout=layout,
     #vertex.shapes=
     vertex.size=2.5,
     vertex.label.cex=0.5,
     vertex.label.dist=0.4,
     vertex.label.color="black",
     vertex.label=NA,
     vertex.frame.width=0.5,
     edge.width=0.3,
     vertex.color=vertex_colors_var

     #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Variant"))
)
legend(x=-1,y=-1, legend = names(color_mapping_var), fill = color_mapping_var, cex = 0.6,bty="n",ncol=6)

dev.off()


pdf("Avrpm17_ibd_all_clusters_supplementary.pdf")
par(oma = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

par(mfrow = c(1, 3))

#plot(g,
#     layout=layout,
#     #vertex.shapes=
#     vertex.size=2.5,
#     vertex.label.cex=0.5,
#     vertex.label.dist=0.4,
#     vertex.label.color="black",
#     vertex.label=NA,
#     vertex.color=vertex_colors_period
#     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Variant"))
#)
#legend("bottomleft", legend = names(color_mapping_period), fill = color_mapping_period, cex = 0.5)

par(mar = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

plot(g,
     layout=layout,
     #vertex.shapes=
     vertex.size=2.5,
     vertex.label.cex=0.5,
     vertex.label.dist=0.4,
     vertex.label.color="black",
     vertex.label=NA,
     vertex.frame.width=0.5,
     edge.width=0.3,
     vertex.color=vertex_colors_pop
     #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Variant"))
)
legend(x=-1,y=-1, legend = names(color_mapping_pop), fill = color_mapping_pop, cex = 0.6, bty="n",horiz= "True")

par(mar = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

plot(g,
     layout=layout,
     #vertex.shapes=
     vertex.size=2.5,
     vertex.label.cex=0.5,
     vertex.label.dist=0.4,
     vertex.label.color="black",
     vertex.label=NA,
     vertex.frame.width=0.5,
     edge.width=0.3,
     vertex.color=vertex_colors_var
     
     #     vertex.color = as.factor(vertex_attr(my_i_clusters$i.network, "Variant"))
)
legend(x=-1,y=-1, legend = names(color_mapping_var), fill = color_mapping_var, cex = 0.6,bty="n",ncol=6)

palette <- colorRampPalette(brewer.pal(min(length(unique(dna_var$dna_var)), 9), "Set1"))(length(unique(dna_var$dna_var)))
vertex_colors <- palette[as.numeric(factor(V(g)$dna_Variant, levels = unique(dna_var$dna_var)))]
plot(g,
     layout=layout,
     #vertex.shapes=
     vertex.size=2.5,
     vertex.label.cex=0.5,
     vertex.label.dist=0.4,
     vertex.label.color="black",
     vertex.label=NA,
     vertex.frame.width=0.5,
     edge.width=0.3,
     vertex.color=vertex_colors
     
#     vertex.color = as.factor(vertex_attr(g, "dna_Variant"))
)

legend(x=-1,y=-1, legend = unique(dna_var$dna_var), fill = palette, cex = 0.6,bty="n",ncol=6)


par(mar = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)


dev.off()




## plot one cluster at the time, graoh and map
world <- map_data("world")
i=1

for (i in 1:10){

  # identify vertices in cluster
subcluster_vertices <- which(V(g)$name %in% my_i_clusters$clusters[[i]])

# Create a subgraph containing only the vertices of the subcluster
subgraph <- induced_subgraph(g, subcluster_vertices)

vertices<-V(subgraph)$name


# fetch metainformations and make attributes for subgraph
sub<-groups_sorted[groups_sorted$ID %in% vertices,]
ver_to_exclude=NULL
ver_to_exclude<-subset(sub,is.na(sub$Latitude))$ID

subgraph <- delete_vertices(subgraph, ver_to_exclude)
vertices<-V(subgraph)$name


sub_sorted <- sub[match(vertices, sub$ID), ]


#vertex_attr(subgraph, "Population", index = V(subgraph)) <- sub_sorted$Population
vertex_attr(subgraph, "Longitude", index = V(subgraph)) <- sub_sorted$Longitude
vertex_attr(subgraph, "Latitude", index = V(subgraph)) <- sub_sorted$Latitude
vertex_attr(subgraph, "Variant", index = V(subgraph)) <- sub_sorted$var



# Step 2: Assign colors to the vertices based on their attributes
vertex_colors_var <- color_mapping_var[V(subgraph)$Variant]

par(mar = c(4,0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)
par(oma = c(2, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

pdf(paste0("Avrpm17_ibd_cluster_",i,".pdf"))

plot(subgraph,
     edge.width=0.2,
     vertex.size=2.5,
     vertex.label.cex=0.5,
     vertex.label.dist=0.4,
     vertex.label.color="black",
     vertex.label=NA,
     vertex.frame.width=0.5,
     edge.width=0.3,
     vertex.color=vertex_colors_var)

legend(x=-1,y=1, legend = names(color_mapping_var), fill = color_mapping_var, cex = 0.6, bty="n",ncol=6)

dev.off()



## plot on map with ggplot
# layout for map plotting


lo <- layout.norm(as.matrix(sub_sorted[,6:7]))
#convert nodes and edges to df
nodes <- data.frame(
  id = V(subgraph)$name,
  lat = V(subgraph)$Latitude,
  lon = V(subgraph)$Longitude
)

nodes$color <- color_mapping_var[V(subgraph)$Variant]


edges <- as.data.frame(get.edgelist(subgraph))
colnames(edges) <- c("from", "to")
edges <- edges %>%
  left_join(nodes, by = c("from" = "id")) %>%
  rename(lat_from = lat, lon_from = lon) %>%
  left_join(nodes, by = c("to" = "id")) %>%
  rename(lat_to = lat, lon_to = lon)


### make plot with nodes

p <- ggplot() +
#  geom_point(data = nodes, aes(x = lon, y = lat), color = vertex_colors, size = 2.5)
  geom_jitter(data = nodes, aes(x = lon, y = lat),color=nodes$color, size = 1.5,width=0.6, height=0.6,show.legend=TRUE)+
  geom_point(data = nodes, aes(x = lon, y = lat), color = "black", size = 1,show.legend=TRUE)+
  scale_color_manual(name='Regression Model',
                     breaks=c('A', 'B', 'C'),
                     values=c('A'='grey', 'B'='blue', 'C'='red'))
  
##add edges
p <- p +
  geom_segment(data = edges, aes(x = lon_from, y = lat_from, xend = lon_to, yend = lat_to), color = "gray50",linewidth=0.15,alpha=0.4)
## add map
p <- p +  geom_path(data = world, aes(x = long, y = lat, group = group), color = "gray50",linewidth=0.2)+
  xlim(c(-12, 45)) +  # Set specific x-axis limits
  ylim(c(20, 65))+
  coord_fixed(ratio = 1)+
  theme_void()  # Ins


p

ggsave(paste0("Avrpm17_ibd_cluster_",i,"_map.pdf"),width=10,height=10)




}




### plot phenotype by variant

groups_sort_p_clean <- groups_sort_p[!is.na(groups_sort_p$Amigo), ]


ggplot(groups_sort_p_clean, aes(x = var, fill = factor(Amigo))) +
  geom_bar(position = "stack") +
  labs(x = "Protein variant", y = "Count", fill = "Phenotype on Amigo") +
  scale_fill_discrete(labels = c("Avirulent", "Intermediate", "Virulent"))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(fill = NA))  # Keep the box around the plot




ggsave(paste0("Avrpm17_pheno_stack_bplot.pdf"),width=10,height=10)



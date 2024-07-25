library(igraph)


#1
df <- read.csv('facebook_combined.txt', sep=' ')

g <- graph_from_data_frame(df, directed=F, vertices=NULL)

plot(g, main = "Plot of the Facebook Network")

is_connected(g)

#2
diameter(g)

#3
deg_dist3 <- degree_distribution((g))
plot(deg_dist3, main = "Degree Distribution of the Facebook Network", xlab = "Degree", ylab = "Distribution")

mean(degree(g))

#4

clean_g_3 <- which(deg_dist3 != 0, arr.ind=TRUE)
x3 <- log(seq(1:length(deg_dist3)))[clean_g_3]
y3 <- log(deg_dist3)[clean_g_3]
lr_3 <- lm(y3~x3)
cf_3 <- coef(lr_3)
print(cf_3[2])
plot(x3, y3, abline(lr_3, col="red"), main = "Log Scale Degree Distribution", xlab = "log(Degree)", ylab = "log(Distribution")

#5
g_id <- neighbors(g, 1)
g_pers <- induced_subgraph(g, c(1, g_id))

print(vcount(g_pers))
print(ecount(g_pers))

#6
diameter(g_pers)

#8
count <- 0
deg_sum <- 0
for (i in c(1:(vcount(g)))) {
  #print(i)
  if (length(neighbors(g, i))>200) {
    count <- count + 1
    deg_sum <- deg_sum + degree(g, i)
  }
}
deg_avg <- deg_sum/count
print(count)
print(deg_avg)


#9
cnl <- c(0, 107, 348, 483, 1086)
#cnl <- c(107, 348)
g_saved <- vector(mode='list', length=5)
counter <- 1

for (node in cnl) {
  node_id <- match(node, V(g)$name)
  g_core <- induced_subgraph(g, c(node_id, neighbors(g, node_id)))
  print(paste("Node ", node, " match ID ", node_id, " Degree ", degree(g, node_id)))
  fast <- cluster_fast_greedy(g_core)
  between <- cluster_edge_betweenness(g_core, directed=F)
  info <- cluster_infomap(g_core)
  
  print(paste("Fast ", modularity(fast), " Between: ", modularity(between), " Info: ", modularity(info)))
  
  plot(fast, g_core, main = paste("Community Structure of Core Node ID ", node+1, " with Fast-Greedy"))
  plot(between, g_core, main = paste("Community Structure of Core Node ID ", node+1, " with Edge-Betweenness"))  
  plot(info, g_core, main = paste("Community Structure of Core Node ID ", node+1, " with Infomap"))
  g_saved[[i]] <- g_core
  counter <- counter + 1
}

#10

for (node in cnl) {
  node_id <- match(node, V(g)$name)
  g_core <- induced_subgraph(g, c(neighbors(g, node_id)))
  print(paste("Node ", node, " match ID ", node_id, " Degree ", degree(g, node_id)))
  fast <- cluster_fast_greedy(g_core)
  between <- cluster_edge_betweenness(g_core, directed=F)
  info <- cluster_infomap(g_core)
  
  print(paste("Fast ", modularity(fast), " Between: ", modularity(between), " Info: ", modularity(info)))
  
  plot(fast, g_core, main = paste("Community Structure of Core Node ID ", node+1, " with Fast-Greedy"))
  plot(between, g_core, main = paste("Community Structure of Core Node ID ", node+1, " with Edge-Betweenness"))  
  plot(info, g_core, main = paste("Community Structure of Core Node ID ", node+1, " with Infomap"))
}

#12
#cnl <- c(0, 107, 348, 483, 1086)

cnl2 <- c(match(0, V(g)$name), match(107, V(g)$name), match(348, V(g)$name), match(483, V(g)$name), match(1086, V(g)$name))
cnl4 <- cnl2
#cnl3 <- c(1, 108, 349, 484, 1087)
g2 <- g
V(g2)$name <- V(g2)


for(node in cnl) {
  node_id <- match(node, V(g)$name)
  g_core <- induced_subgraph(g, c(neighbors(g, node_id)))
  g_core_deg <- degree_distribution(g_core)
  
  plot(g_core_deg, type="h", main=paste("Embeddedness of each neighbor for Core Node ID ", node+1), xlab = "Embeddedness", ylab="Distribution")
}

for(node_id in cnl4) {

  neighb_ego <- unlist(ego(g2, nodes = node_id, mode="all"))
  g_core <- induced_subgraph(g2, vids = as.character(neighb_ego))
  disp <- c()
  rename <- V(g_core)[V(g_core)$name == node_id]
  disp_max <- 0
  disp_max_node <- NA
  
  for (neighbor in neighb_ego) {
    n_dist <- 0
    rename_n <- V(g_core)[V(g_core)$name == neighbor]
    
    if (rename != rename_n) {
      n_neighbors <- unlist(neighbors(g_core, rename))
      neighb_neighbors <- unlist(neighbors(g_core, rename_n))
      n_intersect <- intersect(n_neighbors, neighb_neighbors)
      
      n_disp <- induced_subgraph(g_core, n_intersect)
      
      dists <- distances(n_disp, v=V(n_disp), to=V(n_disp))
      dists[is.infinite(dists)] <- NA
      n_dist <- sum(dists, na.rm = TRUE)/2
    }
    disp <- c(disp, n_dist)
  }
  disp_distribution <- c()
  for (x in unique(disp)) {
    tally <- 0
    for (y in disp) {
      if (x == y) {
        tally <- tally + 1
      }
    }
    disp_distribution <- c(disp_distribution, tally/length(disp))
  }
  plot(x= unique(disp), y=disp_distribution, type="h", main = "Dispersion Distribution", xlab="Dispersion", ylab="Distribution")
}


#13


for (node_id in cnl4) {
  #node_id <- match(node, V(g)$name)
  neighb_ego <- unlist(ego(g2, nodes = node_id, mode="all"))
  g_core <- induced_subgraph(g2, vids = as.character(neighb_ego))
  
  rename <- V(g_core)[V(g_core)$name == node_id]
  disp_max <- 0
  disp_max_node <- NA
  
  for (neighbor in neighb_ego) {
    n_dist <- 0
    rename_n <- V(g_core)[V(g_core)$name == neighbor]
    
    if (rename != rename_n) {
      n_neighbors <- unlist(neighbors(g_core, rename))
      neighb_neighbors <- unlist(neighbors(g_core, rename_n))
      n_intersect <- intersect(n_neighbors, neighb_neighbors)
      
      n_disp <- induced_subgraph(g_core, n_intersect)
      
      dists <- distances(n_disp, v=V(n_disp), to=V(n_disp))
      dists[is.infinite(dists)] <- NA
      n_dist <- sum(dists, na.rm = TRUE)/2
    }
    if (n_dist > disp_max) {
      disp_max <- n_dist
      disp_max_node <- neighbor
    }
  }
  fast <- cluster_fast_greedy(g_core)
  node_colors <- fast$membership + 1
  node_colors[which(V(g_core)$name == disp_max_node)] <- "red"
  V(g_core)$color <- fast$membership + 1
  V(g_core)[which(V(g_core)$name == disp_max_node)]$color <- "red"
  
  node_sizes <- rep(4, length(node_colors))
  node_sizes[which(V(g_core)$name == disp_max_node)] <- 8
  
  edge_colors <- rep("light blue", length(E(g_core)))
  edge_colors[which(get.edgelist(g_core)[,1] == disp_max_node | get.edgelist(g_core)[,2] == disp_max_node)] <- "red"
  
  edge_thickenss <- rep(0.5, length(E(g_core)))
  edge_thickness[which(get.edgelist(g_core)[,1] == disp_max_node | get.edgelist(g_core)[,2] == disp_max_node)] <- 2
  
  plot(g_core, edge.color = edge_colors, vertex.size = node_sizes, edge.width = edge_thickness, mark.border = NA, mark.col = NA, main = paste("Highlighted Personalized Network for Core Node ", node_id), vertex.label = NA)
  
  
}


#14
for (node_id in cnl4) {
  
  max_embed <- 0
  max_embed_id <- NA
  iter <- 1
  
  max_d_e <- 0
  d_e_id <- NA
  ###############
  neighb_ego <- unlist(ego(g2, nodes = node_id, mode="all"))
  g_core <- induced_subgraph(g2, vids = as.character(neighb_ego))
  
  rename <- V(g_core)[V(g_core)$name == node_id]
  disp_list <- c()
  disp_names <- c(V(g_core)[V(g_core)$name != node_id]$name)
  
  for (neighbor in neighb_ego) {

    n_dist <- 0
    n_embed <- 0
    rename_n <- V(g_core)[V(g_core)$name == neighbor]

    iter <- iter + 1
    if (rename != rename_n) {
      n_neighbors <- unlist(neighbors(g_core, rename))
      neighb_neighbors <- unlist(neighbors(g_core, rename_n))
      n_intersect <- intersect(n_neighbors, neighb_neighbors)
      n_embed <- length(n_intersect)
      
      n_disp <- induced_subgraph(g_core, n_intersect)
      
      dists <- distances(n_disp, v=V(n_disp), to=V(n_disp))
      dists[is.infinite(dists)] <- NA
      n_dist <- sum(dists, na.rm = TRUE)/2
      
      if (n_embed > max_embed) {
        max_embed <- n_embed
        max_embed_id <- neighbor
      }
      if (n_embed == 0) {
        div <- 0
      }
      else {
        div <- n_dist / n_embed
      }
      #print(div)
      if (div > max_d_e) {
        max_d_e <- div
        d_e_id <- neighbor
      }

    }
  }
  
  
  fast <- cluster_fast_greedy(g_core)
  node_colors <- fast$membership + 1
  node_colors[which(V(g_core)$name == max_embed_id)] <- "red"
  V(g_core)$color <- fast$membership + 1
  V(g_core)[which(V(g_core)$name == max_embed_id)]$color <- "red"
  node_colors[which(V(g_core)$name == d_e_id)] <- "green"
  V(g_core)[which(V(g_core)$name == d_e_id)]$color <- "green"
  
  
  node_sizes <- rep(4, length(node_colors))
  node_sizes[which(V(g_core)$name == max_embed_id)] <- 8
  node_sizes[which(V(g_core)$name == d_e_id)] <- 8
  
  edge_colors <- rep("light blue", length(E(g_core)))
  edge_colors[which(get.edgelist(g_core)[,1] == max_embed_id | get.edgelist(g_core)[,2] == max_embed_id)] <- "red"
  edge_colors[which(get.edgelist(g_core)[,1] == d_e_id | get.edgelist(g_core)[,2] == d_e_id)] <- "green"
  
  edge_thickness <- rep(0.5, length(E(g_core)))
  edge_thickness[which(get.edgelist(g_core)[,1] == max_embed_id | get.edgelist(g_core)[,2] == max_embed_id)] <- 3
  edge_thickness[which(get.edgelist(g_core)[,1] == d_e_id | get.edgelist(g_core)[,2] == d_e_id)] <- 3
  
  plot(g_core, edge.color = edge_colors, vertex.size = node_sizes, edge.width = edge_thickness, mark.border = NA, mark.col = NA, main = paste("Highlighted Personalized Network for Core Node ", node_id), vertex.label = NA)
  print(paste("Max embed ID: ", max_embed_id, " Max e/d ID: ", d_e_id))
  
}


#16

nr_list <- c()
for (nr_node in V(g)$name) {
  nr_neighbors <- neighbors(g, nr_node)
  if (length(nr_neighbors) == 24) {
    nr_list <- c(nr_list, nr_node)
  }
}
print(length(nr_list))


#17
#common neighbors

common_neighbors = function(g, node_i, node_j) {
  cn_list <- vector(mode='list', length=4)
  neighbors_i <- neighbors(g, node_i)
  cn_list[[1]] <- neighbors_i
  neighbors_j <- neighbors(g, node_j)
  cn_list[[2]] <- neighbors_j
  intersect_ij <- intersection(neighbors_i, neighbors_j)
  cn_list[[3]] <- intersect_ij
  score <- length(intersect_ij)
  cn_list[[4]] <- score
  return (cn_list)
}


#jaccard
jaccard = function(g, node_i, node_j) {
  #cn <- common_neighbors(g, node_i, node_j)
  #si <- cn[[1]]
  #sj <- cn[[2]]
  #i_ij_length <- cn[[4]]
  
  #print(paste("si ", si, " sj ", sj, " inter ", i_ij_length))
  #union_score <- length(si) + length(sj) - i_ij_length
  #print(length(cn[1]))
  #return (i_ij_length/union_score)
  return (similarity(g, c(node_i, node_j))[1,2])
}


#adamic adar
adamic_adar = function(g, node_i, node_j) {
  cn <- common_neighbors(g, node_i, node_j)
  inter_list <- cn[[3]]
  ad_ad_sum <- 0
  for (node in inter_list) {
    temp <- 1/(log(length(neighbors(g, match(node, V(g)$name)))))
    #print(length(neighbors(g, match(node, V(g)$name))))
    ad_ad_sum <- ad_ad_sum + temp
  }
  return (ad_ad_sum)
}

##

n415 <- match(415, V(g)$name)
neighbors_415 <- neighbors(g, n415)
cn_acc <- c()
j_acc <- c()
aa_acc <- c()


nr_length <- length(nr_list)

delete_later <- 1

for (n_node in nr_list) {
  cn_score_list <- c()
  #print(paste("n_ node ", n_node))
  for (i in c(1:1)) {
    reduced_graph <- g
    ri <- c()
    #print(paste("i ", i))
    for (neigh415 in neighbors_415) {
      if (runif(1) <= 0.25) {
        reduced_graph <- delete_edges(reduced_graph, c(n415, neigh415))
        ri <- c(ri, neigh415)
      }
    }
    ri_length <- length(ri)
    
    dists <- distances(reduced_graph, n415, V(reduced_graph)$name)
    non_neighbors <- which(dists==2)
    non_length <- length(non_neighbors)
    
    #new_neigh415 <- neighbors(reduced_graph, n415)
    #new_all_nodes <- V(reduced_graph)$name
    #non_neighbors <- new_all_nodes[! new_all_nodes %in% new_neigh415]
    #non_neighbors <- non_neighbors[non_neighbors != n415]
    
    #non_length <- length(non_neighbors)
    #print(paste("non_length ", non_length))
    
    cn_list <- matrix(NA, nrow=non_length, ncol=2)
    
    cn_friends <- c()
    
    cn_inter <- c()
    
    iter <- 1
    
    for (non_node in non_neighbors) {
      #print(paste("non_node ", non_node))
      #print(delete_later)
      #delete_later <- delete_later + 1
      temp <- common_neighbors(reduced_graph, n415, non_node)
      cn_list[iter, 1] <- non_node
      cn_list[iter, 2] <- temp[[4]]
      iter <- iter + 1
    }
    
    sorted_cn <- cn_list[order(cn_list[,2]),]
    
    for (j in c(1:ri_length)) {
      cn_friends <- c(cn_friends, sorted_cn[j, 1])
    }
    
    cn_inter <- intersect(ri, cn_friends)
    print(ri)
    print(cn_friends)
    print(cn_inter)
    
    cn_score <- length(cn_inter)/ri_length
    
    cn_score_list <- c(cn_score_list, cn_score)
  }
  cn_acc <- c(cn_acc, mean(cn_score_list))
}
print(paste("Common neighbors ", mean(cn_acc)))
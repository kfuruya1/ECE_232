library(igraph)

ma <- as.matrix(read.table("facebook_combined.txt", sep = " ")) +1
g = graph_from_edgelist(ma, directed = F)

#16
g_415 <- induced_subgraph(g, c(415, neighbors(g, 415)))

degree_415 <- degree(g_415)
nr_list <- which(degree_415 == 24)
print(length(nr_list))


#17
#common neighbors

common_neighbors = function(g, node_i, node_j) {

  neighbors_i <- neighbors(g, node_i)

  neighbors_j <- neighbors(g, as.character(node_j))

  intersect_ij <- intersect(neighbors_i, neighbors_j)

  score <- length(intersect_ij)

  return (score)
}


#jaccard
jaccard = function(g, node_i, node_j) {

  return (similarity(g, c(node_i, node_j))[1,2])
}

#adamic adar
adamic_adar = function(g, node_i, node_j) {
  neighbors_i <- neighbors(g, node_i)
  neighbors_j <- neighbors(g, as.character(node_j))
  intersect_ij <- intersect(neighbors_i, neighbors_j)
  ad_ad_sum <- 0
  for (node in intersect_ij) {
    temp <- 1/(log(length(neighbors(g, node))))

    ad_ad_sum <- ad_ad_sum + temp
  }
  return (ad_ad_sum)
}

#######################################
neighbors_415 <- neighbors(g, 415)


nr_length <- length(nr_list)

cn_node_acc <- c()
j_node_acc <- c()
aa_node_acc <- c()

for (n_node in nr_list) {
  cn_acc <- c()
  j_acc <- c()
  aa_acc <- c()
  
  for (i in c(1:10)) {
    reduced_graph <- g_415
    ri <- c()
    n_node_neighbors <- neighbors(reduced_graph, n_node)
    for (n_neigh in n_node_neighbors) {
      if (runif(1) <= 0.25) {
        reduced_graph <- delete_edges(reduced_graph, c(n_node, n_neigh))
        ri <- c(ri, n_neigh)
      }
    }
    ri_length <- length(ri)
    
    #new_neighbors <- neighbors(reduced_graph, n_node)
    #print(new_neighbors)
    #new_all_nodes <- V(reduced_graph)
    #print(new_all_nodes)
    #non_neighbors <- new_all_nodes[! new_all_nodes %in% new_neighbors]
    #non_neighbors <- non_neighbors[non_neighbors != 415]
    #print(non_neighbors)
    new_neighbors <- setdiff(n_node_neighbors, ri)
    new_neighbors <- append(new_neighbors, n_node)
    non_neighbors <- setdiff(V(reduced_graph), new_neighbors)
    
    
    non_length <- length(non_neighbors)
    #print(non_length)
    
    cn_score <- c()
    j_score <- c()
    aa_score <- c()
    
    for (non_node in non_neighbors) {
      cn_score <- c(cn_score, common_neighbors(reduced_graph, n_node, non_node))
      #print(common_neighbors(reduced_graph, n_node, non_node))
      j_score <- c(j_score, jaccard(reduced_graph, n_node, non_node))
      aa_score <- c(aa_score, adamic_adar(reduced_graph, n_node, non_node))
    }
    #print(cn_score)
    #print(j_score)
    #print(aa_score)
    sorted_cn <- non_neighbors[sort(cn_score, decreasing=T, index.return=T)$ix[1:ri_length]]
    sorted_j <- non_neighbors[sort(j_score, decreasing=T, index.return=T)$ix[1:ri_length]]
    sorted_aa <- non_neighbors[sort(aa_score, decreasing=T, index.return=T)$ix[1:ri_length]]
    
    #print(sorted_cn)
    #print(sorted_j)
    #print(sorted_aa)
    #cn_rec <- sorted_cn$ix[1:ri_length]
    #j_rec <- sorted_j$ix[1:ri_length]
    #aa_rec <- sorted_aa$ix[1:ri_length]
    #print(ri)
    cn_intersect <- intersect(sorted_cn, ri)
    j_intersect <- intersect(sorted_j, ri)
    aa_intersect <- intersect(sorted_aa, ri)
    #print(cn_intersect)
    #print(j_intersect)
    #print(aa_intersect)
    
    cn_acc <- c(cn_acc, length(cn_intersect)/ri_length)
    j_acc <- c(j_acc, length(j_intersect)/ri_length)
    aa_acc <- c(aa_acc, length(aa_intersect)/ri_length)
    #print(cn_acc)
  }
  cn_node_acc <- c(cn_node_acc, mean(cn_acc))
  j_node_acc <- c(j_node_acc, mean(j_acc))
  aa_node_acc <- c(aa_node_acc, mean(aa_acc))
  #print(cn_node_acc)
}
final_cn <- mean(cn_node_acc)
final_j <- mean(j_node_acc)
final_aa <- mean(aa_node_acc)

print(paste("cn ", final_cn, " j ", final_j, " aa ", final_aa))

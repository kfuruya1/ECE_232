library("igraph")

file_path ="gplus/" 
edge_files = list.files(path=file_path,pattern="edges") 
circles_files = list.files(path=file_path,pattern="circles")
fts_files = list.files(path=file_path,pattern="feat")
initial_graph = list()
final_graph = list()
graph_circles = list()
ego_nodes = list()


cnt = 0
node_names = c()
for(i in 1:length(edge_files)){
  # get node id
  node = strsplit(edge_files[i],".edges")[[1]]
  node_names <-c(node_names,node)
  #print(node)
  ego_nodes[i] = node
  fc = file(paste(file_path,node,".circles",sep=""),open="r") 
  if(length(fc)>0){
    file_lines <- readLines(fc)
    if(length(file_lines)>0){
      circles =list()
      for(j in 1:length(file_lines)){
        circle_users = strsplit(file_lines[j],"\t")
        circles[[j]] <- circle_users[[1]][-1]
      }
      # find users who have more than 2 circles
      if(length(circles)>2){
        #print(i)
        #print(length(circles))
        cnt = cnt + 1
        initial_graph[[i]] <- read.graph(paste(file_path,edge_files[i],sep=""),format="ncol",directed=TRUE)
        graph_circles[[i]] <- circles
        graph_nodes <- V(initial_graph[[i]])
        #print(length(graph_nodes))
        #print(node)
        # add the core node to his neighbor list and construct the graph
        final_graph[[i]] <- add.vertices(initial_graph[[i]],1,name=node)
        core_index = which(V(final_graph[[i]])$name==node) 
        core_node_edges = list()
        ### add edges connecting to this core node
        for(k in 1:length(graph_nodes)){
          core_node_edges = c(core_node_edges, c(core_index, k))
        }
        final_graph[[i]] <- add.edges(final_graph[[i]],core_node_edges)
      } 
    }
  }
  close(fc)
}

print(graph_circles)
#print(graph_circles[[123]])


#q18
cat("there are ", length(edge_files),"nodes and there are ",cnt,"personal networks" )



#q19 
interest_node = c('109327480479767108490', '115625564993990145546','101373961279443806744')
graph_inds = c()
for (i in 1: length(interest_node)){
  graph_ind <- which(node_names==interest_node[i])
  graph_inds <- c(graph_inds, graph_ind)
  #print(graph_ind)
  tmp_graph = final_graph[[graph_ind]]
  hist(degree(tmp_graph, mode="in"),main = paste("In degree for ", interest_node[i]))
  hist(degree(tmp_graph, mode="out"),main = paste("Out degree for ", interest_node[i]))
}

#20
comm_list <- list()
for (i in 1: length(interest_node)){
  graph_ind <- which(node_names==interest_node[i])
  graph_inds <- c(graph_inds, graph_ind)
  #print(graph_ind)
  tmp_graph = final_graph[[graph_ind]]
  
  walk_comm <- cluster_walktrap(tmp_graph, modularity=T)
  comm_list[[i]] <- walk_comm
  
  mod <- modularity(walk_comm)
  #print(mod)
  
  plot(walk_comm, tmp_graph, vertex.label=NA, vertex.size = 5)
  
}
print(length(comm_list))


#21

find_n = function(node) {
  graph_ind <- which(node_names==node)
  circle_list <- graph_circles[[graph_ind]]
  n_length <- c()
  
  for(i in 1:length(circle_list)) {
    n_length <- c(n_length, circle_list[[i]])
  }
  return(length(n_length))
}

circle_entropy = function (node, n) {
  e_sum <- 0
  graph_ind <- which(node_names==node)
  #print(graph_ind)
  #a <- 0
  
  circle_list <- graph_circles[[graph_ind]]
  
  for (circle in circle_list){
    a_i <- length(circle)
    a_i_n <- a_i/n
    e_sum <- e_sum - (a_i_n * log(a_i_n))
    #a <- a + a_i
  }
  #print(a)
  return(e_sum)
}

k_entropy = function(node, community_info, n) {
  k_sum <- 0
  graph_ind <- which(node_names==node)
  circle_list <- graph_circles[[graph_ind]]
  #k_graph <- final_graph[[graph_ind]]
  #V(k_graph)$community <- membership(community_info)
  
  
  #V(final_graph[[graph_ind]])$community <- membership(temp_com)
  #temp <- final_graph[[graph_ind]]
  #print(temp[which(V(temp)$community == 1)])
  
  #a <- 0
  
  circle_nodes <- c()
  
  for (i in 1:length(circle_list)){
    circle_nodes <- c(circle_nodes, circle_list[[i]])
  }
  #circle_nodes <- circle_nodes
  
  for (j in 1:length(community_info)) {
    #k_i_nodes <- V(k_graph)$name[which(community_info$membership == j)]
    k_i_nodes <- community_info[[j]]
    
    #print(k_i_nodes)
    #print(circle_nodes)
    btemp <- circle_nodes %in% k_i_nodes
    b_i <- length(btemp[which(btemp == T)])
    #b_i <- length(intersect(k_i_nodes, circle_nodes))
    #a <- a + b_i
    #print(paste("b_i", b_i))
    if(b_i != 0) {
    b_i_n <- b_i/n
    #print(paste("b_i/n", b_i_n))
    k_sum <- k_sum - (b_i_n * log(b_i_n))
    #print(paste("math", log(b_i_n)))
    #print(paste("k_sum", k_sum))
    }}
  #print(a)
  return(k_sum)
}

c_cond_k = function(node, community_info, n) {
  graph_ind <- which(node_names==node)
  circle_list <- graph_circles[[graph_ind]]
  #a <- 0
  
  ck_sum <- 0
  circle_nodes <- c()
  for(k in 1:length(circle_list)) {
    circle_nodes <- c(circle_nodes, circle_list[[k]])
  }
  #circle_nodes <- unique(circle_nodes)
  
  for (j in 1:length(community_info)) {
    k_j_nodes <- community_info[[j]]
    btemp <- circle_nodes %in% k_j_nodes
    b_j <- length(btemp[which(btemp == T)])
    #b_j <- length(intersect(k_j_nodes, circle_nodes))
    #a <- a + b_j
    if (b_j != 0) {
    for (i in 1:length(circle_list)) {
      a_ji <- length(intersect(circle_list[[i]], k_j_nodes))
      #a <- a + a_ji
      if(a_ji != 0) {
      a_ji_n <- a_ji/n
      ck_sum <- ck_sum - (a_ji_n * log(a_ji/b_j))
    }}}
  }
  #print(a)
  return (ck_sum)
}

k_cond_c = function(node, community_info, n) {
  graph_ind <- which(node_names == node)
  circle_list <- graph_circles[[graph_ind]]
  circle_nodes <- c()
  kc_sum <- 0
  
  #a <- 0
  
  for(k in 1:length(circle_list)) {
    circle_nodes <- c(circle_nodes, circle_list[[k]])
  }
  
  for (i in 1:length(circle_list)) {
    a_i <- length(circle_list[[i]])
    #print(a_i)
    #print("new circle")
    for (j in 1:length(community_info)) {
      k_j_nodes <- community_info[[j]]
      a_ji <- length(intersect(circle_list[[i]], k_j_nodes))
      #print(length(k_j_nodes))
      #print(a_ji)
      #a <- a + a_ji
      if (a_ji != 0) {
      a_ji_n <- a_ji/n
      #print(a_ji_n)
      #print(a_ji/a_i)
      kc_sum <- kc_sum - (a_ji_n * log(a_ji/a_i))
      #print(kc_sum)
    }}
  }
  #print(a)
  return (kc_sum)
}

h_n = function(node, community_info, n) {
  c_k_cond <- c_cond_k(node, community_info, n)
  print(paste("H(C|K)", c_k_cond))
  h_c <- circle_entropy(node, n)
  print(paste("H(C)", h_c))
  return(1-(c_k_cond/h_c))
}

c_n = function(node, community_info, n) {
  k_c_cond <- k_cond_c(node, community_info, n)
  print(paste("H(K|C)", k_c_cond))
  h_k <- k_entropy(node, community_info, n)
  print(paste("H(K)", h_k))
  return(1-(k_c_cond/h_k))
}

for (n in 1:length(interest_node)) {
  nn <- find_n(interest_node[n])
  print(paste("NN ", nn))
  com <- comm_list[[n]]
  print(paste("Homogeneity", h_n(interest_node[n], com, nn)))
  print(paste("Completeness", c_n(interest_node[n], com, nn)))
}






interest_node = c('109327480479767108490', '115625564993990145546','101373961279443806744')




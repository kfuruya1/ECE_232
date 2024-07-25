library('igraph')
library('Matrix')
library('pracma')
library('resample')
library('gtools')
library('stats')

#1a
g_2_1a <- sample_gnp(n=900, p=0.015)

plot(g_2_1a, main="n=900, p=0.015")


#1b

create_transition_matrix = function (g){
  
  # WARNING: make sure your graph is connected (you might input GCC of your graph)
  
  vs = V(g) #contains all the vertices of g
  n = vcount(g) #the number of vertices in g
  adj = as_adjacency_matrix(g) #matrix indicating if an edge exists between two nodes
  adj[diag(rowSums(adj) == 0)] = 1  # handle if the user is using the function for networks with isolated nodes by creating self-edges
  z = matrix(rowSums(adj, , 1))
  
  transition_matrix = adj / repmat(z, 1, n)  # normalize to get probabilities
  
  return(transition_matrix)
}

random_walk = function (g, num_steps, start_node, transition_matrix = NULL){
  if(is.null(transition_matrix))
    transition_matrix = create_transition_matrix(g)
  end_nodes <- vector(mode='list', length = num_steps)
  v = start_node
  for(i in 1:num_steps){
    #fprintf('Step %d: %d\n', i, v)  # COMMENT THIS
    PMF = transition_matrix[v, ]
    v = sample(1:vcount(g), 1, prob = PMF)
    end_nodes[[i]] <- v
  }
  
  return(end_nodes)
}

random_walker = function (g, t, runs) {
  dist_list <- matrix(NA, nrow=runs, ncol=t) 
  deg_list <- array(1:runs)
  node_list <- vector(mode='list', length = t)
  temp <- vector(mode='list', length=t)
  for (i in 1:runs) {
    start_node = sample(1:vcount(g), 1)
    node_list <- random_walk(g, num_steps = t, start=start_node)
    j <- 1
    assign("temp", NULL)
    for (node in node_list){
      temp[j] <- distances(g, start_node, to=node)
      j <- j + 1
    }
    dist_list[i,] <- temp
    last_node <- node_list[t]
    deg_list[i] <- last_node #degree(g, last_node)
  }

  plot(x=1:t, y=colMeans(dist_list), main = "<s(t)> vs t", xlab = "steps", ylab = "<s(t)>")
  plot(x=1:t, y=colVars(dist_list), main = "<(s(t) - <s(t)>)^2>", xlab = "steps", ylab = "<(s(t) - <s(t)>)^2>")
  return(deg_list)
}

deg_list <- random_walker(g_2_1a, t=100, runs=1000)

#1c
plot(degree_distribution(g_2_1a), main="Degree Distribution", xlab = "Degrees", ylab="Distribution")
plot(degree_distribution(g_2_1a, v=deg_list), main="End Node Degree Distribution", xlab ="Degrees", ylab="Distribution")

#1d
g_2_1d <- sample_gnp(n=9000, p = 0.015)
random_walker(g_2_1d, t=100, runs=1000)
diameter(g_2_1a)
diameter(g_2_1d)

###################################

#2a

g_2_2a <- sample_pa(n=900, m=1, directed=F)

plot(g_2_2a, main = "Preferential Attachment Network with n=900 and m=1")

#2b

b2 <- random_walker(g_2_2a, t=100, runs=1000)

#2c
plot(degree_distribution(g_2_2a), main="Degree Distribution for n=900 and m=1", xlab = "Degrees", ylab = "Distribution")
plot(degree_distribution(g_2_2a, v=b2), main = "Degree Distribution for End Nodes", xlab = "Degrees", ylab = "Distribution")

#2d

g_2_2d <- sample_pa(n=90, m=1, directed=F)
g_2_2d9000 <- sample_pa(n=9000, m=1, directed=F)

d90 <- random_walker(g_2_2d, t=100, runs=1000)
d9000 <- random_walker(g_2_2d9000, t=100, runs=1000)
diameter(g_2_2d)
diameter(g_2_2d9000)

########################################################################

#3a

g_2_3a <- sample_pa(n=900, m=4, directed=T)
g_2_3a_2 <- sample_pa(n=900, m=4, directed=T)

plot(g_2_3a, main = "Directed network of n=900 and m=4 prior to adding additional edges")

edge_list2 <- as_edgelist(g_2_3a_2)

scramble <- edge_list2[sample(nrow(edge_list2)),]

index_1 <- permute(scramble[,1])
index_2 <- permute(scramble[,2])

scramble2 <- cbind(index_1, index_2)

p <- add_edges(graph=g_2_3a, edges=scramble2)#permute(edge_list2)) #scramble)


plot(p, main = "Directed network of n=900 and m=4 after adding edges")

deg_dist_3 <- random_walker(p, t=100, runs=1000)

node_prob = function (deg_dist_3) {
  d_count <- numeric(900)
  d_sum <- 0
  d_prob <- array(1:900)
  for (i in 1:1000) {
    temp3 <- deg_dist_3[[i]]
    d_count[[temp3]] <- d_count[[temp3]] + 1
    d_sum <- d_sum + temp3
  }
  for(i in 1:900) {
    d_prob[[i]] <- d_count[[i]]/900
  }
  return (d_prob)
}

d_prob <- node_prob(deg_dist_3)

plot(x=1:900, y=d_prob, main = "Probability of Walker Visiting a Node", xlab = "Node Index", ylab="Probability")

cust_degree = function(nodes) {
  clean_nodes <- numeric(900)
  for (i in 1:1000) {
    num <- nodes[[i]]
    clean_nodes[[num]] <- degree(p, num)
  }
  return(clean_nodes)
}

d_degrees <- cust_degree(deg_dist_3)

plot(x=unlist(d_degrees), y=d_prob, main="Node Degree vs. Visit Probability", xlab="Degree", ylab="Visit Probability")

#3b

cust_rand_walk = function (g, num_steps, start_node, t_prob = 0, transition_matrix = NULL, pg_type = "1/N"){
  if(is.null(transition_matrix))
    transition_matrix = create_transition_matrix(g)
  end_nodes <- vector(mode='list', length = num_steps)
  v = start_node
  for(i in 1:num_steps){
    #fprintf('Step %d: %d\n', i, v)  # COMMENT THIS
    if (pg_type == "pg") {
      temp <- unlist(page_rank(g, directed=T)[1])#$vector
      PMF <- temp/sum(temp)
    }
    else if (pg_type == "median") {
      temp <- unlist(page_rank(g, directed=T)[1])
      node_ids <- c(1:900)
      temp2 <- cbind(node_ids, temp)
      rank_order <- temp2[order(temp2[,2], decreasing=F),]
      med <- vcount(g)%/%2
      temp2[,2] <- 0
      temp2[med:(med+1),2] <- 0.5
      PMF <- temp2[,2]
    }
    else if (pg_type=="equal") {
      temp <- unlist(page_rank(g, directed=T)[1])
      node_ids <- c(1:900)
      temp2 <- cbind(node_ids, temp)
      temp2[,2] <- 1/900
      med <- vcount(g)%/%2
      temp2[med:(med+1), 2] <- 0.5
      PMF <- temp2[,2]
    }
    else {
      PMF = transition_matrix[v, ]
    }
    if (runif(1) > t_prob) {
        v=sample(1:vcount(g), 1, prob = transition_matrix[v,])

    }
    else {
      if (pg_type == "1/N") {
        v = sample(1:vcount(g), 1)
      }
      else {
        v = sample(1:vcount(g), 1, prob = PMF)
      }
    }
    end_nodes[[i]] <- v
  }
  
  return(end_nodes)
}

random_walker_teleport = function (g, t, runs, t_prob, pg_type="1/N") {
  dist_list <- matrix(NA, nrow=runs, ncol=t) #vector(mode = 'list', length = runs) #
  deg_list <- array(1:runs)
  node_list <- vector(mode='list', length = t)
  temp <- vector(mode='list', length=t)
  for (i in 1:runs) {
    start_node = sample(1:vcount(g), 1)
    node_list <- cust_rand_walk(g, num_steps = t, start=start_node, t_prob=t_prob, pg_type=pg_type)
    j <- 1
    assign("temp", NULL)
    for (node in node_list){
      temp[j] <- distances(g, start_node, to=node)
      j <- j + 1
    }
    dist_list[i,] <- temp
    last_node <- node_list[t]
    deg_list[i] <- last_node
  }
  return(deg_list)
}

end_nodes_d <- random_walker_teleport(g=p, t=100, runs=1000, t_prob=0.2)

db_prob <- node_prob(end_nodes_d)
plot(x=1:900, y=db_prob, main= "Probability of Walker Visiting a Node, Teleportation", xlab="Node Index", ylab="Probability")
plot(x=unlist(cust_degree(end_nodes_d)), y = db_prob, main = "Node Degree vs. Visit Probability", xlab = "Degree", ylab="Visit Probability")

end_nodes_d5 <- random_walker_teleport(g=p, t=100, runs=1000, t_prob=0.5)

db5_prob <- node_prob(end_nodes_d5)
plot(x=1:900, y=db5_prob, main= "Probability of Walker Visiting a Node, Teleportation Probability 0.5", xlab="Node Index", ylab="Probability")
plot(x=unlist(cust_degree(end_nodes_d5)), y = db5_prob, main = "Node Degree vs. Visit Probability", xlab = "Degree", ylab="Visit Probability")

end_nodes_d8 <- random_walker_teleport(g=p, t=100, runs=1000, t_prob=0.8)

db8_prob <- node_prob(end_nodes_d8)
plot(x=1:900, y=db8_prob, main= "Probability of Walker Visiting a Node, Teleportation Probability 0.8", xlab="Node Index", ylab="Probability")
plot(x=unlist(cust_degree(end_nodes_d8)), y = db8_prob, main = "Node Degree vs. Visit Probability", xlab = "Degree", ylab="Visit Probability")

#########################################################################

#4a
end_nodes_4 <- random_walker_teleport(p, t=100, runs=1000, t_prob = 0.2, pg_type="pg")

d4b_prob <- node_prob(end_nodes_4)
plot(x=1:900, y=d4b_prob, main="Probability of Walker Visiting a Node, Normal PageRank", xlab="Node Index", ylab="Probability")
plot(x=unlist(cust_degree(end_nodes_4)), y = d4b_prob, main = "Node Degree vs. Visit Probability", xlab = "Degree", ylab="Visit Probability")


#4b
b4 <- random_walker_teleport(p, t=100, runs=1000, t_prob = 0.2, pg_type="median")

b4_prob <- node_prob(b4)
plot(x=1:900, y=b4_prob, main="Probability of Walker Visiting a Node, Median PageRank", xlab="Node Index", ylab="Probability")
plot(x=unlist(cust_degree(b4)), y = b4_prob, main = "Node Degree vs. Visit Probability", xlab = "Degree", ylab="Visit Probability")



#4c

c4 <- random_walker_teleport(p, t=100, runs=1000, t_prob = 0.2, pg_type ="equal")

c4_prob <- node_prob(c4)
plot(x=1:900, y=c4_prob, main="Probability of Walker Visiting a Node, Personalized PageRank", xlab="Node Index", ylab="Probability")
plot(x=unlist(cust_degree(c4)), y = c4_prob, main = "Node Degree vs. Visit Probability", xlab = "Degree", ylab="Visit Probability")


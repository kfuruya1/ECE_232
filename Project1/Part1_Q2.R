library(igraph)

#2a
g_2a <- sample_pa(n=1050, m=1, directed=F)

for (i in 1:100){
  g <- sample_pa(n=1050, m=1, directed=FALSE)
  if (!is_connected(g)) {
    print("Not connected")
  }
}

#2b
mod_var = function(g) {
  fastgreedycom <- cluster_fast_greedy(g)
  print(modularity(fastgreedycom))
  print(assortativity_degree(g, directed=FALSE))
}

mod_var(g_2a)

#2c
g10500 <- sample_pa(n=10500, m=1, directed=FALSE)

mod_var(g10500)

#2d

plot_deg_dist = function(deg_dist) {
  clean <- which(deg_dist != 0, arr.ind=TRUE)
  ax <- log(seq(1:length(deg_dist)))[clean]
  ay <- log(deg_dist)[clean]
  lr <- lm(ay~ax)
  plot(ax, ay, abline(lr, col="red"), main="Degree distribution in log-log scale", xlab="log(Degrees)", ylab="log(Distribution)")
  cf <- coef(lr)
  print(cf[2])
}

plot_deg_dist(degree_distribution(g_2a))

plot_deg_dist(degree_distribution(g10500))


#2e

rand_deg_dist = function(g) {
  a_degrees <- vector(mode = 'list', length=1000)
  a_size = length(degree_distribution(g))

  a_count <- numeric(a_size)
  a_sum <- 0
  
  for (i in 1:1000){
    
    ia <- sample(V(g), 1)
    a_adj <- adjacent_vertices(g, ia, mode="all")
    un_a_adj <- unlist(a_adj)
    ja <- sample(un_a_adj, 1)
    a_degrees[[i]] <- length(adjacent_vertices(g, ja, mode="all")[[1]])
    a_sum <- a_sum + a_degrees[[i]]
    a_count[[a_degrees[[i]]]] <- a_count[[a_degrees[[i]]]] + 1
  }
  a_dist <- numeric(a_size)
  for (i in 1:a_size){
    a_dist[[i]] <- a_count[[i]]/a_sum
  }
  plot_deg_dist(a_dist)
}

rand_deg_dist(g_2a)
rand_deg_dist(g10500)


#2f
degree_age = function(n, m) {
ba_degrees <- vector(mode='list', length = n)

for (i in 1:n) {
  ba <- sample_pa(n=n, m=m, directed=FALSE)
  ba_degrees[[i]] <- degree(ba)
}

ba_avg_degree <- vector(mode='list', length = n)
for (i in 1:n) {
  j_sum <- 0
  for (j in ba_degrees) {
    j_sum <- j_sum + j[i]
  }
  ba_avg_degree[[i]] <- j_sum/n
}

plot(x=c(n:1), y = ba_avg_degree, main="Expected Degree Based on Node Age", xlab = "Age of node", ylab="Degree")
}

degree_age(1050, 1)
#2g

#a
for (i in 1:100){
  g2 <- sample_pa(n=1050, m=2, directed=FALSE)
  g6 <- sample_pa(n=1050, m=6, directed=FALSE)
  if (!is_connected(g2)) {
    print("m = 2 Not connected")
  }
  if (!is_connected(g6)) {
    print("m = 6 not connected")
  }
}

#b
mod_var(g2)
mod_var(g6)

#c
g10500_2 <- sample_pa(n=10500, m=2, directed=FALSE)
g10500_6 <- sample_pa(n=10500, m=6, directed=FALSE)

mod_var(g10500_2)
mod_var(g10500_6)

plot(g_2a, main="n=1050, m=1")
plot(g10500, main="n=10500, m=1")

plot(g2, main="n=1050, m=2")
plot(g6, main="n=1050, m=6")

plot(g10500_2, main="n=10500, m=2")
plot(g10500_6, main="n=10500, m=6")



#d
plot_deg_dist(degree_distribution(g2))
plot_deg_dist(degree_distribution(g6))

plot_deg_dist(degree_distribution(g10500_2))
plot_deg_dist(degree_distribution(g10500_6))

#e
rand_deg_dist(g2)
rand_deg_dist(g6)

rand_deg_dist(g10500_2)
rand_deg_dist(g10500_6)

#f
degree_age(1050, 2)
degree_age(1050, 6)

degree_age(10500, 2)
degree_age(10500, 6)

#2h

g_h <- sample_pa(n=1050, m=1, directed=F)
g_new <- sample_degseq(degree(g_h), method="simple.no.multiple")
g_new2 <- sample_degseq(degree(g_h), method="vl")

fast_h <- cluster_fast_greedy(g_h)
fast_new <- cluster_fast_greedy(g_new)
fast_new2 <- cluster_fast_greedy(g_new2)

modularity(fast_h)
modularity(fast_new)
modularity(fast_new2)


plot(fast_h, g_h)
plot(fast_new, g_new)
plot(fast_new2, g_new2)




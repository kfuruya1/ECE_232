library(igraph)


#3a

g_3a <- sample_pa_age(n=1050, pa.exp = 1, aging.exp = -1, directed=F)
deg_dist_3a <- degree_distribution(g_3a)
plot(g_3a, main = "Network with New Preferential Attachment Method")
plot(deg_dist_3a, main = "Degree Distribution", xlab = "Degrees", ylab="Distribution")


clean_g_3a <- which(deg_dist_3a != 0, arr.ind=TRUE)
x3a <- log(seq(1:length(deg_dist_3a)))[clean_g_3a]
y3a <- log(deg_dist_3a)[clean_g_3a]
lr_3a <- lm(y3a~x3a)
cf_3a <- coef(lr_3a)
print(cf_3a[2])
plot(x3a, y3a, abline(lr_3a, col="red"), main = "Log Scale Degree Distribution", xlab="log(Degree)", ylab="log(Distribution)")


#3b

fast_3a <- cluster_fast_greedy(g_3a)

modularity(fast_3a)
plot(fast_3a, g_3a, main="Community Structure of Age-Based Preferential Attachment")

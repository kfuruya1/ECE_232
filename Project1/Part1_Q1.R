library(igraph)
############
#1a



#create undirected random networks with n=900 nodes and p= 0.002, 0.006, 0.012, 0.045, and 0.1

g002 <- sample_gnp(n=900, p=0.002, directed = FALSE)
g006 <- sample_gnp(n=900, p=0.006, directed = FALSE)
g012 <- sample_gnp(n=900, p=0.012, directed = FALSE)
g045 <- sample_gnp(n=900, p=0.045, directed = FALSE)
g1 <- sample_gnp(n=900, p=0.1, directed = FALSE)

#plot the degree distribution

plot(degree_distribution(g002), main="Degree Distribution for p=0.002", xlab = "Degrees", ylab = "Distribution")
plot(degree_distribution(g006), main="Degree Distribution for p=0.006", xlab = "Degrees", ylab = "Distribution")
plot(degree_distribution(g012), main="Degree Distribution for p=0.012", xlab = "Degrees", ylab = "Distribution")
plot(degree_distribution(g045), main="Degree Distribution for p=0.045", xlab = "Degrees", ylab = "Distribution")
plot(degree_distribution(g1), main="Degree Distribution for p=0.1", xlab = "Degrees", ylab = "Distribution")


#calculate mean and variance and compare to theoretical values

m002 <- mean(degree(g002))
m006 <- mean(degree(g006))
m012 <- mean(degree(g012))
m045 <- mean(degree(g045))
m1 <- mean(degree(g1))

v002 <- var(degree(g002))
v006 <- var(degree(g006)) 
v012 <- var(degree(g012))
v045 <- var(degree(g045))
v1 <- var(degree(g1))

tm002 <- 899*0.002
tm006 <- 899*0.006
tm012 <- 899*0.012
tm045 <- 899*0.045
tm1 <- 899*0.1

tv002 <- tm002*0.998
tv006 <- tm006*0.994
tv012 <- tm012*0.988
tv045 <- tm045*0.955
tv1 <- tm1*0.9

df <- data.frame(p = c('0.002', '0.006', '0.012', '0.045', '0.1'),
                 Mean = c(m002, m006, m012, m045, m1),
                 Theoretical_Mean = c(tm002, tm006, tm012, tm045, tm1),
                 Variance = c(v002, v006, v012, v045, v1),
                 Theoretical_Variance = c(tv002, tv006, tv012, tv045, tv1))

print(df)

##############
#1b

#are all random instances connected?

check_connected = function (prob) {
  for (i in 1:100) {
    if (!is_connected(sample_gnp(n=900, p=prob, directed=F))) {
      print("Not connected")
      break
    }
  }
}

check_connected(0.002)
check_connected(0.006)
check_connected(0.012)
check_connected(0.045)
check_connected(0.1)

#calculate probability that it is connected

probs = c(0.002, 0.006, 0.012, 0.045, 0.1)
for (p in probs){
  x <- 0
  for (i in 1:100) {
    g <- sample_gnp(n=900, p=p, directed = FALSE)
    if (is_connected(g)){
      x <- x+1
    }
  }
  print(x/100)
}

#find one instance of GCC and its diameter
for (p in probs){
  for (i in 1:100) {
    g <- sample_gnp(n=900, p=p, directed = FALSE)
    if (!is_connected(g)){
      print(diameter(largest_component(g)))
      break
    }
  }
  print(diameter(largest_component(g)))
}

#1ci
#x_val <- c(0, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0055, 0.006, 0.0065, 0.007, 0.0075, 0.008, 0.0085, 0.009, 0.0095, 0.01, 0.0105, 0.011, 0.0115, 0.012, 0.0125, 0.013, 0.0135, 0.014, 0.0145, 0.015)
x_val <- c(seq(from=0, to=0.015, by=0.0005))
a <- 0
y <- 0
gccx <- vector(mode = 'list', length = 3000)
gccy <- vector(mode = 'list', length = 3000)
gcc_avg <- vector(mode = 'list', length = 31)
for (p in seq(from = 0, to = 0.015, by = 0.0005)) {
  a <- a + 1
  x <- 0
  avg <- 0
  for (i in 1:100){
    y <- y + 1
    g <- sample_gnp(900, p, directed = FALSE)
    if (is_connected(g)){
      x <- x + 1
      z <- 1
    }
    else{
      z <- (components(largest_component(g))$csize)/900
    }
    gccx[[y]] <- p
    gccy[[y]] <- z
    avg <- avg + z
  }
  gcc_avg[[a]] <- avg/100
}

plot(x=gccx, y=gccy, main = "Normalized and Average GCC Size", xlab = "Probability (p)", ylab = "Normalized GCC Size", xaxt="n")
lines(x_val, gcc_avg, col = "red")
points(x_val, gcc_avg, col = "red", pch = 1)
axis(1, at = seq(0, 0.015, by = 0.0005), las = 2)
abline(v=(1/900), col="blue")
abline(v=(log(900)/900), col="green")
legend(x="bottomright", legend=c("Norm GCC size", "Avg GCC size", "O(1/N)", "O(ln(n)/n)"), col = c("black", "red", "blue", "green"), lwd = 1, lty = c(NA, 1, 1, 1), pch=1)

#1cii
plot(x=gccx, y=gccy, ylim=c(0.98, 1), xaxt = "n", yaxt="n", main = "Normalized and Average GCC Size Above 0.99", xlab = "Probability (p)", ylab = "Normalized GCC Size")
axis(1, at = seq(0, 0.015, by = 0.001), las = 2)
axis(2, at = seq(0.98, 1, by =0.01), las = 2)
points(x_val, gcc_avg, col="red")
abline(v=(1/900), col="blue")
abline(v=(log(900)/900), col="green")
legend(x="bottomright", legend=c("Norm GCC size", "Avg GCC size", "O(1/N)", "O(ln(n)/n)"), col = c("black", "red", "blue", "green"), lwd = 1, lty = c(NA, 1, 1, 1), pch=1)


#1di
n_array <- vector(mode='list', length = 100)
p1 <- 0
j <- 0
size_array <- vector(mode = 'list', length = 100)
for (n in seq(100, 10000, by=100)) {
  j <- j + 1
  p1 <- 0.5/n
  n_array[[j]] <- n
  gcc_sum <- 0
  for (i in 1:100){
    g <- sample_gnp(n, p1)
    gcc_sum <- gcc_sum + components(largest_component(g))$csize
  }
  size_array[[j]] <- gcc_sum/100
}

plot(x = n_array, y = size_array, main = "Expected GCC Size", xlab = "n", ylab = "Size")

#1dii
n_arrayii <- vector(mode='list', length = 100)
p1 <- 0
j <- 0
size_arrayii <- vector(mode = 'list', length = 100)
for (n in seq(100, 10000, by=100)) {
  j <- j + 1
  p1 <- 1/n
  n_arrayii[[j]] <- n
  gcc_sum <- 0
  for (i in 1:100){
    g <- sample_gnp(n, p1)
    gcc_sum <- gcc_sum + components(largest_component(g))$csize
  }
  size_arrayii[[j]] <- gcc_sum/100
}
plot(x = n_arrayii, y = size_arrayii, main = "Expected GCC Size", xlab = "n", ylab = "Size")


#1diii
n_array1 <- vector(mode='list', length = 100)
p1 <- 0
j <- 0
size_array1 <- vector(mode = 'list', length = 100)
size_array2 <- vector(mode = 'list', length = 100)
size_array3 <- vector(mode = 'list', length = 100)
for (n in seq(100,10000, by = 100)) {
  j <- j + 1
  p1 <- 1.15/n
  p2 <- 1.25/n
  p3 <- 1.35/n
  n_array1[[j]] <- n
  gcc_sum1 <- 0
  gcc_sum2 <- 0
  gcc_sum3 <- 0
  for (i in 1:100){
    g1 <- sample_gnp(n, p1)
    gcc_sum1 <- gcc_sum1 + components(largest_component(g1))$csize
    g2 <- sample_gnp(n, p2)
    gcc_sum2 <- gcc_sum2 + components(largest_component(g2))$csize
    g3 <- sample_gnp(n, p3)
    gcc_sum3 <- gcc_sum3 + components(largest_component(g3))$csize
  }
  size_array1[[j]] <- gcc_sum1/100
  size_array2[[j]] <- gcc_sum2/100
  size_array3[[j]] <- gcc_sum3/100
}

plot(x = n_array1, y = size_array1, main = "Expected GCC Size", xlab = "n", ylab = "Size")
points(x=n_array1, y=size_array2, col="red", pch = 2)
points(x=n_array1, y=size_array3, col="blue", pch = 5)
legend(x="bottomright", legend = c("c=1.15", "c=1.25", "c=1.35"), col=c("black", "red", "blue"), pch = c(1, 2, 5))

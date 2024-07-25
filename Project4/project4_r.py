# -*- coding: utf-8 -*-
"""Project4_R.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1YA4y8eLr8CUkzNwpI9_0KenZUxfKQ6-5
"""

if (!require("igraph")) install.packages("igraph")
library ("igraph")
if (!require("clevr")) install.packages("clevr")
library ("clevr")

unzip('finance_data.zip')

path = "finance_data/data/"
out.file<-""
file.names <- dir(path, pattern =".csv")
m <- matrix(, nrow = 0, ncol = 765)
sectors.table <- read.table("finance_data/Name_sector.csv",header=TRUE, sep=",", stringsAsFactors=TRUE)
sectors = c()
sectors.names <- c()
for(i in 1:length(file.names)){
    file <- read.table(paste("finance_data/data/",file.names[i],sep = ""),header=TRUE, sep=",", stringsAsFactors=FALSE)
    if(length(file$Close)==765){
        m <- rbind(m, matrix(file$Close, nrow=1, ncol=765))
        mystr <- substr(file.names[i], 1, nchar(file.names[i])-4)
        sector <- sectors.table$Sector[which(sectors.table$Symbol == mystr)]
#         print(sector)
#        sectors = c(sectors, as.factor(sector))
        sectors.names = c(sectors.names, toString(sector))
    }
}
com_num = length(sectors.names)
data_num = ncol(m)
sector.set<-as.factor(sectors.names)
sectors.index<-as.numeric(sector.set)
num_sector <- length(unique(sectors.index))
table(sector.set)



qit <- function(p, t){
  q <- (p[t] - p[t-1])/p[t-1]
  return(q)
}

rit <- function(p, t){
  r <- log(1+qit(p, t))
  return(r)
}

pij <- function(i, j){
  ri <- 0
  rj <- 0
  ri_sum <- 0
  rj_sum <- 0
  rij <- 0
  ri2 <- 0
  rj2 <- 0

  for (t in 2:2765){
    ri <- rit(i, t)
    rj <- rit(j, t)
    ri_sum <- ri_sum + ri
    rj_sum <- rj_sum + rj
    rij <- rij + (ri*rj)
    ri2 <- ri2 + ri^2
    rj2 <- rj2 + rj^2
  }
  a <- rij/765
  b <- ri_sum/765
  c <- rj_sum/765
  d <- ri2/765
  e <- rj2/765

  p <- (a - (b*c))/sqrt((d-(a^2))*(e-(b^2)))
  return(p)
}

wij <- function(i, j){
  w <- sqrt(2*(1-pij(i, j)))
  return(w)
}

#correlation
cal_cor<-function(mtrx){
  edge_weights <- list()
  edge_names <- list()
  graph_info <- list()
  lngth <- nrow(mtrx)
  for (rw in 1:(lngth-1)){
    for (rw2 in rw+1:lngth){
      edg <- wij(mtrx[rw], mtrx[rw2])
      edge_weights <- append(edge_weights, edg)
      tple <- c(rw, rw2)
      edge_names <- append(edge_names, tple)
      n_e <- c(rw, rw2, edg)
      graph_info <- append(graph_info, n_e)
    }
  }
  val <- c(edge_weights, edge_names, graph_info)
  return(val)
  }

graph_list <- cal_cor(m)
wghts <- graph_list[1]
nms <- graph_list[2]
g_info <- graph_list[3]

hist(unlist(wghts))
# ADMM for Charging, INFORMs 2018
library(ggplot2)
library(reshape2)
library(data.table)

# clear all list----
rm(list = ls())

plot_iter_proposal <- function(data, yLab){
  quartz(width = 8, height = 5)
  
  ggplot(data=data, aes(x=iteration, y=value, color = variable)) + 
    geom_line(size=0.5, alpha=1) + 
    #theme(legend.position = c(0.8,0.2)) +
    theme(legend.position = "none", axis.title=element_text(size=14,face="bold") , axis.text=element_text(size=14)) +
    coord_cartesian(xlim = c(0, max(data$iteration))) + scale_x_continuous(breaks=seq(0, max(data$iteration), 20), labels = seq(0, max(data$iteration), 20)) +
    #coord_cartesian(ylim = c(-20, 30))  +
    labs(x="Iteration", y=yLab, color = "Scenario & Charging Places")
}

plot_iter <- function(data, yLab){
  quartz(width = 8, height = 5)
  
  ggplot(data=data, aes(x=iteration, y=value, color = variable)) + 
    geom_line(size=0.5, alpha=1) + 
    #theme(legend.position = c(0.8,0.2)) +
    theme(legend.position = 'right') +
    coord_cartesian(xlim = c(0, max(data$iteration))) + scale_x_continuous(breaks=seq(0, max(data$iteration), 20), labels = seq(0, max(data$iteration), 20)) +
    #coord_cartesian(ylim = c(-20, 30))  +
    labs(x="Iteration", y=yLab, color = "Scenario & Charging Places")
  
}

# set working directory ----
setwd("./Results/")
#setwd("./Results_20210510_0033/") #mac

# read in data ----
es <- read.csv("Resulting_exsu.csv", header = FALSE, sep = ",")
price <- read.csv("Resulting_prices.csv", header = FALSE, sep = ",")
sd <- read.csv("Resulting_scendiff.csv", header = FALSE, sep = ",")

# propose data ----
# transpose, so that each column is one line in plot
t_es <- transpose(es)
t_price <- transpose(price)
t_sd <- transpose(sd)
# assign rownames to colnames
temp_1 <- paste(t_es[1,],t_es[2,],sep=",")
temp_2 <- paste(t_price[1,],t_price[2,],sep=",")
temp_3 <- paste(t_sd[1,],t_sd[2,],sep=",")
t_es <- data.frame(apply(t_es[-c(1:2),],2,as.numeric))
t_price<-data.frame(apply(t_price[-c(1:2),], 2,as.numeric))
t_sd<-data.frame(apply(t_sd[-c(1:2),], 2,as.numeric))

colnames(t_es) <- temp_1
colnames(t_price) <- temp_2
colnames(t_sd) <- temp_3
# only keep the rows with non NA
t_es <- t_es[complete.cases(t_es), ]
t_price <- t_price[complete.cases(t_price), ]
t_sd <- t_sd[complete.cases(t_sd), ]
# add the iteration number for x axis
t_es$iteration <- 1:nrow(t_es)
t_price$iteration <- 1:nrow(t_price)
t_sd$iteration <- 1:nrow(t_sd)

# plot ----
data <- melt(t_es, id=c("iteration"))
yLab <- "Excess Supply (kWh)"
plot_iter_proposal(data, yLab)

data <- melt(t_price, id=c("iteration"))
yLab <- "Locational Charging Prices ($/kWh)"
plot_iter_proposal(data, yLab)

data <- melt(t_sd, id=c("iteration"))
yLab <- "Investment Deviation from Scenario Average (kW)"
plot_iter_proposal(data, yLab)


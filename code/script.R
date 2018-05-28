#install.packages(c("ape","mclust"))
library(ape)
library(reshape2)

#################################################################
#  Modify this section to load in different examples. To do so, #
#  put a comment (#) in front of the lines that load in the     #
#  first example (data, tree, out_file). Remove the comments    #
#  next to the lines for the second or third example.           #
#################################################################

#read in data
data <- read.csv("../data/rWGD_8ssp_angiosperms_9655families.csv", row.names = 1) #Example 1
#data <- read.csv("../data/aWGD_8ssp_angiosperms_9603families.csv", row.names = 1) #Example 2
#data <- read.csv("../data/10spp_fungi_12107.csv", row.names = 1) #Example 3

#read in tree
tree <- read.tree(file = "../data/rWGD_8ssp_angiosperms_9655families_tree.txt") #Example 1
#tree <- read.tree(file = "../data/aWGD_8ssp_angiosperms_9603families_tree.txt") #Example 2
#tree <- read.tree(file = "../data/10species_fungi_tree.txt") #Example 3

#specify out file prefix
out_file <- "../results/rWGD_8ssp_angiosperms_9655families_0426_" #Example 1
#out_file <- "../results/aWGD_8ssp_angiosperms_9603families_0426_" #Example 2
#out_file <- "../results/10spp_fungi_12107_" #Example 3


#load load load, if you already ran the loop, then load in the .Rdata file

#Load for Example 1
load("../results/rWGD_8ssp_angiosperms_9655families_0426_dr.Rdata")
load("../results/rWGD_8ssp_angiosperms_9655families_0426_dr_p.Rdata")

#Load for Example 2
#load("../results/aWGD_8ssp_angiosperms_9603families_0426_dr.Rdata")
#load("../results/aWGD_8ssp_angiosperms_9603families_0426_dr_p.Rdata")

#Load for Example 3
#load("../results/10spp_fungi_12107_dr.Rdata")
#load("../results/10spp_fungi_12107_dr_p.Rdata")

#permutation number for p-value calculation
permut <- 100

########################################################
######      No need to modify this section     #########
########################################################

#calculate the length of the columns
len <- length(data[1,])

#calculate standard deviation of each gene family
data_sd <- as.numeric(apply(t(data[,2:len]), 2, sd))

#keep only data that has a sd of gene family size < 3 
data <- data[(data_sd < 3 & data_sd != 0),]

#do a loop and run through all pairwise comparisons

##the commented code below should only be run one time, which saves out the big data matrices
##double comments (##) are just regular comments. Single comments should be removed to run code (#)

##make a variable to hold the correlation coefficients
#dr <- matrix(NA, nrow = length(data[,1]), ncol = length(data[,1])) #distribution of r values 
#dr_p <- matrix(NA, nrow = length(data[,1]), ncol = length(data[,1])) #distribution of r values 

##time the loop
#start.time <- Sys.time()

##start the loop

# for (i in 1:(length(data[,1])-1)){ #length(data[1,])
#   for (j in (i+1):length(data[,1])){ #
# 
#     #calculate independent contrasts
#     x <- pic(data[i,1:len], tree)
#     y <- pic(data[j,1:len], tree)
# 
#     #get the correlation coefficient
#     dr[i,j] <- cor(x,y, method = "spearman")
#     temp <- array(NA, dim = permut)
#     for(k in 1:permut){
#       temp[k] <- cor(x,sample(y), method = "spearman")
#     }
#     dr_p[i,j] <- sum(temp >= dr[i,j])/permut
#   }
# }
# time.taken <- Sys.time() - start.time
# time.taken

##save the matrix, so you don't have to run it in the future and can just load the r matrix
#save(dr, file = paste("../results/", out_file, "dr.Rdata", sep = ""))
#save(dr_p, file = paste("../results/", out_file, "dr_p.Rdata", sep = ""))

######################
######################
######################

##code to run the analysis without phylogenetic independent contrasts

# #do a loop and run through all pairwise comparisons without pic
# dr_no_phylo <- matrix(NA, nrow = length(data[,1]), ncol = length(data[,1])) #distribution of r values 
# 
# start.time <- Sys.time()
# for (i in 1:(length(data[,1])-1)){ #length(data[1,])
#   for (j in (i+1):length(data[,1])){ #
#     
#     #calculate independent contrasts
#     x <- data[i,1:len]
#     y <- data[j,1:len]
#     
#     #get the correlation coefficient
#     dr_no_phylo[i,j] <- cor(t(x),t(y))
#     
#   }
# }
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# save(dr_no_phylo, file = paste("../results/", out_file, "dr_no_phylo.Rdata", sep = ""))
# 
## do a loop and run through all pairwise comparisons with random species
# random_compare <- matrix(NA, nrow = length(data[,1]), ncol = length(data[,1])) #distribution of r values 
# 
# start.time <- Sys.time()
#   for (i in 1:(length(data[,1])-1)){ #length(data[1,])
#     for (j in (i+1):length(data[,1])){ #
#     
#       #calculate independent contrasts
#       x <- sample(data[i,1:len])
#       y <- data[j,1:len]
#     
#       #get the correlation coefficient
#       random_compare[i,j] <- cor(t(x),t(y))
#     }
#   }
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# save(random_compare, file = paste("../results/", out_file, "random_compare.Rdata", sep = ""))

###########################################################

pdf(file = paste("../results/", out_file, "histograms.pdf", sep = ""))
temp <- hist(dr, col = "gray")
hist(dr[which(dr_p < 0.05, arr.ind = T)], col = "red", add = T, breaks = temp$breaks)
hist(dr[which(dr_p > 0.95, arr.ind = T)], col = "red", add = T, breaks = temp$breaks)
#hist(dr, freq = F, ylim = c(0, 3.1), col = "gray", main = out_file, xlab = "Correlation Coefficient")
#hist(dr_no_phylo, freq = F, ylim = c(0, 3.1), add = T, border = "blue")
#hist(random_compare, freq = F, ylim = c(0, 3.1), add = T, border = "red")
#legend("topleft",legend = c("PIC","No PIC","Random"), col = c("black", "blue", "red"), lwd = 1, bty = "n")
dev.off()

# pdf(file = paste("../results/", out_file, "densities.pdf", sep = ""))
# plot(density(dr[!is.na(dr)],bw = 0.07, na.rm = T), ylim = c(0, 3.1), col = "black", lwd = 2, main = out_file, xlab = "Correlation Coefficient")
# lines(density(dr_no_phylo[!is.na(dr_no_phylo)], bw = 0.07), col = "blue", lwd = 2, lty = 2)
# lines(density(random_compare[!is.na(random_compare)], bw = 0.07), col = "red", lwd = 2, lty = 2)
# legend("topleft",legend = c("PIC","No PIC","Random"), col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = 2, bty = "n")
# dev.off()

###########################################################

#read out all pairs

colnames(dr) <- rownames(dr) <- rownames(data)
dr_melt <- setNames(melt(dr), c('fam1', 'fam2', 'r'))
dr_melt <- dr_melt[!is.na(dr_melt[,3]),]

colnames(dr_p) <- rownames(dr_p) <- rownames(data)
dr_p_melt <- setNames(melt(dr_p), c('fam1', 'fam2', 'p'))
dr_p_melt <- dr_p_melt[!is.na(dr_p_melt[,3]),]

pairs_r_p <- data.frame(dr_melt, "p" = dr_p_melt[,3])

p_temp <- dr_p_melt$p
p_temp[dr_p_melt$p > 0.5] <- (1 - dr_p_melt$p[dr_p_melt$p > 0.5])
p_temp_adj <- p.adjust(p_temp, method = "fdr")

sum(p_temp_adj < 0.05)
length(p_temp_adj)
sum(p_temp_adj < 0.05)/length(p_temp_adj) * 100

temp <- hist(dr, col = "gray")
hist(dr[which(dr_p < 0.05, arr.ind = T)], col = "red", add = T, breaks = temp$breaks)
hist(dr[which(dr_p > 0.95, arr.ind = T)], col = "red", add = T, breaks = temp$breaks)

#write all pairs
#write.table(pairs_r_p, paste("../results/", out_file, "all_pairs_r_p.csv", sep = ""), sep = ",")

#write sig pairs
pos <- pairs_r_p[,4] < 0.05
neg <- pairs_r_p[,4] > 0.95
sig_pairs_r_p <- pairs_r_p[pos|neg,]
#write.table(sig_pairs_r_p, paste("../results/", out_file, "SIG_pairs_r_p.csv", sep = ""), sep = ",")

#total number of gene families after big sd removed
length(data[,1])
#total number pairs
length(pairs_r_p[,1])
#total number of sig pairs
length(sig_pairs_r_p[,1])
#percent of significant pairs
(length(sig_pairs_r_p[,1]) / length(pairs_r_p[,1])) * 100


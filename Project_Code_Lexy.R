install.packages("epitools")
library(epitools)
install.packages("PropCIs")
library(PropCIs)

# Just replace the two sets of parentheses in the actual with code with each of these lines
# t1 : (2,48,1,49) | (2, 50, 1, 50, 0.95) #
# t2 : (4,46,2,48) | (4, 50, 2, 50, 0.95) #
# t3 : (6,44,3,47) | (6, 50, 3, 50, 0.95) #
# t4 : (8,42,4,46) | (8, 50, 4, 50, 0.95) #
# t5 : (10,40,5,45) | (10, 50, 5, 50, 0.95) #
# t6 : (12,38,6,44) | (12, 50, 6, 50, 0.95) #
# t7 : (14,36,7,43) | (14, 50, 7, 50, 0.95) #
# t8 : (16,34,8,42) | (16, 50, 8, 50, 0.95) #
# t9 : (18,32,9,41) | (18, 50, 9, 50, 0.95) #
# t10 : (20,30,10,40) | (20, 50, 10, 50, 0.95) #
# t11 : (22,28,11,39) | (22, 50, 11, 50, 0.95) #
# t12 : (24,26,12,38) | (24, 50, 12, 50, 0.95) #
# t13 : (26,24,13,37) | (26, 50, 13, 50, 0.95) #
# t14 : (28,22,14,36) | (28, 50, 14, 50, 0.95) #
# t15 : (30,20,15,35) | (30, 50, 15, 50, 0.95) #
# t16 : (32,18,16,34) | (32, 50, 16, 50, 0.95) #
# t17 : (34,16,17,33) | (34, 50, 17, 50, 0.95) #
# t18 : (36,14,18,32) | (36, 50, 18, 50, 0.95) #
# t19 : (38,12,19,31) | (38, 50, 19, 50, 0.95) #
# t20 : (40,10,20,30) | (40, 50, 20, 50, 0.95) #
# t21 : (42,8,21,29) | (42, 50, 21, 50, 0.95) #
# t22 : (44,6,22,28) | (44, 50, 22, 50, 0.95) #
# t23 : (46,4,23,27) | (46, 50, 23, 50, 0.95) #
# t24 : (48,2,24,26) | (48, 50, 24, 50, 0.95)


table_hw7 <- matrix(c(48,2,24,26),ncol=2,byrow=TRUE)
colnames(table_hw7) <- c("Yes","No")
rownames(table_hw7) <- c("A","B")
table_hw7 <- as.table(table_hw7)
table_hw7
  # wald
oddsratio.wald(table_hw7, y = NULL,
          conf.level = 0.95,
          rev = c("neither"),
          correction = FALSE,
          verbose = FALSE)
  # fisher
oddsratio.fisher(table_hw7, y = NULL,
          conf.level = 0.95,
          rev = c("neither"),
          correction = FALSE,
          verbose = FALSE)
  # score
orscoreci(48, 50, 24, 50, 0.95)

#Lexy: More on the confidence interval graphs: Posted 11/18/18
# m1 : (4,1,2,3) | (4, 5, 2, 5, 0.95) #
# m2 : (8,2,4,6) | (8, 10, 4, 10, 0.95) #
# m3 : (16,4,8,12) | (16, 20, 8, 20, 0.95) #
# m4 : (32,8,16,24) | (32, 40, 16, 40, 0.95) #
# m5 : (64,16,32,48) | (64, 80, 32, 80, 0.95) #
# m6 : (128,32,64,96) | (128, 160, 64, 160, 0.95) #
# m7 : (256,64,128,192) | (256, 320, 128, 320, 0.95) #
# m8 : (512,128,256,384) | (512, 640, 256, 640, 0.95) 

table_hw7 <- matrix(c(512,128,256,384),ncol=2,byrow=TRUE)
colnames(table_hw7) <- c("Yes","No")
rownames(table_hw7) <- c("A","B")
table_hw7 <- as.table(table_hw7)
table_hw7
# wald
oddsratio.wald(table_hw7, y = NULL,
               conf.level = 0.95,
               rev = c("neither"),
               correction = FALSE,
               verbose = FALSE)
# fisher
oddsratio.fisher(table_hw7, y = NULL,
                 conf.level = 0.95,
                 rev = c("neither"),
                 correction = FALSE,
                 verbose = FALSE)
# score
orscoreci(512, 640, 256, 640, 0.95)

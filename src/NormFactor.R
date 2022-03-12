### This script is to measure the normal factor using the ERCC
library(tidyverse)
library(ggpubr)
theme_set(theme_pubr())

#install.packages('reshape2')
library(reshape2)



setwd("~/Documents/Projects/RNA_DDR/Results/")
MT<-read.table("ERCC_ZC12.txt",header = T)

### SampleA SampleB
### 10  20
### 30  40
### 1 5
### 1000  3000
###
attach(MT)
## SampleA vs SampleB

help(geom_step)

library(ggplot2)

## Visualization, create a scatter plot displaying the scale uniqe and add a smooth line for 
ggplot(MT,aes(x=log10(SampleA+0.01),y=log10(SampleB+0.01)))+  
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

##  compute the correlation coefficient between the two variables using the R function cor():
cor(log10(SampleB+0.01),log10(SampleA+0.01))

## Computing the linear regression
model <- lm(log10(SampleB+0.01)~log10(SampleA+0.01), data = MT)
summary(model)

## ----echo=F,knitr-options,message=FALSE, warning=FALSE------------------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')
options(warn=-1)

## -----------------------------------------------------------------------------
library(ISnorm)
library(dbscan)

## -----------------------------------------------------------------------------
# check raw counts matrix
mat[1:3,1:3]
# In order to decrease run time, we randomly select 5000 genes and 300 cells
set.seed(1)
sub.mat=mat[sample(1:nrow(mat), size = 5000),sample(1:ncol(mat), size = 300)]
nrow(sub.mat)

## -----------------------------------------------------------------------------
gene_dis<-calculate.dis(sub.mat,detection_rate=0.90)

## -----------------------------------------------------------------------------
spike_candidate<-dbscan.pick(dis=gene_dis,ngene=(1:floor(nrow(gene_dis)/25))*5)

## -----------------------------------------------------------------------------
str(spike_candidate)

## -----------------------------------------------------------------------------
candidate_res<-candidate.norm(mat=sub.mat,spike_candidate=spike_candidate)

## -----------------------------------------------------------------------------
names(candidate_res)

## -----------------------------------------------------------------------------
candidate_res$sf[1:3,1:3]

## -----------------------------------------------------------------------------
candidate_res$inst[1:3,1:3]

## -----------------------------------------------------------------------------
# plot every candidate geneset instability score
candidate_res$inst[1:3,1:3]

## -----------------------------------------------------------------------------
ISnorm_res<-opt.candidate(mat=mat,candidate_res=candidate_res)

## -----------------------------------------------------------------------------
names(ISnorm_res)


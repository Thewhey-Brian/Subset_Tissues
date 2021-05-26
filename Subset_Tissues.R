# Code for selecting subset of tissues
# --------------------------------------------------------
# Load needed packages
# --------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(qqman)
library(ggrepel)
library(plotly)
library(manhattanly)
library(RColorBrewer)
library(reshape)
library(Matrix)
library(MASS)
library(susieR)

# --------------------------------------------------------
# Load data
# --------------------------------------------------------
load("./cov_matrix.RData")
load("./mean_matrix.RData")
load("./ad_TWAS.RData")

# --------------------------------------------------------
# Needed functions
# --------------------------------------------------------
# The GAUSS function
Cal_S_EqualWeight <- function(Z){
  K<-length(Z)
  Z1<-sort(Z, decreasing=TRUE)
  S<-rep(0,K)
  for(i in 1:K){
    S[i]<-sum(Z1[1:i])/sqrt(i)
  }
  id.max<-which(S == max(S,na.rm =T))
  re<-list(maxS = S[id.max], idx = names(Z1[1:id.max]), S=S)
  return(re)
}

#' Simulate samples from null 
#'
#' @param gene The gene name
#' @param cov_matrix The null covariance matrix for all genes
#' @param mean_matrix The null mean matrix for all genes
#' @param n Number of samples needed
#'
#' @return A matrix of n samples for all existed tissues
#' @export
sim_null <- function(gene, cov_matrix, mean_matrix, n) {
  cov_gene <- cov_matrix[[gene]]
  mea_gene <- mean_matrix[[gene]]
  return(mvrnorm(n, mea_gene, cov_gene))
}

#' Get GDP-estimated p-values
#'
#' @param Nexc Exceedances threshold
#' @param n Number of simulated samples
#' @param list_stats List of statistics from certain test
#' @param sub_out Result from GAUSS function
#'
#' @return An estimated p-value
#' @export
#'
estp_gdp <- function(Nexc, n, list_stats, sub_out) {
  y_star <- sort(list_stats, decreasing = T)
  t <- (y_star[Nexc] + y_star[Nexc + 1]) / 2
  z <- y_star[y_star - t > 0] - t
  x_bar <- mean(z)
  x_var <- var(z)
  alpha <- x_bar * (x_bar ^ 2 / x_var + 1) / 2
  k <- (x_bar ^ 2 / x_var - 1) / 2
  while (sub_out$maxS > alpha / k & k > 0) {
    z <- c(z, sub_out$maxS)
    x_bar <- mean(z)
    x_var <- var(z)
    alpha <- x_bar * (x_bar ^ 2 / x_var + 1) / 2
    k <- (x_bar ^ 2 / x_var - 1) / 2
    cat(sub_out$maxS, ">= a/k")
  }
  if (k == 0) {
    omFz <- exp(-(sub_out$maxS - t) / alpha)
  }
  else {
    omFz <- (1 - k * (sub_out$maxS - t) / alpha) ^ (1 / k)
  }
  p_value <- Nexc / n * omFz
  return(p_value)
}

select <- dplyr::select

# --------------------------------------------------------
# Run through all genes
# --------------------------------------------------------
data_ad <- dat_ad_n[[1]] %>% 
  filter(GENE %in% all_gene)
gene_list <- data_ad$GENE
res <- data.frame()
for (gene in gene_list) {
  twas_z <- data.frame(data_ad) %>% 
    filter(GENE == gene) %>% 
    select_if(~ !is.na(.)) %>% 
    select(-"GENE")
  sub_out <- Cal_S_EqualWeight(twas_z)
  # try 10000 simulations first
  n <- 10000
  tryCatch({
    x <- sim_null(gene, cov_matrix, mean_matrix, n)
  }, error = function(e){cat("ERROR :", gene, conditionMessage(e), "\n")})
  list_stats <- c()
  for (i in 1:n) {
    tem <- Cal_S_EqualWeight(x[i, ])
    list_stats <- c(list_stats, tem$maxS)
  }
  # do more simulation if no touched the tail
  M <- sum(list_stats >= sub_out$maxS)
  if (M <= 10) {
    n <- 1e+06
    tryCatch({
      x <- sim_null(gene, cov_matrix, mean_matrix, n)
    }, error = function(e){cat("ERROR :", gene, conditionMessage(e), "\n")})
    list_stats <- c()
    for (i in 1:n) {
      tem <- Cal_S_EqualWeight(x[i, ])
      list_stats <- c(list_stats, tem$maxS)
    }
  }
  # tail approximation of p-value
  M <- sum(list_stats >= sub_out$maxS)
  if (M >= 10) {
    p_value <- M/n
  }
  else {
    p_value <- estp_gdp(250, n, list_stats, sub_out)
  }
  twas_p <- pnorm(as.numeric(twas_z), lower.tail = F)
  names(twas_p) <- names(twas_z)
  sig_p <- twas_p[order(twas_p[twas_p < 2.5e-6])]
  tem_out <- cbind(Gene = gene, 
                   Subset_Tissue = paste(sub_out$idx, collapse = ", "), 
                   Number_of_Tissues = length(sub_out$idx), 
                   MAX_S = sub_out$maxS, 
                   P_value = p_value, 
                   Sig_Tissues = paste(names(sig_p), sig_p, collapse = "; "), 
                   Number_of_Sig_Tissue = length(sig_p), 
                   Extra_Subset = paste(setdiff(sub_out$idx, names(sig_p)), collapse = ", "), 
                   Extra_TWAS = paste(setdiff(names(sig_p), sub_out$idx), collapse = ", "))
  res <- rbind(res, tem_out)
}


# --------------------------------------------------------
# Run through all genes from SuSiE Approach
# --------------------------------------------------------

#' The help function to run SuSiE
#'
#' @param gene The name of gene
#' @param data The data matrix with z-values from TWAS
#' @param cov_matrix The null covariance matrix for all genes
#' @param l Number of causal tissues expected
#'
#' @return A data frame of SuSiE summary results
#' @export
run_susie <- function(gene, data, cov_matrix, l = 10){
  dat_gene <- data %>% 
    filter(GENE == gene) %>% 
    select_if(!is.na(.))
  dat_gene <- dat_gene[, -1]
  # check whether z value is finite
  dat_gene <- dat_gene[, is.finite(as.numeric(dat_gene))] 
  cov_gene <- cov_matrix[[gene]]
  # unify all tissue names in both dat_gene and cov_gene
  rownames(cov_gene) <- sapply(rownames(cov_gene), 
                               function(x) gsub("[^[:alnum:]\\_]", "_", x), USE.NAMES = F)
  colnames(cov_gene) <- sapply(colnames(cov_gene), 
                               function(x) gsub("[^[:alnum:]\\_]", "_", x), USE.NAMES = F)
  names(dat_gene) <- sapply(names(dat_gene), 
                            function(x) gsub("[^[:alnum:]\\_]", "_", x), USE.NAMES = F)
  if(!setequal(names(dat_gene), rownames(cov_gene))){
    intrst <- intersect(names(dat_gene), rownames(cov_gene))
    dat_gene = dat_gene %>% select(all_of(intrst))
    cov_gene = cov_gene[intrst, intrst]
  }
  dat_input <- as.numeric(dat_gene)
  fitted_rss <- susie_rss(dat_input, cov_gene, L = l)
  nms <- names(dat_gene)
  smy <- summary(fitted_rss)
  out <- smy$vars %>% 
    mutate(tissues = nms[smy$vars$variable], 
           GENE = gene) %>% 
    filter(cs > 0) %>% 
    select(-variable)
  return(out)
}

data_ad <- dat_ad_n[[1]] %>% 
  filter(GENE %in% all_gene)
gene_list <- data_ad$GENE
res_susie <- data.frame()
for (gene in gene_list[1:10]) {
  tryCatch({
    res <- run_susie(gene, data_ad, cov_matrix, 10)
    tem_out <- cbind(Gene = gene, 
                     SuSiE_Tissue = paste0(res$tissues, ": ", 
                                           res$variable_prob, " ", 
                                           "[", res$cs, "]", 
                                           collapse = "; "), 
                     Number_of_SuSiE_Tissues = length(res$tissues))
    res_susie <- rbind(res_susie, tem_out)
  }, error = function(e){cat("ERROR :", gene, conditionMessage(e), "\n")})
}







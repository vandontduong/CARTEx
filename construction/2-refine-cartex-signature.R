# refine cartex signature
# Vandon Duong

setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)

library(glmnet)
library(ggplot2)
library(ggplotify)
library(ggpubr)

####################################################################################################
############################################# Functions ############################################
####################################################################################################

Z=function(s){
  s=as.numeric(s)
  z=(s - mean(s))/sd(s)
  return (z)
}

cartex <- function(expression_df, metadata, weights){
  
  # Load files
  normgenes <- read.csv(expression_df, row.names=1, header=TRUE)
  meta <- read.csv(metadata, row.names=1, header = TRUE)
  cartex <- read.csv(weights, header = TRUE)
  
  cartex_200=cartex$PC1
  names(cartex_200) <- cartex$X
  common <- intersect(names(cartex_200), rownames(normgenes))
  expr=t(normgenes[match(common, rownames(normgenes)),])
  cartex_200=cartex_200[match(common, names(cartex_200))]
  weights = as.numeric(cartex_200)
  
  if(all(colnames(expr)==names(cartex_200)) & all(rownames(meta)==rownames(expr))){
    scores = apply(expr, 1, function(x) x %*% weights)
    scores = data.frame(cbind(scores, samples=rownames(expr), meta))
    rownames(scores) <- rownames(expr)
    scores$cartex_score <- Z(scores$scores)
    scores=scores[,-1]
    return(scores)
  }
  return(paste0("genes or samples don't match"))
}

generate_figs = function(figure_object, file_name, dimensions){
  if(missing(dimensions)){
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object)
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white")
  } else {
    ggsave(filename = gsub(" ", "", paste(file_name,".pdf")), plot = figure_object, width = dimensions[1], height = dimensions[2])
    ggsave(filename = gsub(" ", "", paste(file_name,".jpeg")), plot = figure_object, bg = "white",  width = dimensions[1], height = dimensions[2])
  }
  return (paste("generating figure for ", file_name))
}

integerize = function(score){
  score_mod = round(score)
  score_mod[score_mod < -4] <- -5
  score_mod[score_mod > 4] <- 5
  return (score_mod)
}

prepare_weights = function(cvfit, lambda){
  weights <- as.matrix(coef(cvfit, s = lambda))
  weights <- weights[apply(weights!=0, 1, all),]
  weights <- data.frame(weights[-1]) # remove intercept
  colnames(weights) <- "PC1"
  return(weights)
}

####################################################################################################
######################################## Regress across alpha ######################################
####################################################################################################

out <- cartex("./data/2021-04-28-cartex-data-logTPM.csv","./data/2021-04-28-metadata.csv","../weights/cartex-200-weights.csv")
express_mat <- read.csv("./data/2021-04-28-cartex-data-logTPM.csv", row.names=1, header=TRUE)
CARTEx_weights <- read.csv("../weights/cartex-200-weights.csv", row.names=1, header = TRUE)
CARTEx_sig <- rownames(CARTEx_weights)

express_mat_sub <- t(express_mat[rownames(express_mat) %in% (CARTEx_sig),])

alpha_sequence <- seq(0, 1, length.out = 9)
CARTEx_size <- rep(0,9)
i = 0
CARTEx_genesets <- list()

for (a in alpha_sequence){
  cvfit_glm <- cv.glmnet(express_mat_sub, as.matrix(out$cartex_score), type.measure = "mse", nfolds = 5, alpha = a)
  cv_tmp <- data.frame(cvfit_glm$lambda, cvfit_glm$cvm, cvfit_glm$cvsd, cvfit_glm$cvup, cvfit_glm$cvlo)
  colnames(cv_tmp) <- c("lambda", "mse", "sde", "upper", "lower")
  weights_reg <- prepare_weights(cvfit_glm, "lambda.1se")
  i <- i + 1
  CARTEx_size[i] <- dim(weights_reg)[1]
  key_name <- paste("alpha",as.character(a), sep = "_")
  CARTEx_genesets[[key_name]] <- rownames(weights_reg)
}

explore_space <- data.frame(alpha_sequence, CARTEx_size)

plt_cvglmnet_alpha_size <- ggplot(explore_space, aes(x=alpha_sequence, y=CARTEx_size)) + 
  geom_point(size=2) + geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95)

generate_figs(plt_cvglmnet_alpha_size, "./plots/plt_cvglmnet_alpha_size")

refined_set <- Reduce(union, list(CARTEx_genesets$alpha_1, CARTEx_genesets$alpha_0.825, CARTEx_genesets$alpha_0.75, 
                                  CARTEx_genesets$alpha_0.625, CARTEx_genesets$alpha_0.5, 
                                  CARTEx_genesets$alpha_0.375, CARTEx_genesets$alpha_0.25))

length(refined_set)
CARTEx_refined <- CARTEx_weights[rownames(CARTEx_weights) %in% refined_set,]
names(CARTEx_refined) <- refined_set
write.csv(CARTEx_refined, "./data/cartex-84-weights.csv")


####################################################################################################
####################################### Regress specific alpha #####################################
####################################################################################################


# Lasso regression

cvfit_lasso <- cv.glmnet(express_mat_sub, as.matrix(out$cartex_score), type.measure = "mse", nfolds = 5, alpha = 1)
cv_tmp <- data.frame(cvfit_lasso$lambda, cvfit_lasso$cvm, cvfit_lasso$cvsd, cvfit_lasso$cvup, cvfit_lasso$cvlo)
colnames(cv_tmp) <- c("lambda", "mse", "sde", "upper", "lower")

plt_cvglmnet_lasso <- ggplot(cv_tmp, aes(log(lambda), mse)) + 
  geom_vline(xintercept=log(cvfit_lasso$lambda.min), linetype='dashed') + geom_vline(xintercept=log(cvfit_lasso$lambda.1se), linetype='dashed') +
  geom_errorbar(aes(ymin = lower, ymax = upper), col = "grey") + geom_point(col = "red") + 
  xlab("log(Lambda)") + ylab("Mean-squared error") + theme_bw() + theme(legend.position="none")

generate_figs(plt_cvglmnet_lasso, "./plots/plt_cvglmnet_lasso")

weights_lasso_min <- prepare_weights(cvfit_lasso, "lambda.min")
write.csv(weights_lasso_min, "./data/weights_lasso_min.csv")

weights_lasso_reg <- prepare_weights(cvfit_lasso, "lambda.1se")
write.csv(weights_lasso_reg, "./data/weights_lasso_reg.csv")


# Ridge regression

cvfit_ridge <- cv.glmnet(express_mat_sub, as.matrix(out$cartex_score), type.measure = "mse", nfolds = 5, alpha = 0)
cv_tmp <- data.frame(cvfit_ridge$lambda, cvfit_ridge$cvm, cvfit_ridge$cvsd, cvfit_ridge$cvup, cvfit_ridge$cvlo)
colnames(cv_tmp) <- c("lambda", "mse", "sde", "upper", "lower")

plt_cvglmnet_ridge <- ggplot(cv_tmp, aes(log(lambda), mse)) + 
  geom_vline(xintercept=log(cvfit_ridge$lambda.min), linetype='dashed') + geom_vline(xintercept=log(cvfit_ridge$lambda.1se), linetype='dashed') +
  geom_errorbar(aes(ymin = lower, ymax = upper), col = "grey") + geom_point(col = "red") + 
  xlab("log(Lambda)") + ylab("Mean-squared error") + theme_bw() + theme(legend.position="none")

generate_figs(plt_cvglmnet_ridge, "./plots/plt_cvglmnet_ridge")

weights_ridge_min <- prepare_weights(cvfit_ridge, "lambda.min")
write.csv(weights_ridge_min, "./data/weights_ridge_min.csv")

weights_ridge_reg <- prepare_weights(cvfit_ridge, "lambda.1se")
write.csv(weights_ridge_reg, "./data/weights_ridge_reg.csv")


# Elastic net penalty

cvfit_elastic <- cv.glmnet(express_mat_sub, as.matrix(out$cartex_score), type.measure = "mse", nfolds = 5, alpha = 0.5)
cv_tmp <- data.frame(cvfit_elastic$lambda, cvfit_elastic$cvm, cvfit_elastic$cvsd, cvfit_elastic$cvup, cvfit_elastic$cvlo)
colnames(cv_tmp) <- c("lambda", "mse", "sde", "upper", "lower")

plt_cvglmnet_elastic <- ggplot(cv_tmp, aes(log(lambda), mse)) + 
  geom_vline(xintercept=log(cvfit_elastic$lambda.min), linetype='dashed') + geom_vline(xintercept=log(cvfit_elastic$lambda.1se), linetype='dashed') +
  geom_errorbar(aes(ymin = lower, ymax = upper), col = "grey") + geom_point(col = "red") + 
  xlab("log(Lambda)") + ylab("Mean-squared error") + theme_bw() + theme(legend.position="none")

generate_figs(plt_cvglmnet_elastic, "./plots/plt_cvglmnet_elastic")

weights_elastic_min <- prepare_weights(cvfit_elastic, "lambda.min")
write.csv(weights_elastic_min, "./data/weights_elastic_min.csv")

weights_elastic_reg <- prepare_weights(cvfit_elastic, "lambda.1se")
write.csv(weights_elastic_reg, "./data/weights_elastic_reg.csv")

# Quasi-Ridge regression

cvfit_quasi_ridge <- cv.glmnet(express_mat_sub, as.matrix(out$cartex_score), type.measure = "mse", nfolds = 5, alpha = 0.25)
cv_tmp <- data.frame(cvfit_quasi_ridge$lambda, cvfit_quasi_ridge$cvm, cvfit_quasi_ridge$cvsd, cvfit_quasi_ridge$cvup, cvfit_quasi_ridge$cvlo)
colnames(cv_tmp) <- c("lambda", "mse", "sde", "upper", "lower")

plt_cvglmnet_quasi_ridge <- ggplot(cv_tmp, aes(log(lambda), mse)) + 
  geom_vline(xintercept=log(cvfit_quasi_ridge$lambda.min), linetype='dashed') + geom_vline(xintercept=log(cvfit_quasi_ridge$lambda.1se), linetype='dashed') +
  geom_errorbar(aes(ymin = lower, ymax = upper), col = "grey") + geom_point(col = "red") + 
  xlab("log(Lambda)") + ylab("Mean-squared error") + theme_bw() + theme(legend.position="none")

generate_figs(plt_cvglmnet_quasi_ridge, "./plots/plt_cvglmnet_quasi_ridge")

weights_quasi_ridge_min <- prepare_weights(cvfit_quasi_ridge, "lambda.min")
write.csv(weights_quasi_ridge_min, "./data/weights_quasi_ridge_min.csv")

weights_quasi_ridge_reg <- prepare_weights(cvfit_quasi_ridge, "lambda.1se")
write.csv(weights_quasi_ridge_reg, "./data/weights_quasi_ridge_reg.csv")


####################################################################################################
########################################## Compare subsets #########################################
####################################################################################################

out_cartex_630 <- cartex("./data/2021-04-28-cartex-data-logTPM.csv","./data/2021-04-28-metadata.csv","../weights/cartex-630-weights.csv")
out_cartex_200 <- cartex("./data/2021-04-28-cartex-data-logTPM.csv","./data/2021-04-28-metadata.csv","../weights/cartex-200-weights.csv")
out_cartex_84 <- cartex("./data/2021-04-28-cartex-data-logTPM.csv","./data/2021-04-28-metadata.csv","../weights/cartex-84-weights.csv")

colorRampPalette(c("lightgrey","lightblue","mediumblue"))(4)
# scale_fill_manual(values = c("0" = "#D3D3D3", "11" = "#B9D6DF", "15" = "#7390DD", "21" = "#0000CD"))
# scale_fill_manual(values = colorRampPalette(c("lightgrey","lightblue","mediumblue"))(4))
# scale_fill_manual(values = c("#D3D3D3", "#B9D6DF", "#7390DD", "#0000CD"))

ggplot(out_cartex_630, aes(x=factor(CAR,level=c("Control","CD19","HA")),y=cartex_score)) + 
  geom_point(stat = "identity",aes(fill = Days, color=Days)) + 
  scale_fill_manual(values = c("0" = "#D3D3D3", "11" = "#B9D6DF", "15" = "#7390DD", "21" = "#0000CD"))


# Days = factor(out_cartex_630$Day)
Days = factor(as.character(out_cartex_630$Day))
# out_cartex_630$Day <- factor(as.character(out_cartex_630$Day))
plt_CARTEx_630 <- ggplot(out_cartex_630, aes(x=factor(CAR,level=c("Control","CD19","HA")),y=cartex_score)) + geom_boxplot() + geom_point(aes(color=Days)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Control','CD19')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('HA','CD19')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('HA','Control')), label = "p.signif", label.y = 2.1) + 
  xlab("CAR") + ylab("CARTEx score 630") + theme_bw()

generate_figs(plt_CARTEx_630, "./plots/plt_CARTEx_630", c(5,4))

plt_CARTEx_200 <- ggplot(out_cartex_200, aes(x=factor(CAR,level=c("Control","CD19","HA")),y=cartex_score)) + geom_boxplot() + geom_point(aes(color=Days)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Control','CD19')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('HA','CD19')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('HA','Control')), label = "p.signif", label.y = 2.1) +
  xlab("CAR") + ylab("CARTEx score 200") + theme_bw()

generate_figs(plt_CARTEx_200, "./plots/plt_CARTEx_200", c(5,4))

plt_CARTEx_84 <- ggplot(out_cartex_84, aes(x=factor(CAR,level=c("Control","CD19","HA")),y=cartex_score)) + geom_boxplot() + geom_point(aes(color=Days)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('Control','CD19')), label = "p.signif", label.y = 1.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('HA','CD19')), label = "p.signif", label.y = 1.8) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c('HA','Control')), label = "p.signif", label.y = 2.1) +
  xlab("CAR") + ylab("CARTEx score 84") + theme_bw()

generate_figs(plt_CARTEx_84, "./plots/plt_CARTEx_84", c(5,4))









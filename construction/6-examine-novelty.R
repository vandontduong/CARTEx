
setwd("/oak/stanford/groups/cmackall/vandon/CARTEx/construction")
set.seed(123)
source("/oak/stanford/groups/cmackall/vandon/CARTEx/cartex-utilities.R")


c5_results <- read.csv('./data/cellformatica-C5-result.csv', header = TRUE)


plot(c5_results$effect_score, c5_results$novelty_score)

plot_effect_novelty_density <- ggplot(c5_results, aes(x = effect_score, y = novelty_score) ) +
  geom_density_2d() + geom_point() +
  theme_classic() + labs(x = "Effect Score", y = "Novelty Score")
generate_figs(plot_effect_novelty_density, './plots/plot_effect_novelty_density', c(2.5,2))

plot_effect_novelty <- ggplot(c5_results, aes(x = effect_score, y = novelty_score)) +
  geom_count() + scale_size_area() +
  theme_classic() + labs(x = "Effect Score", y = "Novelty Score")
generate_figs(plot_effect_novelty, './plots/plot_effect_novelty', c(4,2))

plot_effect_novelty_color <- ggplot(c5_results, aes(x = effect_score, y = novelty_score)) +
  geom_count(aes(color = ifelse(effect_score >= 60 & novelty_score >= 50, "purple",
                                ifelse(effect_score >= 60, "blue",
                                       ifelse(novelty_score >= 50, "red", "black"))))) +
  scale_size_area() +  scale_color_identity() +
  theme_classic() + labs(x = "Effect Score", y = "Novelty Score")


# examine genes with novelty score >= 50
subset(c5_results, novelty_score >= 50)$approved_symbol








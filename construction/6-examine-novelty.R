
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

plot_effect_novelty_color <- ggplot(c5_results, aes(x = effect_score, y = novelty_score)) +
  geom_count(aes(color = ifelse(effect_score >= 80 | novelty_score >= 50, "purple","black"))) +
  scale_size_area() +  scale_color_identity() +
  theme_classic() + labs(x = "Effect Score", y = "Novelty Score")

# examine genes with novelty score >= 50
subset(c5_results, novelty_score >= 50)$approved_symbol

# examine genes with effect score >= 50
subset(c5_results, effect_score >= 70)$approved_symbol


data_subset <- c5_results[c("approved_symbol", "effect_score", "confidence_score", "novelty_score")]
filtered_data <- data_subset[data_subset$novelty_score >= 50 | data_subset$effect_score >= 80, ]

plot_effect_novelty_labeled <- ggplot(data_subset, aes(x = effect_score, y = novelty_score)) +
  geom_count() + scale_size_area() + 
  geom_text_repel(data = filtered_data, 
                  aes(label = approved_symbol),
                  size = 2, max.overlaps = 50) +
  theme_classic() + labs(x = "Effect Score", y = "Novelty Score")
generate_figs(plot_effect_novelty_labeled, './plots/plot_effect_novelty_labeled', c(4,2))



plot_effect_novelty_color_labeled <- ggplot(data_subset, aes(x = effect_score, y = novelty_score)) +
  geom_count(aes(color = ifelse(effect_score >= 80 | novelty_score >= 50, "purple","black"))) +
  scale_size_area() +  scale_color_identity() +
  geom_text_repel(data = filtered_data,
                  min.segment.length = 0, segment.curvature = -0.1,
                  aes(label = approved_symbol),
                  size = 2, max.overlaps = 50) +
  theme_classic() + labs(x = "Effect Score", y = "Novelty Score")
generate_figs(plot_effect_novelty_color_labeled, './plots/plot_effect_novelty_color_labeled', c(4,2))






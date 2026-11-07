library(ggplot2)
library(dplyr)
library(tidyr)

urban_long <- read.csv("Urban.csv") %>%
  pivot_longer(cols = -c(Indvi, Group), names_to = "Condition", values_to = "Value")

rural_long <- read.csv("Rural.csv") %>%
  pivot_longer(cols = -c(Indvi, Group), names_to = "Condition", values_to = "Value")

plot_data <- bind_rows(urban_long, rural_long)

plot_data$Condition <- factor(plot_data$Condition,
                              levels = c("Baseline", "OTA", "MT", "Washout"))
plot_data$Group <- factor(plot_data$Group, levels = c("Urban", "Rural"))

p <- ggplot(plot_data, aes(x = Condition, y = Value)) +
  geom_violin(aes(color = Group), fill = NA, trim = FALSE, scale = "width", linewidth = 0.5) +
  geom_line(aes(group = Indvi), color = "grey80", linewidth = 0.3, alpha = 0.7) +
  geom_point(aes(color = Group), size = 1.4, shape = 21, fill = "white", stroke = 0.6) +
  facet_wrap(~ Group, scales = "free_x", drop = TRUE) +
  scale_color_manual(values = c("Urban" = "#FF6000", "Rural" = "#5CA5B4")) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  labs(y = "Firing rate (Hz)", x = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14)
  )

ggsave("AA_Urban_Rural.pdf", p, width = 10, height = 8, units = "cm")
cat("saved AA_Urban_Rural.pdf (10 x 8 cm)\n")

library(ggplot2)
library(latex2exp)

# Load needed library
library(stats)

# Define tau_c computation
compute_tau_c <- function(constant) {
  b=pnorm(constant,0,1)-pnorm(-constant,0,1)
  sigma_psi_squared=((1-(constant*dnorm(constant,0,1)-(-constant)*dnorm(-constant,0,1))/(pnorm(constant,0,1)-pnorm(-constant,0,1)))*(pnorm(constant,0,1)-pnorm(-constant,0,1))+constant^2*(1-b))
  tau=b^2/sigma_psi_squared
  tau
  return(1/tau)
}

# Example: compute for c = 1.5
c_value <- 1.5
result <- compute_tau_c(0.1)







data <- NULL
for(x in  seq(0, 1.5, by=0.01)){
  data <- rbind(data, c(x,(compute_tau_c(x))))
}
data <- data.frame(data)
names(data) <- c("d","ARE")
# Create plot
plot2<-ggplot(data, aes(x = d, y = ARE, group = 1)) +
  geom_line(
    aes(color = "ARE curve"),
    size = 1.5
  ) +
  ylab(TeX("τ$_c^{-1}$")) +
  xlab(TeX("$c$")) +
  ggtitle("") +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22)
  ) +
  scale_y_continuous(limits = c(1, 1.8)) +  # adjust as needed
  scale_x_continuous(
    limits = c(min(data$d), NA),
    breaks = unique(c(min(data$d), pretty(data$d)))
  ) +
  theme(
    legend.position = c(0.80, 0.9),
    legend.text = element_text(size = 18)
  ) +
  scale_color_manual(
    name = "",
    values = c("ARE curve" = "black"),
    labels = c("ARE curve" = "Huber-type function")
  ) + theme(legend.position="none")


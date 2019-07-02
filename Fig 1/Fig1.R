
library(dplyr)
library(ggplot2)


#load the data base
data.plot <- read.csv ("./database_Fig1.csv")

#select the method - interaction 
type.method.interaction <- interaction(data.plot$Ecoli.cutoff, data.plot$index)



fig1 <- ggplot(data = data.plot, aes(x = Host.cutoff, y = Coef.Estimate, colour =type.method.interaction, shape=type.method.interaction)) + 
  geom_point (size=1.5)+
  theme_classic() +
  scale_shape_manual("Similarity cut-off & Index", values = c(17,16,17,16)) +
  scale_colour_discrete("Similarity cut-off & Index") +
  labs(x= "Host association thresholds", y="Regression coefficient") +
  geom_smooth (method = "loess", size = 0.5, se=F, span=0.5)

plot(fig1)

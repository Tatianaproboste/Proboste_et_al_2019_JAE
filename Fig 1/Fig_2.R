
library(dplyr)
library(ggplot2)

load("./test.unweigthed.v3.rda")
test <-test.unweigthed.v3
  
  plot.data <- data.frame(Ecoli.cutoff = as.factor(test$ecoli.cutoff.raw),
                        Host.cutoff = test$host.cutoff.raw,
                        index = test$host.index.raw,
                        Coef = test$degree.coef.raw %>%
                          filter(Parameter=='Ecoli.degree'),
                        Coef.Year.2015 = test$degree.coef.raw %>%
                          filter(Parameter=='Year.2015'),
                        Coef.betweeness = test$betweeness.coef.raw %>%
                          filter(Parameter=='Ecoli.betweenness'),
                        Coef.betweeness.Year.2015 = test$betweeness.coef.raw %>%
                          filter(Parameter=='Ecoli.betweenness*Year.2015'),
                        Coef.eigencent = test$eigencent.coef.raw %>%
                          filter(Parameter=='Ecoli.evcent'),
                        Coef.eigencent.Year.2015 = test$eigencent.coef.raw %>%
                          filter(Parameter=='Year.2015'),
                        Coef.degree.Year.2015 = test$degree.coef.raw %>%
                          filter(Parameter=='Ecoli.degree*Year.2015'))


type.method.interaction <- interaction(plot.data$Ecoli.cutoff, plot.data$index)

# Compute the number of types and methods
nb.cut.off <- nlevels(plot.data$Ecoli.cutoff)
nb.index <- nlevels(plot.data$index)

#camparison betweeen social index   
cols <- c('#A3245B', '#597D79','#E7A582', '#ACC66B')

fig2 <- ggplot(data = plot.data, aes(x = Host.cutoff, y = Coef.Estimate, colour =type.method.interaction, shape=type.method.interaction)) + 
  geom_point (size=1.5)+
  theme_classic() +
  scale_shape_manual("Similarity cut-off & Index", values = c(17,16,17,16)) +
  scale_colour_manual("Similarity cut-off & Index", values = cols) +
  labs(x= "Host association thresholds", y="Regression coefficient") +
  geom_smooth (method = "loess", size = 0.75, se=F, span=0.6)

plot(fig2)


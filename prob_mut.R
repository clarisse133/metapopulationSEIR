#packages
library(ggplot2)

#N = numero de individuos infectados pela cepa
#P = probablidade de mutacao
#I = total de infectados(fixo)

I <- 10000

P <- seq(0, 1, by=0.001)

N <- P*I

dados <- data.frame(P,N)

#plot

ggplot(dados, aes(x=P, y=N)) +
  geom_line(color="red", linewidth = 0.6) +
  labs(title="Numero de infectados em funcao da prob. de mutacao",
  x="probabilidade de mutacao", 
  y="infectados com mutacao") +
  theme_minimal(base_size = 14)


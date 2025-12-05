I <- 1000
P <- seq(0, 1, 0.00001)

prob_1mut <- 1 - (1 - P)^I

data <- data.frame(P, prob_mut)

ggplot(data, aes(P, prob_mut)) +
  geom_line() +
  labs(x = "Prob. mutação",
       y = "Probabilidade de ocorrer pelo menos 1 mutação",
       title = "Probabilidade de haver pelo menos uma mutação") +
  theme_minimal()

#quanto maior for o n de eventos I, mais rapido o y se aproxima de 1

#packages
library(ggplot2)
library(deSolve)


 mut_model0 <- function(tempo, estado, param) {
  with(as.list(c(estado, param)), {
    dS0 = -beta0*S0*I0
    dE0 = beta0*S0*I0 - E0*sigma0
    dI0 = E0*sigma0 - I0*gama0
    dR0 = I0*gama0   
    
    return(list(c(dS0, dE0, dI0, dR0)))
    
  })
}

param0 <- c(
  beta0 = 0.6,
  sigma0 = 0.3,
  gama0 = 0.2,
  N0 = 1000
)

estado0 <- c(
  S0 = 999,
  E0 = 0,
  I0 = 1,
  R0 = 0
)

tempo0 <- seq(0, 30, by = 1)

output0 <- ode(y = estado0, times = tempo0, func = mut_model0, parms = param0)
out_df0 <- as.data.frame(output0)

plot <- ggplot(out_df0, aes(x = time)) +
  geom_line(linewidth = 1.2, aes(y = S0, color = "Susceptible(S0)")) +
  geom_line(linewidth = 1.2, aes(y = E0, color = "Exposed(E0)")) +
  geom_line(linewidth = 1.2, aes(y = I0, color = "Infected(I0)")) +
  geom_line(linewidth = 1.2, aes(y = R0, color = "Recovered(R0)")) +
  labs(
    title = "Modelo para implementar mutação viral",
    x = "Tempo (dias)",
    y = "População",
    color = "Estado"
  ) +
  scale_color_manual(values = c(
    "Susceptible(S0)" = "green",
    "Exposed(E0)" = "purple",
    "Infected(I0)" = "#ff3333",
    "Recovered(R0)" = "black"
  )) +
  theme_minimal()
  
  
  plot
  

rbinom(1, 1000, 0.001)

mut_model1 <- function(tempo1, estado1, param1) {
  with(as.list(c(estado1, param1)), {
    
    dS1 = -S1*(beta0*I0 + beta1*I1)
    dE1 = S1*(beta0*I0 + beta1*I1) - E1*(sigma0 + sigma1)
    dI1 = E1*(sigma0 + sigma1) - I1*(gama0 + gama1)
    dR1 = I1*(gama0 + gama1)
    
    return(list(c(dS1, dE1, dI1, dR1)))
    
  })
}

result <- rbinom(1, 1000, 0.001)

param1 <- c(
  param0,
  beta1 = 0.9,
  sigma1 = 0.7,
  gama1 = 0.3,
  N1 = 1000
)

estado1 <- c(
  S1 = 1000 - result,
  E1 = 0,
  I1 = result,
  R1 = 0
)

tempo1 <- seq(0, 30, by = 1)

output1 <- ode(y = estado1, times = tempo1, func = mut_model1, parms = param1)
out_df1 <- as.data.frame(output1)

plot1 <- ggplot(out_df1, aes(x = time)) +
  geom_line(linewidth = 1.2, aes(y = S1, color = "Susceptible(S1)")) +
  geom_line(linewidth = 1.2, aes(y = E1, color = "Exposed(E1)")) +
  geom_line(linewidth = 1.2, aes(y = I1, color = "Infected(I1)")) +
  geom_line(linewidth = 1.2, aes(y = R1, color = "Recovered(R1)")) +
  labs(
    title = "Modelo para implementar mutação viral",
    x = "Tempo (dias)",
    y = "População",
    color = "Estado"
  ) +
  scale_color_manual(values = c(
    "Susceptible(S1)" = "green",
    "Exposed(E1)" = "purple",
    "Infected(I1)" = "#ff3333",
    "Recovered(R1)" = "black"
  )) +
  theme_minimal()

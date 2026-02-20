#----packages----
library(deSolve) 
library(devtools)
library(ggplot2)
library(patchwork)

#' @param betta transmission rate
#' @param gamma recovery rate
#' @param M migration rate
#' @param sigma transition rate from exposed to infected
#' @param alpha vaccination rate
#' @param d death rate
#' @param a immunity time
#' 
#' @export
#' 

#---function---

#adicionar migração de outra forma

metap_mut <- function(time_mut, state_mut, parameters_mut) {
  with(as.list(c(state_mut, parameters_mut)), {
    
    dS = -S*(beta1*I1 + beta2*I2 + alpha) + a*R
    dE1 = -sigma1*E1 + beta1*S*I1
    dI1 = sigma1*E1 - gamma1*I1 - d1*I1
    dR = -a*R + gamma1*I1 + gamma2*I2
    
    dE2 = -sigma2*E2 + beta2*S*I2
    dI2 = sigma2*E2 - gamma2*I2 - d2*I2
    
    dV = -V*(betav1*I1 + betav2*I2) + gammav1*Iv1 + gammav2*Iv2 + alpha*S
    dEv1 = betav1*V*I1 - sigmav1*Ev1
    dEv2 = betav2*V*I2 - sigmav2*Ev2
    dIv1 = sigmav1*Ev1 - gammav1*Iv1
    dIv2 = sigmav2*Ev2 - gammav2*Iv2 - dv2*Iv2
    
    dO1 = d1*I1
    dO2 = d2*I2
    dOv = dv2*Iv2

      return(list(c(dS, dE1, dI1, dR, dE2, dI2, dV, dEv1, dEv2, dIv1, dIv2, dO1, dO2, dOv)))
    
    
  })
}
#---parameter values---

#adicionar migração(M) depois 
#vacinados possuem mais resistencia ao virus 1 do que o 2

parameters_mut <- c(
  beta1= 0.04, #0.4 - 0.6 
  beta2= 0.06,
  betav1= 0.01,
  betav2= 0.03,
  sigma1= 1/5, #5 - 6 dias
  sigma2= 1/4,
  sigmav1= 1/3,
  sigmav2= 1/3,
  gamma1= 1/6,
  gamma2= 1/7,
  gammav1= 1/5,
  gammav2= 1/5,
  alpha= 0.05,
  d1= 0.02,
  d2= 0.035,
  dv2= 0.01,
  a= 1/150 
  

) 
#inicialmente variar beta e beta2
#analise de parametros que potencialmente geram uma bifurcação

#---initial conditions---

state_mut <- c(
  S=9998,
  E1=0,
  I1=1,
  R=0,

  E2=0,
  I2=1,
  
  V=0,
  Ev1=0,
  Ev2=0,
  Iv1=0,
  Iv2=0,

  O1= 0,
  O2= 0,
  Ov2= 0
  
)


time_mut <- seq(0, 200, by=1)

#--- solve the EDO ---

output_mut <- ode(y = state_mut, times = time_mut, func = metap_mut, parms = parameters_mut)
head(output_mut)

#--- data frame ---
output_mutdf <- as.data.frame(output_mut)

g1 <- ggplot(output_mutdf, aes(x = time)) +
  geom_line(linewidth = 0.75, aes(y = S, color = "Susceptible")) +
  geom_line(linewidth = 0.75, aes(y = E1, color = "Exposed 1")) +
  geom_line(linewidth = 0.75, aes(y = I1, color = "Infected 1")) +
  geom_line(linewidth = 0.75, aes(y = R, color = "Recovered")) +
  geom_line(linewidth = 0.75, aes(y = O1, color = "Dead 1")) +
  labs(x = "Time(days)", y = "Number os individuals", title = "SEIRS Dynamics for virus 1") +
  scale_color_manual(values = c("Susceptible"="#2f5aa8", "Exposed 1"="#f6ed24", "Infected 1"="#ed1b25", "Recovered"="#60af46", "Dead 1"="black")) +
  theme_minimal() +
  theme(legend.title = element_blank())

g2 <- ggplot(output_mutdf, aes(x = time)) +
  geom_line(linewidth = 0.75, aes(y = S, color = "Susceptible")) +
  geom_line(linewidth = 0.75, aes(y = E2, color = "Exposed 2")) +
  geom_line(linewidth = 0.75, aes(y = I2, color = "Infected 2")) +
  geom_line(linewidth = 0.75, aes(y = O2, color = "Dead 2")) +
  geom_line(linewidth = 0.75, aes(y = R, color = "Recovered")) +
  labs(x = "Time(days)", y = "Number of population", title = "SEIRS Dynamics for virus 2") +
  scale_color_manual(values = c("Susceptible"="#2f5aa8", "Exposed 2"="#f6ed24", "Infected 2"="#ed1b25", "Recovered"="#60af46", "Dead 2"="black")) +
  theme_minimal() +
  theme(legend.title = element_blank())

g3 <- ggplot(output_mutdf, aes(x = time)) +
  geom_line(linewidth = 0.75, aes(y = V, color = "Vaccinated")) +
  geom_line(linewidth = 0.75, aes(y = Ev1, color = "Exposed v1")) +
  geom_line(linewidth = 0.75, aes(y = Iv1, color = "Infected v1")) +
  labs(x = "Time(days)", y = "Number of population", title = "SEIRS Dynamics for vaccinated virus 1") +
  scale_color_manual(values = c("Vaccinated"="#a41a44", "Exposed v1"="#f6ed24", "Infected v1"="#ed1b25")) +
  theme_minimal() +
  theme(legend.title = element_blank())

g4 <- ggplot(output_mutdf, aes(x = time)) +
  geom_line(linewidth = 0.75, aes(y = V, color = "Vaccinated")) +
  geom_line(linewidth = 0.75, aes(y = Ev2, color = "Exposed v2")) +
  geom_line(linewidth = 0.75, aes(y = Iv2, color = "Infected v2")) +
  geom_line(linewidth = 0.75, aes(y = Ov2, color = "Dead v2")) +
  labs(x = "Time(days)", y = "Number of population", title = "SEIRS Dynamics for vaccinated virus 2") +
  scale_color_manual(values = c("Vaccinated"="#a41a44", "Exposed v2"="#f6ed24", "Infected v2"="#ed1b25", "Dead v2"="black")) +
  theme_minimal() +
  theme(legend.title = element_blank())

g1/g2/g3/g4

  #Dividir em 4 gráficos:Um para não vacinados e outro para vacinados
#Devo criar compartimento de n obitos para cada virus? ex: O1 O2 Ov2
#revisar os valores dos parametros
#repensar a ideia da migração entre A e B
#posteriormente: usar alphas distintos para A e B
#usar parametros para experimentar e observar, mas tambem usar valores reais
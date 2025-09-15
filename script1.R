
#----packages----
library(deSolve)
library(devtools)
library(ggplot2)

#' @param betta transmission rate
#' @param gamma recovery rate
#' @param un birth rate
#' @param um mortality rate
#' @param M migration rate
#' @param sigma transition rate from exposed to infected
#' @param alpha vaccination rate
#' 
#' @export

#--- Ma = Mb = M ---
#--- un = um ---

#---function---

metap <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
  
  dSa = -betta*Sa*Ia + un*(Sa+Ea+Ia+Ra) - Sa*um + M*(Sb - Sa)
  dEa = betta*Sa*Ia - Ea*(um + sigma) + M*(Eb - Ea)
  dIa = sigma*Ea - Ia*(um + gamma) + M*(Ib - Ia)
  dRa = gamma*Ia - Ra*um + M*(Rb- Ra)
     
  dSb = -betta*Sb*Ib + + un*(Sb+Eb+Ib+Rb) - Sb*um + M*(Sa - Sb)
  dEb = betta*Sb*Ib - Eb*(um + sigma) + M*(Ea - Eb)
  dIb = sigma*Eb - Ib*(um + gamma) + M*(Ia - Ib)
  dRb = gamma*Ib - Rb*um + M*(Ra- Rb)
  
   return(list(c(dSa, dEa, dIa, dRa, dSb, dEb, dIb, dRb)))


  })
}
#---parameter values---

parameters <- c(
  betta=0.4,
  gamma=0.1,
  un=0.02,
  um=0.02,
  M=0.4,
  sigma=0.8
) 

#---initial conditions---

state <- c(
  Sa=0.9,
  Ea=0,
  Ia=0.1,
  Ra=0,
  
  Sb=0.9,
  Eb=0,
  Ib=0.1,
  Rb=0
  
)

# Experimentando com outro estado inicial
state <- c(
  Sa=0.1-0.01,
  Ea=0,
  Ia=0.01,
  Ra=0.9,
  
  Sb=0.1-0.01,
  Eb=0,
  Ib=0.01,
  Rb=0.9
  
)

time <- seq(0, 100, by=1)
#time <- seq(0, 1000, by=1)
time <- seq(0, 250, by=1)

#--- solve the EDO ---

output <- ode(y = state, times = time, func = metap, parms = parameters)
head(output)

#---plot A---

out_df <- as.data.frame(output)
out_df$sum_a = out_df$Sa + out_df$Ea + out_df$Ia + out_df$Ra
out_df$sum_b = out_df$Sb + out_df$Eb + out_df$Ib + out_df$Rb

ggplot(out_df, aes(x=time)) +
  geom_line(size = 0.7, aes(y = sum_a, color = "Sum A")) +
  geom_line(size = 0.7, aes(y = sum_b, color = "Sum B")) 
  

ggplot(out_df, aes(x = time)) +
  geom_line(size = 0.7, aes(y = Sa, color = "Susceptible A")) +
  geom_line(size = 0.7, aes(y = Ea, color = "Exposed A")) +
  geom_line(size = 0.7, aes(y = Ia, color = "Infected A")) +
  geom_line(size = 0.7, aes(y = Ra, color = "Recovered A")) +
  labs(x = "Time(days)", y = "Proportion of population", title = "SEIR Dynamics for the population A") +
  scale_color_manual(values = c("Susceptible A" = "blue", "Exposed A" = "orange", "Infected A" = "red", "Recovered A" = "green")) +
  theme_minimal() +
  theme(legend.title = element_blank())

#---plot B---

ggplot(out_df, aes(x = time)) +
  geom_line(size = 0.7, aes(y = Sb, color = "Susceptible B")) +
  geom_line(size = 0.7, aes(y = Eb, color = "Exposed B")) +
  geom_line(size = 0.7, aes(y = Ib, color = "Infected B")) +
  geom_line(size = 0.7, aes(y = Rb, color = "Recovered B")) +
  labs(x = "Time(days)", y = "Proportion of population", title = "SEIR Dynamics for the population B") +
  scale_color_manual(values = c("Susceptible B" = "blue", "Exposed B" = "orange", "Infected B" = "red", "Recovered B" = "green")) +
  theme_minimal() +
  theme(legend.title = element_blank())


### Introduzir Vacina leaky

metap_leaky <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dSa = -betta*Sa*Ia + un*(Sa+Ea+Ia+Ra) - Sa*um + M*(Sb - Sa) - Sa*alpha
    dEa = betta*Sa*Ia - Ea*(um + sigma) + M*(Eb - Ea)
    dIa = sigma*Ea - Ia*(um + gamma) + M*(Ib - Ia)
    dRa = gamma*Ia - Ra*um + M*(Rb- Ra)
    
    dVa = Sa*alpha - bettav*Va*Ia - um*Va + M*(Vb - Va)
    dEav = 
    dIav =
    dRav = 
    
    dSb = -betta*Sb*Ib + + un*(Sb+Eb+Ib+Rb) - Sb*um + M*(Sa - Sb)
    dEb = betta*Sb*Ib - Eb*(um + sigma) + M*(Ea - Eb)
    dIb = sigma*Eb - Ib*(um + gamma) + M*(Ia - Ib)
    dRb = gamma*Ib - Rb*um + M*(Ra- Rb)
    
    return(list(c(dSa, dEa, dIa, dRa, dSb, dEb, dIb, dRb)))
    
    
  })
}

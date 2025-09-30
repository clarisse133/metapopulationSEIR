
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
     
  dSb = -betta*Sb*Ib + un*(Sb+Eb+Ib+Rb) - Sb*um + M*(Sa - Sb)
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

metap_leaky <- function(time_leaky, state_leaky, parameters_leaky) {
  with(as.list(c(state_leaky, parameters_leaky)), {
    
    dSa = - betta*Sa*Ia - betta_vn * Sa*Iav + un*(Sa+Ea+Ia+Ra+ Va+Eav+Iav+Rav) - Sa*um + M*(Sb - Sa) - Sa*alpha_a
    dEa = betta*Sa*Ia + betta_vn * Sa*Iav - Ea*(um + sigma) + M*(Eb - Ea)
    dIa = sigma*Ea - Ia*(um + gamma) + M*(Ib - Ia)
    dRa = gamma*Ia - Ra*um + M*(Rb- Ra)
    
    dVa = Sa*alpha_a - betta_v * Va*Ia - betta_vi * Va* Iav - um*Va + M*(Vb - Va)
    dEav =  betta_v * Va*Ia + betta_vi * Va* Iav - Eav *(um+sigma_v) + M*(Ebv - Eav) 
    dIav =  sigma_v * Eav - Iav *(um+ gamma_v) + M*(Ibv - Iav)
    dRav =  gamma_v * Iav - Rav*um + M*(Rbv - Rav)

    dSb = - betta*Sb*Ib - betta_vn * Sb*Ibv + un*(Sb+Eb+Ib+Rb+Vb+Ebv+Ibv+Rbv) - Sb*um + M*(Sa - Sb) - Sb*alpha_b
    dEb = betta*Sb*Ib + betta_vn * Sb*Ibv - Eb*(um + sigma) + M*(Ea - Eb)
    dIb = sigma*Eb - Ib*(um + gamma) + M*(Ia - Ib)
    dRb = gamma*Ib - Rb*um + M*(Ra- Rb)
    
    dVb = Sb*alpha_b - betta_v * Vb*Ib - betta_vi * Vb* Ibv - um*Vb + M*(Va - Vb)
    dEbv =  betta_v * Vb*Ib + betta_vi * Vb* Ibv - Ebv *(um+sigma_v) + M*(Eav - Ebv) 
    dIbv =  sigma_v * Ebv - Ibv *(um+ gamma_v) + M*(Iav - Ibv)
    dRbv =  gamma_v * Ibv - Rbv*um + M*(Rav - Rbv)
    
#    dSb = -betta*Sb*Ib + + un*(Sb+Eb+Ib+Rb) - Sb*um + M*(Sa - Sb)
#    dEb = betta*Sb*Ib - Eb*(um + sigma) + M*(Ea - Eb)
#    dIb = sigma*Eb - Ib*(um + gamma) + M*(Ia - Ib)
#    dRb = gamma*Ib - Rb*um + M*(Ra- Rb)
    
    return(list(c(dSa, dEa, dIa, dRa, dSb, dEb, dIb, dRb, dVa, dEav, dIav, dRav, dVb, dEbv, dIbv, dRbv)))
    
    
  })
}

parameters_leaky <- c(
  
  #---vaccinated---
 alpha_a=0.75,
 alpha_b=0.75,
 betta_v=3.0,
 betta_vi=2.5,
 sigma_v=0.4,
 gamma_v=0.2,
 
 #---unvaccinated---
 betta=3.0,
 betta_vn=2.5,
 sigma=0.4,
 gamma=0.2,
 M=0.1,
 un=0.02,
 um=0.02
  
)

state_leaky <- c(
  Sa=0.8,
  Ea=0,
  Ia=0.2,
  Ra=0,
  
  Sb=0.8,
  Eb=0,
  Ib=0.2,
  Rb=0,
  
  Va=0,
  Eav=0,
  Iav=0,
  Rav=0,
  
  Vb=0,
  Ebv=0,
  Ibv=0,
  Rbv=0
)

time_leaky <- seq(0, 500, by=1)

#---solve the edo---
output_leaky <- ode(y = state_leaky, times = time_leaky, func = metap_leaky, parms = parameters_leaky)
head(output_leaky)

#---plot A---
outdf_leaky <- as.data.frame(output_leaky)
ggplot(outdf_leaky, aes(x = time)) +
  geom_line(size = 0.7, aes(y = Sa, color = "dSa")) +
  geom_line(size = 0.7, aes(y = Ea, color = "dEa")) +
  geom_line(size = 0.7, aes(y = Ia, color = "dIa")) +
  geom_line(size = 0.7, aes(y = Ra, color = "dRa")) +
  labs(x = "Time(days)", y = "Proportion of population", title = "SEIR Dynamics for the population A novaccinated") +
  scale_color_manual(values = c("dSa"="blue", "dEa"="orange", "dIa"= "red","dRa"="green" )) +
  theme_minimal() +
  theme(legend.title = element_blank())

ggplot(outdf_leaky, aes(x = time)) +
  geom_point(size = 0.6, aes(y = Va, color = "dVa")) +
  geom_point(size = 0.6, aes(y = Eav, color = "dEav")) +
  geom_point(size = 0.6, aes(y = Iav, color = "dIav")) +
  geom_point(size = 0.6, aes(y = Rav, color = "dRav")) +
  labs(x = "Time(days)", y = "Proportion of population", title = "SEIR Dynamics for the population A vaccinated") +
  scale_color_manual(values = c("dVa"="blue","dEav"="orange", "dIav"="red", "dRav"="green" )) +
  theme_minimal() +
  theme(legend.title = element_blank())

library(tidyverse)
outdf_leaky %>%
  mutate(Na = Sa+Ea+Ia+Ra+Va+Eav+Iav+Rav) %>%
  mutate(Nb = Sb+Eb+Ib+Rb+Vb+Ebv+Ibv+Rbv) -> df2
  

#---plot B---
outdf_leaky <- as.data.frame(output_leaky)
ggplot(outdf_leaky, aes(x = time)) +
  geom_line(size = 0.7, aes(y = Sb, color = "dSb")) +
  geom_line(size = 0.7, aes(y = Eb, color = "dEb")) +
  geom_line(size = 0.7, aes(y = Ib, color = "dIb")) +
  geom_line(size = 0.7, aes(y = Rb, color = "dRb")) +
  geom_point(size = 0.6, aes(y = Vb, color = "dVb")) +
  geom_point(size = 0.6, aes(y = Ebv, color = "dEbv")) +
  geom_point(size = 0.6, aes(y = Ibv, color = "dIbv")) +
  geom_point(size = 0.6, aes(y = Rbv, color = "dRbv")) +
  labs(x = "Time(days)", y = "Proportion of population", title = "SEIR Dynamics for the population B") +
  scale_color_manual(values = c("dSb"="blue", "dEb"="orange", "dIb"= "red","dRb"="green", "dVb"="blue","dEbv"="orange", "dIbv"="red", "dRbv"="green" )) +
  theme_minimal() +
  theme(legend.title = element_blank())

# 1 - comecar com os valores de beta iguais aos valores utilizados anteriormente
# 1 - betta = betta_v = betta_vi = betta_vn
# gamma_v = gamma
# alpha_a = alpha_b

# 2 - betta_vi e betta_vn menores que betta_v
# Fazer gráficos. O que acontece, o que vc observa?

# 3 - A partir do 1, alterar o valor de gamma_v para gamma_v > gamma
# Fazer graficos. O que vc observa?

# 4 - A partir do 1, experimentar alpha_a < alpha_b
# Fazer gráficos

# Em aberto: Podemos pensar ema avaliação de 
# vacinas com bloqueio de infecção, bloqueio de transmissao e redução de morbi-mortalidade - pensar mais a frente



# Paa introduzir desigualdade, pode ser com alphas distintos ou condicoes iniciais distintas

# alphas distintos

parameters_leaky <- c(
  
  #---vaccinated---
  alpha_a=0.75,
  alpha_b=0.05,
  betta_v=3.0,
  betta_vi=2.5,
  sigma_v=0.4,
  gamma_v=0.2,
  
  #---unvaccinated---
  betta=3.0,
  betta_vn=2.5,
  sigma=0.4,
  gamma=0.2,
  M=0.1,
  un=0.02,
  um=0.02
  
)

state_leaky <- c(
  Sa=0.8,
  Ea=0,
  Ia=0.2,
  Ra=0,
  
  Sb=0.8,
  Eb=0,
  Ib=0.2,
  Rb=0,
  
  Va=0,
  Eav=0,
  Iav=0,
  Rav=0,
  
  Vb=0,
  Ebv=0,
  Ibv=0,
  Rbv=0
)

time_leaky <- seq(0, 100, by=1)

#---solve the edo---
output_leaky <- ode(y = state_leaky, times = time_leaky, func = metap_leaky, parms = parameters_leaky)
head(output_leaky)

#---plot A---
outdf_leaky <- as.data.frame(output_leaky)
ggplot(outdf_leaky, aes(x = time)) +
  geom_line(size = 0.7, aes(y = Sa, color = "dSa")) +
  geom_line(size = 0.7, aes(y = Ea, color = "dEa")) +
  geom_line(size = 0.7, aes(y = Ia, color = "dIa")) +
  geom_line(size = 0.7, aes(y = Ra, color = "dRa")) +
  labs(x = "Time(days)", y = "Proportion of population", title = "SEIR Dynamics for the population A novaccinated") +
  scale_color_manual(values = c("dSa"="blue", "dEa"="orange", "dIa"= "red","dRa"="green" )) +
  theme_minimal() +
  theme(legend.title = element_blank())

ggplot(outdf_leaky, aes(x = time)) +
  geom_point(size = 0.6, aes(y = Va, color = "dVa")) +
  geom_point(size = 0.6, aes(y = Eav, color = "dEav")) +
  geom_point(size = 0.6, aes(y = Iav, color = "dIav")) +
  geom_point(size = 0.6, aes(y = Rav, color = "dRav")) +
  labs(x = "Time(days)", y = "Proportion of population", title = "SEIR Dynamics for the population A vaccinated") +
  scale_color_manual(values = c("dVa"="blue","dEav"="orange", "dIav"="red", "dRav"="green" )) +
  theme_minimal() +
  theme(legend.title = element_blank())

library(tidyverse)
outdf_leaky %>%
  mutate(Na = Sa+Ea+Ia+Ra+Va+Eav+Iav+Rav) %>%
  mutate(Nb = Sb+Eb+Ib+Rb+Vb+Ebv+Ibv+Rbv) -> df2


#---plot B---
outdf_leaky <- as.data.frame(output_leaky)
ggplot(outdf_leaky, aes(x = time)) +
  geom_line(size = 0.7, aes(y = Sb, color = "dSb")) +
  geom_line(size = 0.7, aes(y = Eb, color = "dEb")) +
  geom_line(size = 0.7, aes(y = Ib, color = "dIb")) +
  geom_line(size = 0.7, aes(y = Rb, color = "dRb")) +
  geom_point(size = 0.6, aes(y = Vb, color = "dVb")) +
  geom_point(size = 0.6, aes(y = Ebv, color = "dEbv")) +
  geom_point(size = 0.6, aes(y = Ibv, color = "dIbv")) +
  geom_point(size = 0.6, aes(y = Rbv, color = "dRbv")) +
  labs(x = "Time(days)", y = "Proportion of population", title = "SEIR Dynamics for the population B") +
  scale_color_manual(values = c("dSb"="blue", "dEb"="orange", "dIb"= "red","dRb"="green", "dVb"="blue","dEbv"="orange", "dIbv"="red", "dRbv"="green" )) +
  theme_minimal() +
  theme(legend.title = element_blank())


# probabilidade de mutacao

prob_mut = .0001

# I = 1000

# Numero esperado de infectados com nova mutacao: I*prob_mut

I = 1000

N=I*prob_mut
N

I = 1000000

N=I*prob_mut
N


# 1
# considerar tamanhos de populcao distintos: N=10000, N=1 000 000, N=1 000 000 (ou outros valores)
# encontrar o que se espera em termos de individuos com as novas mutacoes

# 2
# variar a prob_mut e construir grafico do numero esperado em funcao da prob de mutacao


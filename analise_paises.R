library(tidyverse)

df <- read_csv("./vaccinations.csv")

df %>%
  filter(location == "Afghanistan" | location == "Brazil" | location == "United States") %>%
  ggplot(aes(x=date, y=people_fully_vaccinated_per_hundred)) + 
  geom_step() +
  facet_grid(location ~ .)

df %>%
  filter(people_fully_vaccinated_per_hundred <=100) %>%
  #filter(location == "Afghanistan" | location == "Brazil" | location == "United States") %>%
  ggplot(aes(x=date, y=location, fill = people_fully_vaccinated_per_hundred)) + 
  geom_raster() 

# possibilidades:
# 1 - preencher os espaços vazios com o número anterior
# 2 - acertar o percentual para maximo de 100% (ou entender a variavel)

df %>%
  filter(people_fully_vaccinated_per_hundred>100) -> xx

# Définir les paramètres
alpha <- 0.05
n <- 50
s <- 0.161 # Écart type de l'échantillon
mean <- 0.876 # Moyenne de l'échantillon

# Calculer la valeur critique de la t-distribution
t_critical <- qt(1 - alpha / 2, df = n - 1)

# Calculer l'intervalle de confiance
margin_of_error <- t_critical * (s / sqrt(n))
lower_bound <- mean - margin_of_error
upper_bound <- mean + margin_of_error

intervalle_de_confiance <- c(lower_bound, upper_bound)
print(intervalle_de_confiance)
print(margin_of_error)

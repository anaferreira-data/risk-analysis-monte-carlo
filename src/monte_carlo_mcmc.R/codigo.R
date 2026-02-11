GBrowniano <- function(T, n, seed, mu, sigma, x0) {
  # Define o tamanho do passo de tempo
  h <- T / n
  r <- mu - 0.5 * sigma^2
  sh <- sqrt(h)
  mh <- r * h
  
  # Inicializa o processo
  X <- x0
  P <- matrix(0, nrow = n+1, ncol = 2)
  P[1, ] <- c(0, x0)
  
  # Define semente para reprodutibilidade
  set.seed(seed)
  
  # Loop para gerar o movimento browniano geométrico
  for (j in 1:n) {
    z <- rnorm(1)  # valor aleatório padrão normal
    X <- X * exp(mh + sh * sigma * z)
    P[j+1, ] <- c(j * h, X)
  }
  
  # Retorna os resultados
  colnames(P) <- c("time", "X")
  return(P)
}

t <- GBrowniano(10, 2000, 159284, 0.1, 0.15, 100)
#convertemos para dataframe por facilidade
dft <- as.data.frame(t)
plot(dft$time, dft$X, type = "l", col = 'blue', lwd = 2,
     xlab = "t", ylab = "X(t)",
     main = "Realização do movimento Browniano Geométrico\n com taxa de crescimento esperado 0.1, volatilidade 0.15, preço inicial = 100",
     ylim=c(40,300),cex.main=0.8)
grid()

#Etapa 1 - simulação monte carlo
set.seed(123)

N <- 10000     # número de simulações
X0 <- 100      # preço inicial
sigma <- 0.15  # volatilidade
T <- 10        # tempo em anos
r <- 0.05      # taxa livre de risco
K <- 97        # preço de exercício da opção
# -------------------------------
# Simulação do preço futuro 
# -------------------------------
Z <- rnorm(N)  # ruído normal padrão
XT <- X0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)
# ------------------------------
# Payoff da opção de venda (put)
# ------------------------------
payoff_put <- pmax(K - XT, 0)
V_put <- exp(-r * T) * mean(payoff_put)  # valor estimado da put
se_put <- exp(-r * T) * sd(payoff_put) / sqrt(N)  # erro padrão

cat("Valor estimado da put:", round(V_put, 4), "\n")
## Valor estimado da put: 2.3804
cat("Erro padrão:", round(se_put, 4), "\n")
## Erro padrão: 0.0626

BS_put_MC <- function(m, t, te, sigma, r, k, x0) {
  
  m1 <- (r - 0.5 * sigma^2) * (te - t)
  m2 <- sigma * sqrt(te - t)
  disc <- exp(-r * (te - t))
  
  s1 <- 0
  s2 <- 0
  
  for (j in 1:m) {
    z <- rnorm(1)  
    
    # Variáveis antitéticas
    xa <- x0 * exp(m1 + m2 * z)
    xb <- x0 * exp(m1 - m2 * z)
    
    # Payoff da PUT (diferença principal em relação à CALL)
    payoffa <- max(k - xa, 0)
    payoffb <- max(k - xb, 0)
    
    payoff_avg <- (payoffa + payoffb) / 2
    
    s1 <- s1 + payoff_avg
    s2 <- s2 + payoff_avg^2
  }
  
  put <- disc * s1 / m
  std <- disc * sqrt((s2 - s1^2 / m) / (m - 1) / m)
  
  cat("Put price =", round(put,4), "\n")
  cat("Standard error =",round(std,4), "\n")
  cat("95% confidence interval =", round(put - 1.96*std,4), "to", round(put + 1.96*std,4), "\n")
}
BS_put_MC(m = 10000, t = 0, te = 10, sigma = 0.15, r = 0.05, k = 97, x0 = 100)
## Put price = 2.3797 
## Standard error = 0.0411 
## 95% confidence interval = 2.2991 to 2.4603

BS_put <- function(X, K, r, sigma, T) {
  d1 <- (log(X/K) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  put <- K*exp(-r*T)*pnorm(-d2) - X*pnorm(-d1)
  return(round(put,4))
}

X <- 100; K <- 97; r <- 0.05; sigma <- 0.15; T <- 10
BS_put(X, K, r, sigma, T)
## [1] 2.3634

#1.1 medidas de risco 
# Valor da carteira = ações + opção de venda
# ------------------------------
# O investidor tem 100 ações e 1 put
# ------------------------------
p_carteira <- 100 * XT + exp(-r * T) * payoff_put

# Retorno da carteira (logarítmico)
p_retornos <- log(p_carteira / (100 * X0))

alpha <- 0.95
p_VaR_95 <- quantile(p_retornos, probs = 1 - alpha)
p_ES_95 <- mean(p_retornos[p_retornos < p_VaR_95])

cat("VaR 95%:", round(p_VaR_95, 4), "\n")
## VaR 95%: -0.3886
cat("ES 95%:", round(p_ES_95, 4), "\n")
## ES 95%: -0.5862
# ------------------------------
#Carteira sem put (sem 'seguro')
# ------------------------------
carteira <- 100 * XT
retornos <- log(carteira / (100 * X0))
alpha <- 0.95

VaR_95 <- quantile(retornos, probs = 1 - alpha)
ES_95 <- mean(retornos[retornos < VaR_95])

cat("VaR 95%:", round(VaR_95, 4), "\n")
## VaR 95%: -0.3913
cat("ES 95%:", round(ES_95, 4), "\n")
## ES 95%: -0.5909
# ------------------------------
# Comparação
# ------------------------------
par(mfrow=c(1,2))
hist(p_retornos, breaks = 50, col = "lightblue", border = "White",
     main = "Distribuição dos Retornos da Carteira com put",
     xlab = "Retorno logarítmico",cex.main=0.8)
abline(v = p_VaR_95, col = "red", lwd = 2, lty = 2)
abline(v = p_ES_95, col = "darkred", lwd = 2, lty = 3)
legend("topright", legend = c("VaR 95%", "ES 95%"),
       col = c("red", "darkred"), lty = c(2, 3), lwd = 2, bty = "n")
hist(retornos, breaks = 50, col = "lightblue", border = "White",
     main = "Distribuição dos Retornos da Carteira sem put",
     xlab = "Retorno logarítmico",cex.main=0.8)
abline(v = VaR_95, col = "red", lwd = 2, lty = 2)
abline(v = ES_95, col = "darkred", lwd = 2, lty = 3)
legend("topright", legend = c("VaR 95%", "ES 95%"),
       col = c("red", "darkred"), lty = c(2, 3), lwd = 2, bty = "n")


#etapa.2 - monte carlo via cadeias de markov
set.seed(123)
# ----------------------------
# Simulação de dados verdadeiros
# ----------------------------
Nsim <- 10000
X0 <- 100
mu_true <- 0.1
sigma_true <- 0.15
T <- 10
r<-0.05

Z <- rnorm(Nsim)
XT <- X0 * exp((mu_true - 0.5 * sigma_true^2) * T + sigma_true * sqrt(T) * Z)
Y <- log(XT / X0)   # log-retorno TOTAL
# ----------------------------
# Funções de verossimilhança e prior (ajustadas para Y)
# ----------------------------
loglik_total <- function(mu, sigma, Y, T) {
  mean_theoretical <- (mu - 0.5 * sigma^2) * T
  sd_theoretical <- sigma * sqrt(T)
  sum(dnorm(Y, mean = mean_theoretical, sd = sd_theoretical, log = TRUE))
}

prior <- function(mu, sigma) {
  # priors leves: mu ~ N(0,1), sigma ~ Gamma(2,10) (em escala positiva)
  dnorm(mu, 0, 1, log = TRUE) + dgamma(sigma, shape = 2, rate = 10, log = TRUE)
}

# ----------------------------
# MCMC (Metropolis-Hastings)
# ----------------------------
Niter <- 10000
mu_chain <- numeric(Niter)
sigma_chain <- numeric(Niter)
mu_chain[1] <- 0.08
sigma_chain[1] <- 0.25

for (i in 2:Niter) {
  # proposta para mu e log(sigma)
  mu_prop <- rnorm(1, mu_chain[i-1], 0.01)
  log_sigma_prop <- rnorm(1, log(sigma_chain[i-1]), 0.05)
  sigma_prop <- exp(log_sigma_prop)
  
  # evitar sigma <= 0
  if (sigma_prop <= 0) {
    mu_chain[i] <- mu_chain[i-1]
    sigma_chain[i] <- sigma_chain[i-1]
    next
  }
  
  logA <- (loglik_total(mu_prop, sigma_prop, Y, T) + prior(mu_prop, sigma_prop)) -
    (loglik_total(mu_chain[i-1], sigma_chain[i-1], Y, T) + prior(mu_chain[i-1], sigma_chain[i-1]))
  
  if (log(runif(1)) < logA) {
    mu_chain[i] <- mu_prop
    sigma_chain[i] <- sigma_prop
  } else {
    mu_chain[i] <- mu_chain[i-1]
    sigma_chain[i] <- sigma_chain[i-1]
  }
}
#Análise da Convergência
cum_mean<-cumsum(mu_chain)/seq_along(mu_chain)
cum_var<-cumsum((mu_chain-cum_mean)^2)/pmax(seq_along(mu_chain)-1,1)
cum_error<-1.96*sqrt(cum_var/seq_along(mu_chain))
dados<-data.frame(
  iter=seq_along(mu_chain),
  mean=cum_mean,
  upper=cum_mean+cum_error,
  lower=cum_mean-cum_error
)

library(ggplot2)
## Warning: package 'ggplot2' was built under R version 4.5.2
pm<-ggplot(dados,aes(x=iter))+
  geom_line(aes(y=mean),col='blue',lwd=1)+
  geom_line(aes(y=upper),col='black')+
  geom_line(aes(y=lower),col='black')+
  geom_hline(yintercept=0.1)+
  labs(title = expression("Convergência da CM para " * mu),
       y=expression("média de " * mu_true))

cum_mean<-cumsum(sigma_chain)/seq_along(sigma_chain)
cum_var<-cumsum((sigma_chain-cum_mean)^2)/pmax(seq_along(sigma_chain)-1,1)
cum_error<-1.96*sqrt(cum_var/seq_along(sigma_chain))
dados<-data.frame(
  iter=seq_along(sigma_chain),
  mean=cum_mean,
  upper=cum_mean+cum_error,
  lower=cum_mean-cum_error
)

ps<-ggplot(dados,aes(x=iter))+
  geom_line(aes(y=mean),col='blue',lwd=1)+
  geom_line(aes(y=upper),col='black')+
  geom_line(aes(y=lower),col='black')+
  geom_hline(yintercept=0.15)+
  labs(title= expression("Convergência CM para " * sigma),
       y=expression("média de " * sigma_true))


library(gridExtra)
grid.arrange(pm, ps, ncol = 2) 

burn_in <- 2500
mu_post <- mu_chain[2501:Niter]
sigma_post <- sigma_chain[(burn_in+1):Niter]

cat("Posterior mean mu:", round(mean(mu_post),4), "\n")
## Posterior mean mu: 0.0999
cat("Posterior mean sigma:", round(mean(sigma_post),4), "\n")
## Posterior mean sigma: 0.1499

#trace plots e histogramas
par(mfrow = c(1, 2))
# ------------------------------
# Trace Plots
# ------------------------------
plot(2501:Niter,mu_post, type = "l", col = "steelblue", lwd = 1,
     main = expression("Trace Plot para " * mu), xlab = "Iteração", ylab = expression(mu))
plot(2501:Niter,sigma_post, type = "l", col = "darkorange", lwd = 1,
     main = expression("Trace Plot para " * sigma), xlab = "Iteração", ylab = expression(sigma))

# ------------------------------
# Distribuições posteriores
# ------------------------------
h1<-hist(mu_post,freq=F,breaks = 50, col = "steelblue", border = "white",
         main = expression("Distribuição Posterior de " * mu), xlab = expression(mu))
h2<-hist(sigma_post,freq=F, breaks = 50, col = "darkorange", border = "white",
         main = expression("Distribuição Posterior de " * sigma), xlab = expression(sigma))
# ------------------------------
# Confirmação as densidades
# ------------------------------
sum(h1$density * diff(h1$breaks))
## [1] 1
sum(h2$density * diff(h2$breaks))
## [1] 1

#Etapa 3- reavaliação das medidas de risco
# -------------------------------
# Cálculo das medidas de risco
# -------------------------------
set.seed(5)
N_scenarios <- length(mu_chain)  # cada amostra da posterior

# Inicializar vetor de retornos da carteira
retornos_post <- numeric(N_scenarios)

for (i in 1:N_scenarios) {
  # Simula XT da carteira com cada par (mu, sigma) da posterior
  XT_sim <- X0 * exp((r - 0.5*sigma_chain[i]^2)*T + sigma_chain[i]*sqrt(T)*rnorm(1))
  payoff_post<-pmax(K-XT_sim,0)
  # Valor da carteira: 100 ações + 1 put
  carteira_post <- 100 * XT_sim + exp(-r*T)*payoff_post
  
  # Retorno logarítmico da carteira
  retornos_post[i] <- log(carteira_post / (100*X0))
}

# VaR e ES 95%
par(mfrow=c(1,1))
alpha <- 0.95
VaR_95_post <- quantile(retornos_post, probs = 1 - alpha)
ES_95_post <- mean(retornos_post[retornos_post < VaR_95_post])

cat("VaR 95% da carteira:", round(VaR_95_post,4), "\n")
## VaR 95% da carteira: -0.3919
cat("ES 95% da carteira:", round(ES_95_post,4), "\n")
## ES 95% da carteira: -0.6046
# ------------------------------
# Comparação
# ------------------------------
df <- data.frame(
  retorno = c(p_retornos, retornos_post),
  type = factor(rep(c("Monte Carlo", "Posterior MCMC"), times = c(N, N_scenarios)))
)

library(ggplot2)
linhas <- data.frame(
  x = c(p_VaR_95, p_ES_95),
  tipo = c("VaR", "ES")
)
ggplot(df, aes(x = retorno, fill = type)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = p_VaR_95, color = "red", linetype = "dashed") +
  geom_vline(xintercept = p_ES_95, color = "darkred", linetype = "dashed") +
  geom_vline(xintercept = VaR_95_post, color = "red") +
  geom_vline(xintercept = ES_95_post, color = "darkred") +
  geom_vline(data = linhas, aes(xintercept = x, color = tipo), linetype = "dashed",lwd=0) +
  scale_color_manual(
    name = "Linhas de risco",
    values = c("VaR" = "red", "ES" = "darkred")
  ) +
  labs(
    title = "Comparativo VaR/ES: Monte Carlo vs Posterior MCMC",
    x = "Retorno logarítmico",
    y = "Densidade",
    fill = "Tipo de simulação"
  ) +
  
  theme_minimal()


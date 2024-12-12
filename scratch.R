# repeated sims # N # p

hist(rbinom(10000, 100, 0.99), xlim = c(0, 100), 
     ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.975), xlim = c(0, 100),
     ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.95), xlim = c(0, 100),
  ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.90), xlim = c(0, 100),
     ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.75), xlim = c(0, 100),
     ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.50), xlim = c(0, 100),
     ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.25), xlim = c(0, 100),
     ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.1), xlim = c(0, 100),
     ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.05), xlim = c(0, 100),
     ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.025), xlim = c(0, 100),
     ylim = c(0, 4000), breaks=101)

hist(rbinom(10000, 100, 0.01), xlim = c(0, 100), 
     ylim = c(0, 4000), breaks=101)


# for box 1
#par(mar=c(b,l,t,r))
par(mar=c(4,5,0,0))

n_sims <- 1000
N <- 30

hist(rbinom(n_sims, N, 0.90), 
     xlim = c(0, 35), ylim = c(0, 350),
     cex.lab=1.75, cex.axis=1.75,
     xlab = bquote(y[i]), ylab = "frequency",
     main = "",
     breaks = seq(from=0, to=40, by=1))
abline(v = N, col = "skyblue4", lty = "dashed", lwd = 4)

hist(rbinom(n_sims, N, 0.50), 
     xlim = c(0, 35), ylim = c(0, 350),  
     cex.lab=1.75, cex.axis=1.75,
     xlab = bquote(y[i]), ylab = "frequency",
     main = "",
     breaks = seq(from=0, to=40, by=1))
abline(v = N, col = "skyblue4", lty = "dashed", lwd = 4)

hist(rbinom(n_sims, N, 0.10), 
     xlim = c(0, 35), ylim = c(0, 350),  
     cex.lab=1.75, cex.axis=1.75,
     xlab = bquote(y[i]), ylab = "frequency",
     main = "",
     breaks = seq(from=0, to=40, by=1))
abline(v = N, col = "skyblue4", lty = "dashed", lwd = 4)

dev.off()


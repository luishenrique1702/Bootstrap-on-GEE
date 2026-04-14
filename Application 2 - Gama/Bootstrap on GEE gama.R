  library(dplyr)
	library(glmtoolbox)
	library(gamlss)
	set.seed(1234)
	
	dados <- as_tibble(Milk)
	cluster <- unique(dados$Cow)
	R <- 1000 
	
	n<- length(cluster)
	model <- glmgee(protein ~ Diet + Time + I(Time^2), 
	data = Milk, 
	corstr = "AR-M-dependent(2)", 
	family = Gamma("log"), 
	id = Cow)
	C <- ncol(model.matrix(model))
	observed_coef <- coef(model)
	observed_phi <- model$phi
	observed_rho1 <- model$corr[1, 2]
	observed_rho2 <- model$corr[1, 3]
	
	boot_coefs <- matrix(NA, nrow = R, ncol = C)
	boot_phi <- numeric(R)
	boot_rho1 <- numeric(R)
	boot_rho2 <- numeric(R)
	
	for (i in 1:R) {
		amostra_clusters <- sample(cluster, size = n, replace = TRUE)
		dados_boot <- lapply(1:n, function(j) {
			dados %>% filter(Cow == amostra_clusters[j]) %>% mutate(new_id = j)
		}) %>% bind_rows()
		
		mod_boot <- glmgee(protein ~ Diet + Time + I(Time^2), 
		data = dados_boot,
		corstr = "AR-M-dependent(2)", 
		family = Gamma("log"), 
		id = new_id)
		
		boot_coefs[i, ] <- coef(mod_boot)
		boot_phi[i] <- mod_boot$phi
		boot_rho1[i] <- mod_boot$corr[1, 2]
		boot_rho2[i] <- mod_boot$corr[1, 3]
	}
	
	
	results_df <- as.data.frame(boot_coefs)
	names(results_df) <- names(observed_coef)
	
	jack_coefs <- matrix(NA, nrow = n, ncol = C)
	jack_phi <- numeric(n)
	jack_rho1 <- numeric(n)
	jack_rho2 <- numeric(n)
	
	for (i in 1:n) {
		dados_jack <- dados %>% filter(Cow != cluster[i])
		mod_jack <- glmgee(protein ~ Diet + Time + I(Time^2), 
		data = dados_jack,
		corstr = "AR-M-dependent(2)", 
		family = Gamma("log"), 
		id = Cow)
		jack_coefs[i, ] <- coef(mod_jack)
		jack_phi[i] <- mod_jack$phi
		jack_rho1[i] <- mod_jack$corr[1, 2]
		jack_rho2[i] <- mod_jack$corr[1, 3]
	}
	
	BCa_interval <- function(boot_data, original_value, jack_data, conf = 0.95) {
		z0 <- qnorm(sum(boot_data <= original_value, na.rm = TRUE) / (length(boot_data) + 1))
		
		jack_mean <- mean(jack_data, na.rm = TRUE)
		num <- sum((jack_mean - jack_data)^3, na.rm = TRUE)
		den <- 6 * (sum((jack_data - jack_mean)^2, na.rm = TRUE))^(3/2)
		a <- ifelse(den == 0, 0, num / den)
		
		alpha <- 1 - conf
		z_critical <- qnorm(1 - alpha / 2)
		p1 <- pnorm(z0 + (z0 - z_critical) / (1 - a * (z0 - z_critical)))
		p2 <- pnorm(z0 + (z0 + z_critical) / (1 - a * (z0 + z_critical)))
		
		interval <- quantile(boot_data, probs = c(p1, p2), type = 6, na.rm = TRUE)
		return(interval)
	}
	
	
	pdf("bootstrap_plots.pdf")
	
	par(mfrow = c(3, 3), mar = c(5, 4, 2, 1))
	
	
	for (i in 1:ncol(results_df)) {
		coef_data <- results_df[, i]
		mean_coef <- mean(coef_data, na.rm = TRUE)
		sd_coef <- sd(coef_data, na.rm = TRUE)
		
		hist_data <- hist(coef_data, plot = FALSE)
		dens_kernel <- density(coef_data, na.rm = TRUE)
		max_y <- max(max(hist_data$density), max(dens_kernel$y)) * 1.4
		
		hist(coef_data, freq = FALSE, col = "gray", border = "black", main = "",
		xlab = c(expression(hat(beta)[0]~"(Intercept)"), expression(hat(beta)[1]~"(Barley+Lupins)"), 
		expression(hat(beta)[2]~"(Lupins)"), expression(hat(beta)[3]~"(Time)"), expression(hat(beta)[4]~(Time^2)))[i], 
		ylab = "Density", ylim = c(0, max_y))
		
		lines(dens_kernel$x, dens_kernel$y, col = "blue", lwd = 2)
		curve(dnorm(x, mean = mean_coef, sd = sd_coef), add = TRUE, col = "magenta", lwd = 2, lty = 2)
		
		abline(v = observed_coef[i], col = "black", lwd = 2, lty = 2)
		
		
		abline(v = median(coef_data, na.rm = TRUE), col = "green", lwd = 2, lty = 1)
		
		bca_ci <- BCa_interval(coef_data, observed_coef[i], jack_coefs[, i])
		segments(x0 = bca_ci[1], x1 = bca_ci[2], y0 = max_y*0.05, y1 = max_y*0.05, col = "black", lwd = 3)
		box()
	}
	
	mean_phi <- mean(boot_phi, na.rm = TRUE)
	sd_phi <- sd(boot_phi, na.rm = TRUE)
	hist_phi <- hist(boot_phi, plot = FALSE)
	dens_phi <- density(boot_phi, na.rm = TRUE)
	max_y_phi <- max(max(hist_phi$density), max(dens_phi$y)) * 1.4
	hist(boot_phi, freq = FALSE, col = "gray", border = "black", main = "",
	xlab = expression(hat(phi)), ylab = "Density", ylim = c(0, max_y_phi))
	lines(dens_phi$x, dens_phi$y, col = "blue", lwd = 2)
	curve(dnorm(x, mean = mean_phi, sd = sd_phi), add = TRUE, col = "magenta", lwd = 2, lty = 2)
	abline(v = observed_phi, col = "black", lwd = 2, lty = 2)
	
	
	abline(v = median(boot_phi, na.rm = TRUE), col = "green", lwd = 2, lty = 1)
	
	bca_ci_phi <- BCa_interval(boot_phi, observed_phi, jack_phi)
	segments(x0 = bca_ci_phi[1], x1 = bca_ci_phi[2], y0 = max_y_phi*0.05, y1 = max_y_phi*0.05, col = "black", lwd = 3)
	box()
	
	mean_rho1 <- mean(boot_rho1, na.rm = TRUE)
	sd_rho1 <- sd(boot_rho1, na.rm = TRUE)
	hist_rho1 <- hist(boot_rho1, plot = FALSE)
	dens_rho1 <- density(boot_rho1, na.rm = TRUE)
	max_y_rho1 <- max(max(hist_rho1$density), max(dens_rho1$y)) * 1.4
	hist(boot_rho1, freq = FALSE, col = "gray", border = "black", main = "",
	xlab = expression(hat(rho)[1]), ylab = "Density", ylim = c(0, max_y_rho1))
	lines(dens_rho1$x, dens_rho1$y, col = "blue", lwd = 2)
	curve(dnorm(x, mean = mean_rho1, sd = sd_rho1), add = TRUE, col = "magenta", lwd = 2, lty = 2)
	abline(v = observed_rho1, col = "black", lwd = 2, lty = 2)
	
	
	abline(v = median(boot_rho1, na.rm = TRUE), col = "green", lwd = 2, lty = 1)
	
	bca_ci_rho1 <- BCa_interval(boot_rho1, observed_rho1, jack_rho1)
	segments(x0 = bca_ci_rho1[1], x1 = bca_ci_rho1[2], y0 = max_y_rho1*0.05, y1 = max_y_rho1*0.05, col = "black", lwd = 3)
	box()
	
	mean_rho2 <- mean(boot_rho2, na.rm = TRUE)
	sd_rho2 <- sd(boot_rho2, na.rm = TRUE)
	hist_rho2 <- hist(boot_rho2, plot = FALSE)
	dens_rho2 <- density(boot_rho2, na.rm = TRUE)
	max_y_rho2 <- max(max(hist_rho2$density), max(dens_rho2$y)) * 1.4
	hist(boot_rho2, freq = FALSE, col = "gray", border = "black", main = "",
	xlab = expression(hat(rho)[2]), ylab = "Density", ylim = c(0, max_y_rho2))
	lines(dens_rho2$x, dens_rho2$y, col = "blue", lwd = 2)
	curve(dnorm(x, mean = mean_rho2, sd = sd_rho2), add = TRUE, col = "magenta", lwd = 2, lty = 2)
	abline(v = observed_rho2, col = "black", lwd = 2, lty = 2)
	
	# --- LINHA ADICIONADA ---
	abline(v = median(boot_rho2, na.rm = TRUE), col = "green", lwd = 2, lty = 1)
	
	bca_ci_rho2 <- BCa_interval(boot_rho2, observed_rho2, jack_rho2)
	segments(x0 = bca_ci_rho2[1], x1 = bca_ci_rho2[2], y0 = max_y_rho2*0.05, y1 = max_y_rho2*0.05, col = "black", lwd = 3)
	box()
	
	plot.new()
	
	legend("center", legend = c("Kernel Density", "Adjusted Normal", "Observed Value",  "Bootstrap Median", "CI 95% BCa"),
	col = c("blue", "magenta", "black", "green", "black"), 
	lty = c(1, 2, 2, 1, 1), 
	lwd = c(2, 2, 2, 2, 4), bty = "n")
	
	
	# Fecha o dispositivo PDF
	dev.off()
	
	
	bias_coef <- observed_coef - apply(boot_coefs, 2, mean, na.rm = TRUE)
	se_coef <- apply(boot_coefs, 2, sd, na.rm = TRUE)
	
	bias_phi <- observed_phi - mean(boot_phi, na.rm = TRUE)
	se_phi <- sd(boot_phi, na.rm = TRUE)
	
	bias_rho1 <- observed_rho1 - mean(boot_rho1, na.rm = TRUE)
	se_rho1 <- sd(boot_rho1, na.rm = TRUE)
	
	bias_rho2 <- observed_rho2 - mean(boot_rho2, na.rm = TRUE)
	se_rho2 <- sd(boot_rho2, na.rm = TRUE)
	
	nomes_dos_parametros <- c(
	"Intercepto",
	"Diet(Lupins)",
	"Diet(Barley+Lupins)",
	"Time",
	"Time^2",
	"phi",
	"rho1",
	"rho2"
	)
	
	coef_table <- data.frame(
	Parametro = nomes_dos_parametros,
	
	Observado = c(observed_coef, observed_phi, observed_rho1, observed_rho2),
	
	Media_Bootstrap = c(apply(boot_coefs, 2, mean, na.rm = TRUE), 
	mean(boot_phi, na.rm = TRUE), 
	mean(boot_rho1, na.rm = TRUE), 
	mean(boot_rho2, na.rm = TRUE)),
	
	Vies_Bootstrap = c(bias_coef, bias_phi, bias_rho1, bias_rho2),
	
	Erro_Padrao_Bootstrap = c(se_coef, se_phi, se_rho1, se_rho2)
	)
	
	coef_table[, -1] <- round(coef_table[, -1], 5)
	
	print(coef_table)
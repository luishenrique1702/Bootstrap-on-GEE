	library(dplyr)
	library(glmtoolbox)
	# install.packages("aplore3")
	library(aplore3)
	library(tidyverse)
	
	set.seed(1234)
	
	
	poly <- polypharm
	poly$resp <- as.numeric(poly$polypharmacy) - 1
	poly$tage <- log(poly$age / 100)
	dados <- as_tibble(poly)
	
	cluster <- unique(dados$id)
	n <- length(cluster)
	R <- 1000 
	model <- glmgee(resp ~ mhv4 + gender + tage, id = id, data = dados,
	corstr = "AR-M-dependent(2)", family = binomial)
	
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
		if (i %% 100 == 0) 
		
		
		amostra_clusters <- sample(cluster, size = n, replace = TRUE)
		
		
		dados_boot <- lapply(1:n, function(j) {
			dados %>% filter(id == amostra_clusters[j]) %>% mutate(new_id = j)
		}) %>% bind_rows()
		
		
		try({
			mod_boot <- glmgee(resp ~ mhv4 + gender + tage, data = dados_boot,
			corstr = "AR-M-dependent(2)", family = binomial, id = new_id)
			
			boot_coefs[i, ] <- coef(mod_boot)
			boot_phi[i] <- mod_boot$phi
			boot_rho1[i] <- mod_boot$corr[1, 2]
			boot_rho2[i] <- mod_boot$corr[1, 3]
		}, silent = TRUE)
	}
	
	
	results_df <- as.data.frame(boot_coefs)
	names(results_df) <- names(observed_coef)
	
	
	jack_coefs <- matrix(NA, nrow = n, ncol = C)
	jack_phi <- numeric(n)
	jack_rho1 <- numeric(n)
	jack_rho2 <- numeric(n)
	
	
	for (i in 1:n) {
		if (i %% 50 == 0) 
		
		
		dados_jack <- dados %>% filter(id != cluster[i])
		
		try({
			mod_jack <- glmgee(resp ~ mhv4 + gender + tage, data = dados_jack,
			corstr = "AR-M-dependent(2)", family = binomial, id = id)
			
			jack_coefs[i, ] <- coef(mod_jack)
			jack_phi[i] <- mod_jack$phi
			jack_rho1[i] <- mod_jack$corr[1, 2]
			jack_rho2[i] <- mod_jack$corr[1, 3]
		}, silent = TRUE)
	}
	
	
	BCa_interval <- function(boot_data, original_value, jack_data, conf = 0.95) {
		
		boot_data <- boot_data[!is.na(boot_data)]
		jack_data <- jack_data[!is.na(jack_data)]
		
		if (length(boot_data) < 10 || length(jack_data) < 3) return(c(NA, NA))
		
		
		z0 <- qnorm(sum(boot_data <= original_value) / (length(boot_data) + 1))
		
		
		jack_mean <- mean(jack_data)
		num <- sum((jack_mean - jack_data)^3)
		den <- 6 * (sum((jack_data - jack_mean)^2))^(3/2)
		a <- ifelse(den == 0, 0, num / den)
		
		
		alpha <- 1 - conf
		z_critical <- qnorm(1 - alpha / 2)
		p1 <- pnorm(z0 + (z0 - z_critical) / (1 - a * (z0 - z_critical)))
		p2 <- pnorm(z0 + (z0 + z_critical) / (1 - a * (z0 + z_critical)))
		
		
		p1 <- max(0, min(1, p1))
		p2 <- max(0, min(1, p2))
		
		interval <- quantile(boot_data, probs = c(p1, p2), type = 6, na.rm = TRUE)
		return(interval)
	}
	
	
	pdf("bootstrap_plots_polypharm.pdf")
	
	
	layout(matrix(c(1, 2, 3, 
	4, 5, 6, 
	7, 8, 9, 
	10, 10, 10), nrow = 4, ncol = 3, byrow = TRUE),
	heights = c(3, 3, 3, 1))
	
	par(mar = c(4, 4, 2, 1) + 0.1)
	
	
	xlabels_coef <- c(expression(hat(beta)[0]~"(Intercept)"), 
	expression(hat(beta)[1]~"(1 to 5)"), 
	expression(hat(beta)[2]~"(6 to 14)"), 
	expression(hat(beta)[3]~"(More than 14)"), 
	expression(hat(beta)[4]~"(Male)"),
	expression(hat(beta)[5]~"(log(age/100))"))
	
	
	for (i in 1:ncol(results_df)) {
		coef_data <- results_df[, i]
		obs_val <- observed_coef[i]
		jack_data <- jack_coefs[, i]
		
		mean_coef <- mean(coef_data, na.rm = TRUE)
		sd_coef <- sd(coef_data, na.rm = TRUE)
		
		hist_data <- hist(coef_data, plot = FALSE)
		dens_kernel <- density(coef_data, na.rm = TRUE)
		max_y <- max(max(hist_data$density), max(dens_kernel$y)) * 1.4
		
		hist(coef_data, freq = FALSE, col = "gray", border = "black", main = "",
		xlab = xlabels_coef[i], ylab = "Density", ylim = c(0, max_y))
		
		lines(dens_kernel$x, dens_kernel$y, col = "blue", lwd = 2)
		curve(dnorm(x, mean = mean_coef, sd = sd_coef), add = TRUE, col = "magenta", lwd = 2, lty = 2)
		abline(v = obs_val, col = "black", lwd = 2, lty = 2)
		abline(v = median(coef_data, na.rm = TRUE), col = "green", lwd = 2, lty = 1)
		
		# Chama BCa com os dados corretos
		bca_ci <- BCa_interval(coef_data, obs_val, jack_data)
		if (!any(is.na(bca_ci))) {
			segments(x0 = bca_ci[1], x1 = bca_ci[2], y0 = max_y * 0.05, y1 = max_y * 0.05, col = "black", lwd = 3)
		}
		box()
	}
	
	
	mean_phi <- mean(boot_phi, na.rm = TRUE)
	sd_phi <- sd(boot_phi, na.rm = TRUE)
	dens_phi <- density(boot_phi, na.rm = TRUE)
	max_y_phi <- max(dens_phi$y) * 1.4
	hist(boot_phi, freq = FALSE, col = "gray", border = "black", main = "",
	xlab = expression(hat(phi)), ylab = "Density", ylim = c(0, max_y_phi))
	lines(dens_phi, col = "blue", lwd = 2)
	curve(dnorm(x, mean = mean_phi, sd = sd_phi), add = TRUE, col = "magenta", lwd = 2, lty = 2)
	abline(v = observed_phi, col = "black", lwd = 2, lty = 2)
	abline(v = median(boot_phi, na.rm = TRUE), col = "green", lwd = 2, lty = 1)
	bca_phi <- BCa_interval(boot_phi, observed_phi, jack_phi)
	if (!any(is.na(bca_phi))) segments(bca_phi[1], max_y_phi*0.05, bca_phi[2], max_y_phi*0.05, col="black", lwd=3)
	box()
	
	
	mean_rho1 <- mean(boot_rho1, na.rm = TRUE)
	sd_rho1 <- sd(boot_rho1, na.rm = TRUE)
	dens_rho1 <- density(boot_rho1, na.rm = TRUE)
	max_y_rho1 <- max(dens_rho1$y) * 1.4
	hist(boot_rho1, freq = FALSE, col = "gray", border = "black", main = "",
	xlab = expression(hat(rho)[1]), ylab = "Density", ylim = c(0, max_y_rho1))
	lines(dens_rho1, col = "blue", lwd = 2)
	curve(dnorm(x, mean = mean_rho1, sd = sd_rho1), add = TRUE, col = "magenta", lwd = 2, lty = 2)
	abline(v = observed_rho1, col = "black", lwd = 2, lty = 2)
	abline(v = median(boot_rho1, na.rm = TRUE), col = "green", lwd = 2, lty = 1)
	bca_rho1 <- BCa_interval(boot_rho1, observed_rho1, jack_rho1)
	if (!any(is.na(bca_rho1))) segments(bca_rho1[1], max_y_rho1*0.05, bca_rho1[2], max_y_rho1*0.05, col="black", lwd=3)
	box()
	
	
	mean_rho2 <- mean(boot_rho2, na.rm = TRUE)
	sd_rho2 <- sd(boot_rho2, na.rm = TRUE)
	dens_rho2 <- density(boot_rho2, na.rm = TRUE)
	max_y_rho2 <- max(dens_rho2$y) * 1.4
	hist(boot_rho2, freq = FALSE, col = "gray", border = "black", main = "",
	xlab = expression(hat(rho)[2]), ylab = "Density", ylim = c(0, max_y_rho2))
	lines(dens_rho2, col = "blue", lwd = 2)
	curve(dnorm(x, mean = mean_rho2, sd = sd_rho2), add = TRUE, col = "magenta", lwd = 2, lty = 2)
	abline(v = observed_rho2, col = "black", lwd = 2, lty = 2)
	abline(v = median(boot_rho2, na.rm = TRUE), col = "green", lwd = 2, lty = 1)
	bca_rho2 <- BCa_interval(boot_rho2, observed_rho2, jack_rho2)
	if (!any(is.na(bca_rho2))) segments(bca_rho2[1], max_y_rho2*0.05, bca_rho2[2], max_y_rho2*0.05, col="black", lwd=3)
	box()
	
	
	par(mar = c(0, 0, 0, 0))
	plot.new()
	legend("left", legend = c("Kernel Density", "Adjusted Normal", "Observed Value", "Bootstrap Median", "CI 95% BCa"),
	col = c("blue", "magenta", "black", "green", "black"),
	lty = c(1, 2, 2, 1, 1), lwd = c(2, 2, 2, 2, 4), bty = "n", cex = 1)
	
	dev.off()
	
	
	bias_coef <- observed_coef - apply(boot_coefs, 2, mean, na.rm = TRUE)
	se_coef <- apply(boot_coefs, 2, sd, na.rm = TRUE)
	
	
	param_names <- c("Intercept", "MHV4(1-5)", "MHV4(6-14)", "MHV4(>14)", "Male", "log(age/100)", "phi", "rho1", "rho2")
	
	coef_table <- data.frame(
	Parametro = param_names,
	Observado = c(observed_coef, observed_phi, observed_rho1, observed_rho2),
	Media_Boot = c(colMeans(boot_coefs, na.rm=T), mean(boot_phi, na.rm=T), mean(boot_rho1, na.rm=T), mean(boot_rho2, na.rm=T)),
	Vies = c(bias_coef, observed_phi - mean(boot_phi, na.rm=T), observed_rho1 - mean(boot_rho1, na.rm=T), observed_rho2 - mean(boot_rho2, na.rm=T)),
	SE = c(se_coef, sd(boot_phi, na.rm=T), sd(boot_rho1, na.rm=T), sd(boot_rho2, na.rm=T))
	)
	
	coef_table[, -1] <- round(coef_table[, -1], 5)
	print(coef_table)
	summary(model)
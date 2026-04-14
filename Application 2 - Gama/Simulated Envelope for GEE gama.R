library(glmtoolbox)
	library(mvnfast)
	library(ggplot2)
	library(gamlss)
	
	
	model <- glmgee(protein ~ Diet + Time + I(Time^2), 
	data = Milk, 
	corstr = "AR-M-dependent(2)", 
	id = Cow, 
	family = Gamma("log"))
	
	
	mu_hat <- model$fitted.values
	y_obs <- model$y
	id <- model$id
	n <- length(y_obs)
	phi <- model$phi
	cluster_ids <- unique(id)
	cluster_indices <- split(seq_along(id), id)
	
	
	set.seed(123)
	r_quant <- numeric(n)
	for (i in 1:n) {
		shape <- 1 / phi
		scale <- mu_hat[i] * phi
		u <- pgamma(y_obs[i], shape = shape, scale = scale)
		r_quant[i] <- qnorm(u)
	}
	
	
	B <- 300
	sim_r_quant <- matrix(NA, nrow = B, ncol = n)
	set.seed(123)
	
	for (b in 1:B) {
		y_sim <- numeric(n)
		
		for (i in seq_along(cluster_indices)) {
			idx <- cluster_indices[[i]]
			n_i <- length(idx)
			R_i <- model$corr[1:n_i, 1:n_i]  
			
			Z <- mvnfast::rmvn(1, mu = rep(0, n_i), sigma = R_i)
			U <- pnorm(Z)
			
			shape <- 1 / phi
			scale_i <- mu_hat[idx] * phi
			y_i <- qgamma(U, shape = shape, scale = scale_i)
			
			y_sim[idx] <- y_i
		}
		
		r_b <- numeric(n)
		for (i in 1:n) {
			shape <- 1 / phi
			scale <- mu_hat[i] * phi
			u <- pgamma(y_sim[i], shape = shape, scale = scale)
			r_b[i] <- qnorm(u)
		}
		
		sim_r_quant[b, ] <- sort(r_b)
	}
	
	
	lower <- apply(sim_r_quant, 2, quantile, probs = 0.025, na.rm = TRUE)
	upper <- apply(sim_r_quant, 2, quantile, probs = 0.975, na.rm = TRUE)
	med   <- apply(sim_r_quant, 2, quantile, probs = 0.5, na.rm = TRUE)
	
	
	qq_data <- data.frame(
	theoretical = sort(qnorm(ppoints(n))),
	sample = sort(r_quant),
	lower = lower,
	upper = upper,
	med = med
	)
	
	qq_data$color <- ifelse(qq_data$sample >= qq_data$lower & qq_data$sample <= qq_data$upper,
	"Dentro", "Fora")
	
	pdf("envelope_gama.pdf")
	ggplot(qq_data, aes(x = theoretical, y = sample)) +
	geom_point(aes(color = color), size = 0.8) +
	geom_line(aes(y = med), color = "black", linetype = "dashed") +
	geom_line(aes(y = lower), color = "black") +
	geom_line(aes(y = upper), color = "black") +
	scale_color_manual(values = c("Dentro" = "black", "Fora" = "red")) +
	guides(color = "none") +
	labs(x = "Theoretical Quantile N(0,1)",
	y = "Quantile Residuals") +
	theme_test()
	dev.off()
	
	pdf("qqnorm_gama.pdf")
	qqnorm(r_quant,main = "")
	qqline(r_quant, col = "red")
	dev.off()
	fake_model <- list(
	residuals = r_quant,
	na.action = NULL  
	)
	
	class(fake_model) <- "fake_residuals"
	
	
	resid.fake_residuals <- function(object, ...) object$residuals
	
	pdf("wp_gama.pdf")
	wp(fake_model,ylim.all = .4)
	dev.off()
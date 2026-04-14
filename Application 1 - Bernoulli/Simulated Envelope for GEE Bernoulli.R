library(glmtoolbox)
	library(ggplot2)
	library(aplore3)
	library(mvnfast) 
	
	Poly <- polypharm
	Poly$tage <- log(Poly$age / 100)
	Poly$resp <- as.numeric(Poly$polypharmacy) - 1
	
	
	modelBer <- glmgee(resp ~ mhv4 + gender + tage, id = id,
	family = binomial, corstr = "AR-M-dependent(2)", data = Poly)
	
	
	mu_hat <- modelBer$fitted.values
	y_obs <- modelBer$y
	n <- length(y_obs)
	
	id <- modelBer$id
	cluster_ids <- unique(id)
	cluster_indices <- split(seq_along(id), id)
	
	set.seed(123)
	r_quant <- numeric(n)
	for (i in 1:n) {
		if (y_obs[i] == 1) {
			a <- pbinom(0, size = 1, prob = mu_hat[i])
			b <- 1
		} else {
			a <- 0
			b <- pbinom(0, size = 1, prob = mu_hat[i])
		}
		
		u <- runif(1, min = a, max = b)
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
			R_i <- modelBer$corr[1:n_i, 1:n_i]
			
			Z <- mvnfast::rmvn(1, mu = rep(0, n_i), sigma = R_i)
			
			U <- pnorm(Z)
			mu_i <- mu_hat[idx]
			
			
			y_i <- qbinom(U, size = 1, prob = mu_i)
			
			y_sim[idx] <- y_i 
		}
		
		r_b <- numeric(n)
		
		for (i in 1:n) {
			if (y_sim[i] == 1) {
				a <- pbinom(0, size = 1, prob = mu_hat[i])
				b_val <- 1
			} else {
				a <- 0
				b_val <- pbinom(0, size = 1, prob = mu_hat[i])
			}
			
			u <- runif(1, min = a, max = b_val)
			r_b[i] <- qnorm(u)
		}
		
		sim_r_quant[b, ] <- sort(r_b)
	}
	
	
	lower <- apply(sim_r_quant, 2, quantile, probs = 0.025, na.rm = TRUE)
	upper <- apply(sim_r_quant, 2, quantile, probs = 0.975, na.rm = TRUE)
	med <- apply(sim_r_quant, 2, quantile, probs = 0.5, na.rm = TRUE)
	
	
	qq_data <- data.frame(
	theoretical = sort(qnorm(ppoints(n))),
	sample = sort(r_quant),
	lower = lower,
	upper = upper,
	med = med
	)
	
	
	qq_data$color <- ifelse(qq_data$sample >= qq_data$lower & qq_data$sample <= qq_data$upper,
	"Dentro", "Fora")
	
	
	ggplot(qq_data, aes(x = theoretical, y = sample)) +
	geom_point(aes(color = color),size=0.8) +
	geom_line(aes(y = med), color = "black", linetype = "dashed") +
	geom_line(aes(y = lower), color = "black") +
	geom_line(aes(y = upper), color = "black") +
	scale_color_manual(values = c("Dentro" = "black", "Fora" = "red")) +
	guides(color = "none") + # REMOVE A LEGENDA
	labs(x = "Theoretical quantile N(0,1)",
	y = "Resíduos quantílicos aleatorizados") +theme_classic()
	
	pdf("envelope_Bernoulli.pdf")
	ggplot(qq_data, aes(x = theoretical, y = sample)) +
	geom_point(aes(color = color), size = 0.8) +
	geom_line(aes(y = med), color = "black", linetype = "dashed") +
	geom_line(aes(y = lower), color = "black") +
	geom_line(aes(y = upper), color = "black") +
	scale_color_manual(values = c("Dentro" = "black", "Fora" = "red")) +
	guides(color = "none") + 
	labs(x = "Theoretical Quantile N(0,1)",  # <-- Eixo X (Já estava em inglês)
	y = "Randomized Quantile Residuals") + # <-- Eixo Y (Traduzido)
	theme_test()
	dev.off()
	
	
	fake_model <- list(
	residuals = r_quant,
	na.action = NULL
	)
	pdf("qqnorm_bernoulli.pdf")
	qqnorm(r_quant, 
	cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, 
	main = "", # <-- Título em Inglês
	xlab = "Theoretical Quantiles",                     
	ylab = "Randomized Quantile Residuals")              
	
	qqline(r_quant, col = "red", lwd = 1.5)
	dev.off()
	class(fake_model) <- "fake_residuals"
	
	
	library(gamlss)
	wp(fake_model,ylim.all = .5)
	
	
	plot_worm <- function(r) {
		fake_model <- list(residuals = r, na.action = NULL)
		class(fake_model) <- "fake_residuals"
		resid.fake_residuals <- function(object, ...) object$residuals
		assign("resid.fake_residuals", resid.fake_residuals, envir = .GlobalEnv)
		
	}
	
	plot_worm(r_quant)
	
	custom_rqres_plot <- function(residuals, howmany = 8, ylim_all = 0.6) {
		library(gamlss)
		
		
		n <- length(residuals)
		theoretical <- qnorm(ppoints(n))
		breaks <- qnorm(seq(0, 1, length.out = howmany + 1))
		
		
		df <- data.frame(
		residuals = residuals,
		theoretical = theoretical,
		group = cut(theoretical, breaks = breaks, include.lowest = TRUE)
		)
		
		
		op <- par(mfrow = c(ceiling(howmany / 2), 2), mar = c(4, 4, 2, 1))
		on.exit(par(op))
		
		levels_group <- levels(df$group)
		
		for (g in levels_group) {
			res_sub <- df$residuals[df$group == g]
			if (length(res_sub) >= 5) {
				wp(list(residuals = res_sub, na.action = NULL),
				main = paste("Grupo:", g),
				ylim.all = ylim_all)
			} else {
				plot.new()
				text(0.5, 0.5, paste("Group", g, "with few observations"))
			}
		}
	}
	
	pdf("wp_bernoulli.pdf")
	custom_rqres_plot(r_quant, howmany = 8, ylim_all = 0.8)
	dev.off()
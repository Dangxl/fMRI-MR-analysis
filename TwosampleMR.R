### R scripts
library(mr.raps)
library(TwoSampleMR)
library(MRPRESSO)
library(RadialMR)
library(ggplot2)

setwd("/rsfMRI-MR")
exposure_data <- read_exposure_data("exposure.txt",
	snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1", 
	other_allele_col="A2", pval_col = "P",  sep="\t")

outcome_data <-  read_outcome_data("outcome.txt", 
	snp_col="SNP", beta_col="BETA", se_col="SE", effect_allele_col="A1", 
	other_allele_col="A2", pval_col = "P", sep="\t")

# Clumping
exp <- clump_data(exposure_data, clump_kb = 1000, clump_r2 = 0.001)

# Harmonization
dat <- harmonise_data(exp, outcome_data)
dat <- dat[which(dat$mr_keep==TRUE),]

# test outliers by using ivw_radial
outlier <- ivw_radial(dat)

# remove outliers
outliers <- outlier$outliers
dat <- dat[which(!(dat$SNP %in% outliers$SNP)),]

## Horizontal_pleiotropy by using MR_PRESSO
dat_pleio <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",
        SdOutcome="se.outcome",SdExposure="se.exposure",data = dat,OUTLIERtest = TRUE)
# remove pleiotropic SNP
dat_pleio$`MR-PRESSO results`$`Outlier Test`$name <- rownames(dat_pleio$`MR-PRESSO results`$`Outlier Test`)
data1 <- subset(dat_pleio$`MR-PRESSO results`$`Outlier Test`, dat_pleio$`MR-PRESSO results`$`Outlier Test`$Pvalue>0.05)
data <- dat[data1$name,]

# F-statistic
# r2 is the proportion of variance in the phenotype explained by the genetic variants
# N is the sample size
data$r2 <- data$BETA*data$BETA/(data$BETA*data$BETA+data$SE*data$SE*data$N)
F_data <- merge(dat, data, by = "SNP")
k <- nrow(F_data)		# the number of instruments

# Calculating the F-statistic
F_data$f <- sum(F_data$r2)*(F_data$N-k-1)/k*(1-sum(F_data$r2))

# Performing Mendelian randomization analyses
res <- mr(dat, method_list = c("mr_ivw_mre", "mr_raps", 
	"mr_weighted_median","mr_weighted_mode", "mr_egger_regression"))

# Transform beta as odd ratios
res <- generate_odds_ratios(res)

het <- mr_heterogeneity(dat)		#heterogeneity test

pleio <- mr_pleiotropy_test(dat)	#pleiotropy test

plot1 <- mr_scatter_plot(res, dat)	#scatter plot
res_single <- mr_singlesnp(dat)		#The effect of single SNP on outcome
plot2 <- mr_forest_plot(res_single)	#Forest plot of single SNP effect size
plot3 <- mr_funnel_plot(res_single)	#Funnel plot
single <- mr_leaveoneout(dat)		#leave-one-out sensitivity analysis
plot4 <- mr_leaveoneout_plot(single)	#Forest plot of leave-one-out sensitivity analysis

# Result graph presentation
data1 <- res
annatation_P  <-  data.frame(x_P = c(rep(4.5,5)))		#The position of P_value 
y <- 1:5
Label_P <- as.character(format(data1$pval))  
a <- paste(data1$or, data1$or_lci95, sep = "[")
b <- paste(a, data1$or_uci95, sep = ",")
c <- paste(b, seq = "]")           
d <- gsub(" ", "", c)              
annatation_OR <- data.frame(x_OR = c(rep(3.2, 5)))		#The position of OR[95% CIs] 
annatation_method <- data.frame(x_method = c(rep(-1, 5)))		#The position of methods
annatation_SNP <- data.frame(x_SNP = c(rep(-2, 5)))	#The position of SNP number
annatation_outcome <- data.frame(x_outcome = c(rep(-3, 5)))		#The position of outcome
annatation_exposure <- data.frame(x_exposure = c(rep(-4, 5)))		#The position of exposure

p<-ggplot(data1, aes(or, method)) +
	geom_point(size=3, shape = 15,colour = "black") +
	geom_errorbarh(aes(xmax =or_uci95, xmin = or_lci95), height = 0.15) +
	scale_x_continuous(limits= c(-4.25, 5.4), breaks= seq(0, 2.5, 0.5)) + 
	geom_vline(aes(xintercept = 1), linetype = "dashed") + 
	xlab('OR') +
	ylab(' ') +
	theme_classic() +
	theme(legend.position = "none") +
	geom_text(data=annatation_P,aes(x=x_P,y=y,label=Label_P)) +
	geom_text(data=annatation_OR,aes(x=x_OR,y=y,label=d)) +
	geom_text(data=annatation_SNP,aes(x=x_SNP,y=y,label=data1$nsnp)) +
	geom_text(data=annatation_method,aes(x=x_method,y=y,label=data1$method)) +
	geom_text(data=annatation_exposure,aes(x=x_exposure,y=y,label=data1$exposure)) +
	geom_text(data=annatation_outcome,aes(x=x_outcome,y=y,label=data1$outcome))+
	theme(axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
		axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +  
	geom_segment(aes(x=-4.25, xend=5, y=-0.6, yend=-0.6), colour="grey") +
	geom_segment(aes(x=-4.25, xend=5, y=5.9, yend=5.9)) + 
	geom_segment(aes(x=0, xend=2.5, y=-0.6, yend=-0.6)) 

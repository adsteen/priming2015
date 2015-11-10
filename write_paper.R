###################
# Script to reproduce data analysis and figures from Steen et al paper, 2015ish, about priming in an estuary
###################

#############################
### Load necessary packages and functions and set options
#############################
require(ggplot2)
library(extrafont)
require(plyr)
require(lubridate)
library(grid)
library(RColorBrewer)
require(reshape2)
require(gridExtra)
require(scales)
require(propagate)

# Load custom functions
source("lm_stats.R")
source("perm_priming.R")
source("get_nls_coefs.R")
source("deterministic_jitter.R")

# Functions to process EEM data. These will go into their own package, someday
source("read_EEM_dat.R")
source("read_Raman.R")
source("f4correct.R")
source("read_UV.R")
source("generate_mask.R")
source("plot_EEM.R")
source("trim_EEM.R")
source("check_EEM.R")
source("calc_indices.R")

# Set graphical theme
theme_set(theme_bw() +
            theme(panel.grid=element_blank(),
                  text=element_text(size=9)))
ps <- 2 # Point size in some plots

# Set plotting parameters
onecol <- 3.34646 #inches; Frontiers Marine Science
twocol <- 7.08661 # inches; Frontiers Marine Science
latexcol <- 5.4 #inches, roughly
global.dpi <- 900 

###################################################
### Process data from phytoplankton day pre-incubation
###################################################

# Read source 14C data
phyto <- read.csv("phytoDecay.csv", stringsAsFactors=FALSE)
phyto$elapsed <- as.numeric(mdy(phyto$date) - min(mdy(phyto$date)))/(60*60*24) # Calculate elapsed time in days (convert from seconds)

# Normalize decay data to initial average value
phyto$C.norm <- phyto$bulkOC / mean(phyto$bulkOC[phyto$elapsed==0], na.rm=TRUE)

# Create a linear least-squares model - this is equation 1 of the ms
phyto_model <- nls(C.norm ~ exp(-1*(log(2)/t12)*elapsed), phyto, start=list(t12=45))

# Try a 2-g model - but discard it, because it predicts a negative R
# phyto_model_2g <- nls(C.norm ~ exp(-1*(log(2)/t12)*elapsed) + R, phyto, start=list(t12=30, R=0.3))

# Define the domain for the model fit
phyto_pred_domain <- 0:max(phyto$elapsed)
phyto_pred <- data.frame(elapsed=phyto_pred_domain, predicted=predict(phyto_model, newdata=data.frame(elapsed=phyto_pred_domain)))

# Determine how much was left at the end
mean_bulk_remaining <- mean(phyto$C.norm[phyto$elapsed==max(phyto$elapsed)])
sd_bulk_remaining <- sd(phyto$C.norm[phyto$elapsed==max(phyto$elapsed)])

# Create Figure 1
p_phyto <- ggplot(phyto, aes(x=elapsed, y=C.norm)) + 
  geom_point() +
  #geom_line(data=phyto_pred, aes(x=elapsed, y=predicted)) +
  xlab("incubation time, days") +
  ylab(expression(paste(phantom(0)^{14}, "OC relative to initial"))) +
  scale_y_continuous(labels=percent, breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  expand_limits(y=0) 
print(p_phyto)
# Note: arguments to ggsave can be platform specific
# ggsave("Fig1.pdf", p_phyto, height=2, width=onecol, units="in", dpi=global.dpi)


###################################################
### Load and pre-process data from main experiment
###################################################

# Load the data
cc <- c("character", "numeric", "character", "factor", "factor", "factor", "factor", "numeric", "numeric", "numeric")
all_d <- read.csv("master_priming_dataset.csv", colClasses=cc)
all_d$Treatment <- factor(all_d$Treatment, labels=c("+N+P", "+acetate", "+BSA+P", "+N", "control"))
all_d$Treatment <- factor(all_d$Treatment, levels=c("+acetate", "+BSA+P", "+N", "+N+P", "control"), ordered=TRUE) # Change order

# Calculate "elapsed days"
all_d$Date <- mdy(all_d$Date)
all_d$Day <- as.numeric(all_d$Date - min(all_d$Date))/(3600*24)

# Calibrate DPM to C units (uM C)
all_d$uM.C <- all_d$DPM1 * 2000 * 1000 / 9.272e6
attr(all_d$uM.C, "units") <- "umol per L"

#######################
# Process 14C data
#######################

# Create rad-only data frame
rad <- all_d[!is.na(all_d$uM.C), 
             c("Date", "Day", "Time", "Treatment", "Rep", "Test", "uM.C")]

### Eliminate bad datapoints
# Calculate means of reps
rad <- ddply(rad, c("Day", "Treatment", "Test"), transform,
                mean.uM.C = mean(uM.C, na.rm=TRUE))

# Detrend
rad$detrended <- rad$uM.C - rad$mean.uM.C

# Calculate groupwise standard deviations
rad <- ddply(rad, c("Test"), transform, groupwise.sd=sd(detrended, na.rm=TRUE))

# Determine which points >2.5 standard deviations from the detrended mean
rad$eliminate <- abs(rad$detrended) >= 2.5*rad$groupwise.sd

# Calculate how many points were dropped
rad_to_drop <- ddply(rad, .(Test), summarise, to.drop=sum(eliminate, na.rm=TRUE))

# Actually drop the bad points, using needlessly opaque syntax
rad <- rad[!rad$eliminate | is.na(rad$eliminate), ]
rad <- rad[-which(rad$uM.C > 1200 & (rad$Test=="bulk" | rad$Test=="poc")), ]

# Create CO2-only data frame (with filtered points)
rad_CO2 <- rad[rad$Test=="co2", ]

### Fit the 14C data to theoretical equations

# Identify emperical initial POC and bulk OC concentrations (in paper, I call bulk "total")
initial_POC <- mean(subset(rad, Test=="poc" & Day==0)$uM.C)
initial_bulk <- mean(subset(rad, Test=="bulk" & Day==0)$uM.C)

# Calculate the nls fits
bulk_fits <- dlply(subset(rad, Test=="bulk"), 
                   c("Test", "Treatment"),
                  function(x) nls(formula=formula(I(uM.C ~ (initial_bulk-recalcitrant)*exp(-1*k*Day) + recalcitrant)), 
                                  data=x,start=list(k=0.25, recalcitrant=500)))

poc_fits <- dlply(subset(rad, Test=="poc"), 
                  c("Test", "Treatment"),
                  function(x) nls(formula=formula(I(uM.C ~ (initial_POC-recalcitrant)*exp(-1*k*Day) + recalcitrant)), data=x,
                                  start=list(k=0.25, recalcitrant=500)))
# Concatenate the bulk fits with the poc fits into one data frame
all_fits <- c(bulk_fits, poc_fits)

# Get the coefficients out of the nls fits object
fit_coefs <- ldply(all_fits, get_nls_coefs)
ids <- ldply(strsplit(fit_coefs$.id, ".", fixed=TRUE), identity)
names(ids) <- c("Test", "Treatment")
fit_coefs <- cbind(fit_coefs, ids)

# Make a data-frame format grid for creating model predictions
pred_x <- data.frame(Day=0:max(rad$Day))

# Function to make the predictions
make_nls_prediction <- function(x) {
  pred_vec <- predict(x, newdata=pred_x)
  data.frame(Day=pred_x, pred=pred_vec)
}

# Put all the fits into one data frame, and parse the .id column into treatment and test
pred_data <- ldply(all_fits, make_nls_prediction)
# Rename the .id columsn
ids2 <- ldply(strsplit(pred_data$.id, ".", fixed=TRUE), identity)
names(ids2) <- c("Test", "Treatment")
pred_data <- cbind(pred_data, ids2)

# Make separate data frames for:
#  -actual data, with controls stripped out
#  -actual data, controls only ($Treatment column stripped out)
#  -prediction lines, with controls stripped out
#  -prediction lines, controls only ($Treatment column stripped out)
OC_control <- subset(rad, (Test=="bulk" | Test=="poc") & Treatment=="control")
OC_control <- OC_control[, -which(names(OC_control) %in% "Treatment")]
pred_control <- subset(pred_data, Treatment=="control")
pred_control <- pred_control[ , -which(names(pred_control) %in% "Treatment")]
OC_treat <- subset(rad, (Test=="bulk" | Test=="poc") & Treatment!="control")
pred_treat <- subset(pred_data, Treatment != "control")

# Rename "bulk" to "total" in each of the component data frames of Fig 2
OC_control_to_print <- OC_control
OC_control_to_print$Test <- revalue(OC_control_to_print$Test, c("bulk"="total", "poc" = "POC"))
OC_treat_to_print <- OC_treat
OC_treat_to_print$Test <- revalue(OC_treat_to_print$Test, c("bulk" = "total", "poc" = "POC"))
pred_treat_to_print <- pred_treat
pred_treat_to_print$Test <- revalue(pred_treat_to_print$Test, c("bulk" = "total", "poc" = "POC"))
pred_control_to_print <- pred_control
pred_control_to_print$Test <- revalue(pred_control_to_print$Test, c("bulk" = "total", "poc" = "POC"))

### Make Figure 2
p_OC <- ggplot() + 
  geom_point(data=OC_control_to_print, aes(x=Day, y=uM.C), shape=1) +
  geom_point(data=OC_treat_to_print, aes(x=Day, y=uM.C), shape=16) +
  geom_line(data=pred_treat_to_print, aes(x=Day, y=pred), linetype=1) + 
  geom_line(data=pred_control_to_print, aes(x=Day, y=pred), linetype=2) +
  xlab("elapsed time, days") +
  ylab(expression(paste(mu, "M C"))) +
  expand_limits(y=0) + 
  facet_grid(Test~Treatment) 
print(p_OC)
# ggsave("Fig2.pdf", p_OC, height=3, width=latexcol, units="in", dpi=global.dpi)

# Create Table 1
fit_coefs_tab <- rename(fit_coefs, c(".id"="treatment"))[ , c("Treatment", "k", "k.se", "Test", "recalcitrant", "recalcitrant.se")]
fit_coefs_tab$k <- paste(signif(fit_coefs_tab$k, 2), "+/-", signif(fit_coefs_tab$k.se, 2))
fit_coefs_tab$R <- paste(signif(fit_coefs_tab$recalcitrant, 2), "+/-", signif(fit_coefs_tab$recalcitrant.se, 2))
ftcrm <- melt(fit_coefs_tab[ , c("Treatment", "k", "Test", "R")], id.vars=c("Treatment", "Test"))
ftcrc <- dcast(ftcrm, formula=Test+Treatment ~ variable)
print("Table 1")
print(ftcrc)

###################################################
### Calculate the extent of priming
###################################################
# Correct offset in CO2 data, evident from the start
offset_guesses <- ddply(subset(rad_CO2, Day<=5), c("Treatment"), summarise, intercept=coefficients(lm(uM.C ~ Day))[1])
offset_guess <- offset_guesses$intercept[offset_guesses$Treatment=="control"] #About 14.4
# Subtract offset value
CO2_offset <- rad_CO2
CO2_offset$uM.C <- CO2_offset$uM.C - offset_guess

# Model CO2 prodution in each treatment
exp_form <- formula(I(
  uM.C ~ A * (1 - exp(-1 * lambda * Day))
  ))

# Create nls estimates with a for loop. Trying it with dlply throws an error I don't understand
nls_estimates <- data.frame(Treatment=NULL, A=NULL, lambda=NULL, A.sd=NULL, lambda.sd=NULL)
nls_objects <- list()
CO2_pred <- list()

# Define domain for predictions
pred_domain <- data.frame(Day=0:36)

###
### THIS BLOCK TAKES 2-ISH MINUTES TO EXECUTE ON MY LAPTOP
###
for(i in levels(CO2_offset$Treatment)) {# The issue (I think) is whether I use package predictNLS or my modified verison
  # Create NLS model
  m <- nls(uM.C ~ A * (1 - exp(-1 * lambda * Day)), 
           subset(CO2_offset, Treatment==i), start=list(A=300, lambda=1/20))
  
  # Pull out summary statistics (including standard errors) in a better format than coef() gives
  sm <- summary(m)
  nls_estimates <- rbind(nls_estimates, data.frame(Treatment=i, 
                             A=summary(m)$coefficients[1, 1],
                             lambda=summary(m)$coefficients[2, 1],
                             A.sd=summary(m)$coefficients[1, 2],
                             lambda.sd=summary(m)$coefficients[2, 2]))
  
  # Place the nls objects in a list and name them appropriately
  nls_objects[[i]] <- m
  names(nls_objects[i]) <- i
  
  # use predictNLS to get prediction intervals
  CO2_predict_df <- predictNLS(m, newdata=pred_domain)$summary
  #CO2_predict_df <- predictNLS(m, )$summary
  CO2_pred[[i]] <- cbind(pred_domain, CO2_predict_df)
}

# Make plots of the nls parameter estimates - in the manuscript I put these in a table
p_lambda <- ggplot(nls_estimates, aes(x=Treatment, y=lambda)) + 
  geom_pointrange(aes(ymin=lambda-lambda.sd, ymax=lambda+lambda.sd))
p_A <- ggplot(nls_estimates, aes(x=Treatment, y=A)) + 
  geom_pointrange(aes(ymin=A-A.sd, ymax=A+A.sd))

### Table 2: Make printable table of nls estimates
nls_est_toprint <- nls_estimates
nls_est_toprint$A.toprint <- paste(signif(nls_est_toprint$A, 3), "+/-", signif(nls_est_toprint$A.sd, 3))
nls_est_toprint$lambda.toprint <- paste(signif(nls_est_toprint$lambda, 3), "+/-", signif(nls_est_toprint$lambda.sd, 3))
nls_est_toprint <- nls_est_toprint[ , c(1, 6, 7)]
nls_est_toprint <- rename(nls_est_toprint, c("A.toprint" = "A", "lambda.toprint"="lambda"))
print("Table 2:")
print(nls_est_toprint)

# Collapse the predictions data frame into a single data frame
CO2_pred_df <- ldply(CO2_pred, identity)
CO2_pred_df <- rename(CO2_pred_df, c(".id"="Treatment", "Sim.Mean"="uM.C.pred", "Sim.sd"="uM.C.pred.sd"))

# Create data frames with just control data and no "Treatment" column
CO2_control <- subset(CO2_offset, Treatment=="control")[ , !(names(CO2_offset) %in% "Treatment")]
control_fit <- subset(CO2_pred_df, Treatment=="control")[ , !(names(CO2_pred_df) %in% "Treatment")]

# Make 'raw' plot of CO2 in treatment vs control
# ps <- 2
# p_CO2_production <- ggplot() +
#   geom_ribbon(data=control_fit, aes(x=Day, ymin=uM.C.pred-uM.C.pred.sd, ymax=uM.C.pred+uM.C.pred.sd), alpha=0.2) + #predicted error ribbon for control
#   geom_ribbon(data=subset(CO2_pred_df, Treatment!="control"), aes(x=Day, ymin=uM.C.pred-uM.C.pred.sd, ymax=uM.C.pred+uM.C.pred.sd), alpha=0.2) +
#   geom_line(data=control_fit, aes(x=Day, y=uM.C.pred), linetype=2) + #control modeled line
#   geom_point(data=CO2_control, aes(x=Day, y=uM.C), shape=1, size=ps) + # control points
#   geom_point(data=subset(CO2_offset, Treatment!="control"), aes(x=Day, y=uM.C, shape=Treatment, colour=Treatment), size=ps) + #data points
#   geom_line(data=subset(CO2_pred_df, Treatment!="control"), aes(x=Day, y=uM.C.pred)) + #data modeled line
#   scale_colour_brewer(palette="Dark2") +
#   xlab("time elapsed, days") +
#   ylab(expression(paste(CO[2], " evolved, ", mu, "M C"))) +
#   facet_wrap(~Treatment, nrow=1) +
#   theme(legend.position="top")
# print(p_CO2_production)
# Calculate priming
# Summary of offset-corrected CO2 data
CO2_summ <- ddply(CO2_offset, c("Day", "Treatment"), summarise,
                  mean.uM.C = mean(uM.C), 
                  sd.uM.C = sd(uM.C))
# Drop days with no control values
CO2_summ <- ddply(CO2_summ, c("Day"), transform, 
                  n.control = length(mean.uM.C[Treatment=="control"]))
CO2_summ <- subset(CO2_summ, n.control > 0)
CO2_summ <- ddply(CO2_summ, c("Day"), transform,
                  control.CO2 = mean.uM.C[Treatment=="control"],
                  control.CO2.sd = sd.uM.C[Treatment=="control"])

# Actually calculate priming
CO2_summ$priming <- CO2_summ$mean.uM.C / CO2_summ$control.CO2 - 1

# Calculate priming based on smoothed data
CO2_pred_slim <- CO2_pred_df[ , c("Treatment", "Day", "uM.C.pred", "uM.C.pred.sd")]
CO2_pred_slim <- ddply(CO2_pred_slim, c("Day"), transform, modeled.priming=(uM.C.pred/uM.C.pred[Treatment=="control"])-1)
CO2_pred_slim <- subset(CO2_pred_slim, Treatment!="control")

#modeled_priming <- merge(CO2_pred_df, control_fit, by="Day") # merge summary data with control fitted priming data
#modeled_priming$mod.prim <- modeled_priming$uM.C.x / modeled_priming$uM.C.y - 1
#modeled_priming <- subset(modeled_priming, Treatment!="control")

# # # Make plot of priming
# tiff("plots/front_mar_sci_figs/priming_plot_for_print.tiff", height=4, width=twocol, units="in", res=300, compression="lzw", type="cairo")
# print(p_CO2_production, vp=vp_prim1)
# print(p_priming, vp=vp_prim2)
# dev.off()


### Simulate power of priming calculation: how much priming is statistically significant?
### (Uses permutation approach; takes about a minute on my machine)
# number of bootstrap reps
n <- 1000 #You'd get a decent estimate with a few hundred

# initialize a list
sim_fits <- list()
set.seed(212)
# Populate it with data frames, each with a different set of simulations
system.time({
  for (i in 1:n) {
    #sim_fits[[i]] <- calc_mean_diffs(CO2_offset) 
    sim_fits[[i]] <- perm_priming(CO2_offset)
  }
})

# Put the simulations back together into a data frame
all_sim_fits <- ldply(sim_fits, identity) # There's a lone NA in there somewhere

# Determine 95% confidence intervals for variation between spline of treatment and control
sim_probs <- ddply(all_sim_fits, c("Day", "Treatment"), function(x) quantile(x$sim.priming, c(0.025, 0.5, 0.975), na.rm=TRUE))
sim_probs <- rename(sim_probs, c("2.5%"="low", "50%"="median", "97.5%"="high"))
sim_probs_m <- melt(sim_probs, id.vars=c("Day", "Treatment"), variable.name="quantile", value.name="cutoff")

# Where do the actual spline values sit with respect to the "confidence intervals"
all_the_priming <- merge(CO2_pred_slim, sim_probs, by=c("Day", "Treatment"))

# Make a column showing whether the modeled priming is significantly different from zero
all_the_priming$is.signif <- (all_the_priming$modeled.priming < all_the_priming$low) |
  (all_the_priming$modeled.priming > all_the_priming$high)

# Manually look at when priming was distinguishable from zero
subset(all_the_priming, Treatment=="+acetate")[order(subset(all_the_priming, Treatment=="+acetate")$Day, decreasing=FALSE), ]
subset(all_the_priming, Treatment=="+BSA+P")[order(subset(all_the_priming, Treatment=="+BSA+P")$Day, decreasing=FALSE), ]
subset(all_the_priming, Treatment=="+N+P")[order(subset(all_the_priming, Treatment=="+N+P")$Day, decreasing=FALSE), ]
subset(all_the_priming, Treatment=="+N")[order(subset(all_the_priming, Treatment=="+N")$Day, decreasing=FALSE), ]

##### Print Fig 3
p_CO2_bw <- ggplot() +
  geom_ribbon(data=control_fit, aes(x=Day, ymin=uM.C.pred-uM.C.pred.sd, ymax=uM.C.pred+uM.C.pred.sd), alpha=0.4) + #predicted error ribbon for control
  geom_ribbon(data=subset(CO2_pred_df, Treatment!="control"), aes(x=Day, ymin=uM.C.pred-uM.C.pred.sd, ymax=uM.C.pred+uM.C.pred.sd), alpha=0.4) +
  geom_line(data=control_fit, aes(x=Day, y=uM.C.pred), linetype=2) + #control modeled line
  geom_point(data=CO2_control, aes(x=Day, y=uM.C), shape=1, size=ps) + # control points
  #geom_point(data=subset(CO2_offset, Treatment!="control"), aes(x=Day, y=uM.C, shape=Treatment), size=ps) + #data points
  geom_point(data=subset(CO2_offset, Treatment!="control"), aes(x=Day, y=uM.C), size=ps) + #data points
  geom_line(data=subset(CO2_pred_df, Treatment!="control"), aes(x=Day, y=uM.C.pred)) + #data modeled line
  #scale_colour_brewer(palette="Dark2") +
  xlab("incubation time, days") +
  ylab(expression(paste(CO[2], " evolved, ", mu, "M C"))) +
  facet_wrap(~Treatment, nrow=1) +
  theme(legend.position="top")

p_priming_bw <- ggplot() + 
  #geom_point(aes(shape=Treatment), size=ps) + 
  geom_point(data=subset(CO2_summ, Treatment != "control"), aes(x=Day, y=priming), size=ps) + 
  geom_line(data=CO2_pred_slim, aes(x=Day, y=modeled.priming), linetype=1) +
  geom_hline(yintercept=0, linetype=2) +
  geom_ribbon(data=subset(sim_probs, Treatment != "control"), aes(x=Day, ymin=low, ymax=high), fill="black", alpha=0.4) +
  scale_y_continuous(labels=percent) + 
  xlab("incubation time, days") +
  facet_wrap(~Treatment, nrow=1) +
  theme(legend.position="none")

# Define subplot window size
lower_height=0.5
vp_prim1 <- viewport(x=0.5, y=lower_height+(1-lower_height)/2, width=1, height=(1-lower_height))
vp_prim2 <- viewport(x=0.5, y=lower_height/2, width=1, height=lower_height)

# Print the plots
# tiff("Fig3.tiff", height=4, width=twocol, units="in", res=global.dpi, compression="lzw", type="cairo")
print(p_CO2_bw, vp=vp_prim1)
print(p_priming_bw, vp=vp_prim2)
# dev.off()

# ##### For talks
# p_CO2_talk <- p_CO2_production +
#   theme(text=element_text(size=18)) + 
#   scale_colour_brewer(palette="Dark2")
# p_priming_talk <- p_priming +
#   theme(text=element_text(size=18)) +
#   scale_colour_brewer(palette="Dark2")
# png("../plots/2014_05_19_priming_for_talk.png", height=5, width=8, units="in", res=300)
# print(p_CO2_talk, vp=vp_prim1)
# print(p_priming_talk, vp=vp_prim2)
# dev.off()


###################################################
### Cell counts
###################################################

# Load and analyze cell data (not a lot of analysis going on here)
cells <- read.csv("priming_ms_cell_cts_w_error.csv")
cells$Date <- mdy(cells$Date)
cells$Day <- as.numeric(cells$Date-min(cells$Date))/(3600*24) + 1 
cells$error <- cells$Cell.Numbers*(cells$Percent.Error/100)

# This makes more sense in color
cols <- c(brewer.pal(4, "Dark2"), "#000000")
p_cells_color <- ggplot(cells, aes(x=Day, y=Cell.Numbers, shape=Treatment, colour=Treatment)) + 
  #geom_vline(xintercept=enz_days, colour="gray50") +
  geom_pointrange(size=ps, aes(ymin=Cell.Numbers-error, ymax=Cell.Numbers+error), linetype=1)  +
  geom_line() +
  xlab("incubation time, days") +
  ylab(expression(paste("cells ", ml^{-1}))) +
  scale_x_continuous(breaks=seq(from=0, to=50, by=10)) +
  #scale_y_log10(limits=c(5e5, 5e7), breaks=c(5e5, 1e6, 1e7, 5e7)) +
  scale_y_continuous(breaks=c(0, 2.5e6, 5e6, 7.5e6, 1e7)) +
  scale_colour_manual(values=cols) +
  scale_shape_manual(values=c(16, 17, 15, 3, 1)) +
  #annotation_logticks(sides="l") +
  theme(legend.position="bottom") +
  guides(shape=guide_legend(nrow=2),
         color=guide_legend(nrow=2))
print(p_cells_color)
# ggsave("Fig4.tiff", p_cells_color, height=2.5, width=onecol, units="in", dpi=global.dpi, compression="lzw", type="cairo")

###################################################
### Enzymes
###################################################
# Define fl-only data frame
fl <- all_d[all_d$Test=="enzyme" & !is.na(all_d$Test), ]

# Add a column for true time
fl$Rtime <- fl$Date + hm(fl$Time)

# Determine elapsed time for each incubation
fl <- ddply(fl, c("Date", "Treatment", "Subs"), transform, elapsed=as.numeric(Rtime-min(Rtime))/3600)

# Make plots of raw fluorescence data - this is a QC check for me
p_enz_raw <- ggplot(fl, aes(x=elapsed, y=FL, colour=Subs, shape=Rep)) +
  geom_point() + 
  geom_smooth(method="lm", se=FALSE) + 
  facet_grid(Day~Treatment)
#print(p_enz_raw)

# Calculate uncalibrated activities
slopes <- ddply(fl, c("Date", "Treatment","Subs", "Rep"), function(x) lm_stats(x, xvar="elapsed", yvar="FL"))
slopes$incubation.time <- as.numeric(slopes$Date-min(all_d$Date))/(3600*24) 
attr(slopes$incubation.time, "units") <- "days"

# Drop outlier in BSA data
slopes <- subset(slopes, !((incubation.time==7 | incubation.time==16) & Treatment=="+BSA+P" & Subs=="leu-AP" & slope < 0))

# Calibrate activities
AMC_calib <- read.csv("fl_calibration_data_AMC.csv")
AMC_calib$std <- "AMC"
MUB_calib <- read.csv("fl_calibration_data_MUB.csv")
MUB_calib$std <- "MUB"
calib_raw <- rbind(AMC_calib, MUB_calib)

# Plot calibration curves, for reference
p_calib <- ggplot(calib_raw, aes(x=conc.uM, y=fl)) + 
  geom_point() + 
  geom_smooth(data=subset(calib_raw, conc.uM<=10), aes(x=conc.uM, y=fl), method="lm") + #That's quite good.
  facet_wrap(~std)
# print(p_calib)

calib_coefs <- ddply(subset(calib_raw, conc.uM <= 10), c("std"), summarise, slope=coefficients(lm(fl~conc.uM))[2])
slopes$calib.slope <- NA
slopes$calib.slope[slopes$Subs=="b-glu" | slopes$Subs=="PO4"] <- slopes$slope[slopes$Subs=="b-glu" | slopes$Subs=="PO4"]/calib_coefs$slope[2] * 1000 #MUB
slopes$calib.slope[slopes$Subs=="leu-AP"] <- slopes$slope[slopes$Subs=="leu-AP"]/calib_coefs$slope[1] * 1000 #AMC
attr(slopes$calib.slope, "units") <- "nm per hr"
#

# Calculate summarises of calibration slopes
sum_slopes <- ddply(slopes, c("incubation.time", "Treatment", "Subs"), summarise, 
                    mean.calib.slope=mean(calib.slope, na.rm=TRUE),
                    sd.calib.slope=sd(calib.slope, na.rm=TRUE))

# Slightly jitter points
sum_slopes2 <- sum_slopes
wd=0.15

# The "deterministic jitter" here keeps the error bars off from on top of each other
sum_slopes2$incubation.time <- sum_slopes2$incubation.time - wd*(as.numeric(sum_slopes2$Treatment)-(max(as.numeric(sum_slopes2$Treatment))+1) / 2)
sum_slopes2$incubation.time <- sum_slopes$incubation.time - deterministic_jitter(sum_slopes$Treatment, width=0.15)

# Make plot of calibrated enzyme activities
p_enz <- ggplot(sum_slopes2, aes(x=incubation.time, y=mean.calib.slope, colour=Treatment, shape=Treatment)) + 
  geom_point() + 
  geom_line() +
  #geom_point(data=slopes, aes(x=incubation.time, y=slope, colour=Treatment)) +
  geom_errorbar(aes(ymin=mean.calib.slope-sd.calib.slope, ymax=mean.calib.slope+sd.calib.slope), width=0) +
  #geom_errorbar(aes(ymin=slope-slope.se, ymax=slope+slope.se)) + 
  scale_colour_manual(values=cols) +
  xlab("incubation time, days") +
  ylab(expression(paste(V[max], ", nM ", hr^{-1}))) +
  facet_wrap(~Subs, nrow=1, scales="free")
print(p_enz)
# ggsave("Fig5.tiff", height=2., width=twocol, units="in", dpi=global.dpi, compression="lzw", type="cairo")

##############
# Fluorescence data processing (EEMs)
##############

# Need slightly different graphical theme for EEMs data
theme_set(theme_bw() + 
            theme(panel.grid=element_blank()))

###### NEED TO FILTER EEMS SO THAT WAVELENGTHS > 280 (or so; check)


base_path <- "EEMs/"

# Create a generic sample with which to generate a mask
#samp <- read_EEM_dat("L1AfinalD50 (01)_Graph_S1_R1.dat") # THis yields a data frame correct row names

# Trim the sample due to absorbance of methacrylate below a set limit
#low.lim <- 300 # This seems conservative, 295 or 290 would also probably be OK
#samp_trimmed <- trim_EEM(samp, low.ex=low.lim, low.em=low.lim)

# Generate mask
#mask <- generate_mask(samp)


###################################################
### code chunk number 3: read_Raman
###################################################
# read Raman file
Raman <- read_Raman(paste(base_path, "FLRaman/DfltEm (01)_Graph.dat", sep=""))

### Read UV files
# list of UV filenames
UV_list <- list(init3A=read_UV(paste(base_path, "UVFiles/UV1.csv", sep="")),
                init3A_d50=read_UV(paste(base_path, "UVFiles/UV2.csv", sep="")), #Actual Title is 3Ainitial2d, which I THINK means 2x diluted
                init5A_d50=read_UV(paste(base_path, "UVFiles/UV3.csv", sep="")),
                final1A_d50=read_UV(paste(base_path, "UVFiles/UV4.csv", sep="")),
                final2A_d50=read_UV(paste(base_path, "UVFiles/UV5.csv", sep="")),
                final3A_d50=read_UV(paste(base_path, "UVFiles/UV6.csv", sep="")),
                final4A_d50=read_UV(paste(base_path, "UVFiles/UV7.csv", sep="")),
                final5A_d50=read_UV(paste(base_path, "UVFiles/UV8.csv", sep="")),
                Lblank=read_UV(paste(base_path, "UVFiles/UV9.csv", sep="")),
                final1B_d50=read_UV(paste(base_path, "UVFiles/UV10.csv", sep="")),
                final2B_d50=read_UV(paste(base_path, "UVFiles/UV11.csv", sep="")),
                final3B_d50=read_UV(paste(base_path, "UVFiles/UV12.csv", sep="")),
                final4B_d50=read_UV(paste(base_path, "UVFiles/UV13.csv", sep="")),
                final5B_d50=read_UV(paste(base_path, "UVFiles/UV14.csv", sep=""))
)

# Read the UV files
UV_df <- ldply(UV_list, identity)

# Plot the UV_files
p_UV <- ggplot(UV_df, aes(x=wavelength, y=abs)) +
  geom_line() +
  geom_vline(xintercept=295) +
  ylim(c(0, 0.01)) +
  facet_wrap(~.id)
print(p_UV)
# 5BFinal must be bad. (Big hump above 500 nm)
# Fisher says these are transparent above 285 nm. I'll remove everything below 290 to be safe.

### Read all the sample (fluorescence) files
# Trim files to below 290 nm
low.lim <- 290 #300
trim_UV <- function(UV, low_limit=low.lim) {
  UV <- UV[UV$wavelength >= low_limit, ]
}

# Actually trim the UV files
UV_list <- llply(UV_list, trim_UV)

# Read the blank for each sample
# Remove final5B_d50
raw_path <- paste(base_path, "RawEEMs/", sep="")
blank_list <- list(init5A = read_EEM_dat("B5Ainitial (01)_Graph_S1_R1.dat", path=raw_path),
                   final1A = read_EEM_dat("B1Afinal (01)_Graph_S1_R1.dat", path=raw_path),
                   final2A = read_EEM_dat("B2Afinal (01)_Graph_S1_R1.dat", path=raw_path),
                   final2B = read_EEM_dat("B2Bfinal (01)_Graph_S1_R1.dat", path=raw_path),
                   final3A = read_EEM_dat("B3Afinal (01)_Graph_S1_R1.dat", path=raw_path),
                   final3B = read_EEM_dat("B3Bfinal (01)_Graph_S1_R1.dat", path=raw_path), 
                   final4A = read_EEM_dat("B4Afinal (01)_Graph_S1_R1.dat", path=raw_path),
                   final4B = read_EEM_dat("B4Bfinal (01)_Graph_S1_R1.dat", path=raw_path),
                   final5A = read_EEM_dat("B5Afinal (01)_Graph_S1_R1.dat", path=raw_path),
                   final5B = read_EEM_dat("B5Bfinal (01)_Graph_S1_R1.dat", path=raw_path)
)

# Trim the blanks to the low limit or greater
trimmed_blanks <- llply(blank_list, trim_EEM, low.ex=low.lim, low.em=low.lim)


# Mask the blanks to see what the blanks look like by themselves
masked_blanks <- llply(trimmed_blanks, function(x) x*generate_mask(trimmed_blanks[[1]], mask.Raman=TRUE))

# Make a data frame of masked blanks
masked_blanks_df <- ldply(masked_blanks, melt, value.name="fl")

# Plot of the blanks
p_all_blanks <- ggplot(masked_blanks_df, aes(x=ex, y=em, fill=fl)) + 
  geom_raster() + 
  stat_contour(aes(z=fl)) + 
  scale_fill_gradientn(colours=rainbow(7)) + 
  xlab(expression(paste(lambda[ex]))) +
  ylab(expression(paste(lambda[em]))) +
  coord_cartesian() +
  facet_wrap(~.id)
print(p_all_blanks)


# Read the actual data
raw_sample_EEMS <- list(init3A = read_EEM_dat("L3Ainitiald50 (01)_Graph_S1_R1.dat", path=raw_path),
                        init5A = read_EEM_dat("L5AinitialD50 (01)_Graph_S1_R1.dat", path=raw_path),
                        final1A = read_EEM_dat("L1AfinalD50 (01)_Graph_S1_R1.dat", path=raw_path),
                        final1B = read_EEM_dat("L1Bfinal (01)_Graph_S1_R1.dat", path=raw_path),
                        final2A = read_EEM_dat("L2Afinal (01)_Graph_S1_R1.dat", path=raw_path),
                        final2B = read_EEM_dat("L2Bfinal (01)_Graph_S1_R1.dat", path=raw_path),
                        final3A = read_EEM_dat("L3Afinal (01)_Graph_S1_R1.dat", path=raw_path),
                        final3B = read_EEM_dat("L3Bfinal (01)_Graph_S1_R1.dat", path=raw_path),
                        final4A = read_EEM_dat("L4Afinal (01)_Graph_S1_R1.dat", path=raw_path),
                        final4B = read_EEM_dat("L4Bfinal (01)_Graph_S1_R1.dat", path=raw_path),
                        final5A = read_EEM_dat("L5Afinal (01)_Graph_S1_R1.dat", path=raw_path),
                        final5B = read_EEM_dat("L5Bfinal (01)_Graph_S1_R1.dat", path=raw_path),
                        process_blank = read_EEM_dat("Lblank (01)_Graph_S1_R1.dat", path=raw_path))
trimmed_raw_sample_EEMS <- llply(raw_sample_EEMS, trim_EEM, low.ex=low.lim, low.em=low.lim)
raw_sample_EEMS_mask <- llply(trimmed_raw_sample_EEMS, function(x) x * generate_mask(trimmed_raw_sample_EEMS[[1]], mask.Raman=TRUE))
samp_raw_melt <- ldply(raw_sample_EEMS_mask, melt, value.name="fl")

p_raw_samples <- ggplot(samp_raw_melt, aes(x=ex, y=em, fill=fl)) +
  geom_raster() + 
  scale_fill_gradientn(colours=rainbow(7)) +
  stat_contour(aes(z=fl)) +
  facet_wrap(~.id)
# print(p_raw_samples)

### Process the EEMs
ma <- generate_mask(trimmed_raw_sample_EEMS$init3A, mask.Raman=TRUE)
proc_samp <- list(
  init3A = f4correct(samp=trimmed_raw_sample_EEMS$init3A, 
                     blank=trimmed_blanks$init5A, Raman=Raman, # NOTE I"M USING THE WRONG BLANK HERE
                     UV=UV_list$init3A_d50, mask=ma),
  init5A = f4correct(samp=trimmed_raw_sample_EEMS$init5A, 
                     blank=trimmed_blanks$init5A, Raman=Raman, 
                     UV=UV_list$init5A_d50, mask=ma),
  final1A = f4correct(samp=trimmed_raw_sample_EEMS$final1A, 
                      blank=trimmed_blanks$final1A, Raman=Raman, 
                      UV=UV_list$final1A_d50, mask=ma),
  final2A = f4correct(samp=trimmed_raw_sample_EEMS$final2A, 
                      blank=trimmed_blanks$final2A, Raman=Raman, 
                      UV=UV_list$final2A_d50, mask=ma),
  final2B = f4correct(samp=trimmed_raw_sample_EEMS$final2B, 
                      blank=trimmed_blanks$final2B, Raman=Raman, 
                      UV=UV_list$final2B_d50, mask=ma),
  final3A = f4correct(samp=trimmed_raw_sample_EEMS$final3A, 
                      blank=trimmed_blanks$final3A, Raman=Raman, 
                      UV=UV_list$final3A_d50, mask=ma),
  final3B = f4correct(samp=trimmed_raw_sample_EEMS$final3B, 
                      blank=trimmed_blanks$final3B, Raman=Raman, 
                      UV=UV_list$final3B_d50, mask=ma),
  final4A = f4correct(samp=trimmed_raw_sample_EEMS$final4A, 
                      blank=trimmed_blanks$final4A, Raman=Raman, 
                      UV=UV_list$final4A_d50, mask=ma),
  final4B = f4correct(samp=trimmed_raw_sample_EEMS$final4B, 
                      blank=trimmed_blanks$final4B, Raman=Raman, 
                      UV=UV_list$final4B_d50, mask=ma),
  final5A = f4correct(samp=trimmed_raw_sample_EEMS$final5A, 
                      blank=trimmed_blanks$final5A, Raman=Raman, 
                      UV=UV_list$final5A_d50, mask=ma),
  final5B = f4correct(samp=trimmed_raw_sample_EEMS$final5B, 
                      blank=trimmed_blanks$final5B, Raman=Raman, 
                      UV=UV_list$final5B_d50, mask=ma)
)
proc_samp_m <- ldply(proc_samp, melt, value.name="fl") # L1 is due to the fact that f4correct returns a list (at the moment)

###### Define sample treatment, timepoint and replicate
# Note: I should do this with regular expressions, but http://xkcd.com/1171/
# Parse timepoint
proc_samp_m$timepoint <- NA
proc_samp_m$timepoint[grep("final", proc_samp_m$.id)] <- "final"
proc_samp_m$timepoint[grep("init", proc_samp_m$.id)] <- "init"
proc_samp_m$timepoint <- factor(proc_samp_m$timepoint, levels=c("init", "final"), labels=c("initial", "final"), ordered=TRUE)

# Parse treatment
proc_samp_m$treat.num <- NA
proc_samp_m$treat.num[proc_samp_m$timepoint=="initial"] <- substr(as.character(proc_samp_m$.id[proc_samp_m$timepoint=="initial"]), start=5, stop=5)
proc_samp_m$treat.num[proc_samp_m$timepoint=="final"] <- substr(proc_samp_m$.id[proc_samp_m$timepoint=="final"], start=6, stop=6)
proc_samp_m$treat.num <- as.numeric(proc_samp_m$treat.num)
proc_samp_m$treatment <- factor(proc_samp_m$treat.num, levels=1:5, labels=c("+N+P", "+acetate", "+BSA+P", "+N", "control"))
proc_samp_m$treatment <- factor(proc_samp_m$treatment, levels=c("+acetate", "+BSA+P", "+N", "+N+P", "control", ordered=TRUE))

# Parse replicate
proc_samp_m$replicate <- NA
proc_samp_m$replicate[grep("A", proc_samp_m$.id)] <- "A"
proc_samp_m$replicate[grep("B", proc_samp_m$.id)] <- "B"

# Omit bad spectrum final5B
proc_samp_m <- subset(proc_samp_m, !(timepoint=="final" & treat.num==5 & replicate=="B"))

# Fig 6: Processed EEMs for paper
p_proc <- ggplot(proc_samp_m, aes(x=em, y=ex, fill=fl)) +
  geom_raster() + 
  scale_fill_gradientn(colours=rainbow(7),
                       limits=c(0, 15000)) +
  stat_contour(aes(z=fl)) +
  coord_equal() +
  xlab(expression(paste(lambda[em]))) + 
  ylab(expression(paste(lambda[ex]))) +
  facet_grid(replicate+timepoint~treatment) +
  theme(axis.text.x=element_text(angle=-45, hjust=0),
        text=element_text(size=10))
print(p_proc)
#ggsave("Fig6.pdf", p_proc, height=4, width=latexcol, units="in", dpi=global.dpi)

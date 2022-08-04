# ---------------------------------------------------------------------------- #
#                                                                              #
# R code for the article: "A note on synthetic data for replication purposes   #
# in agricultural economics"                                                   #
#                                                                              #
# Authors: Stefan Wimmer and Robert Finger                                     #
#                                                                              #
# Corresponding author: Stefan Wimmer (swimmer@ethz.ch)                        #
#                                                                              #
# Citation: Wimmer, S., Finger, R. (2023). "A note on synthetic data for       #
# replication purposes in agricultural economics." Journal of Agricultural     #
# Economics.                                                                   #
#                                                                              #
# The program stores all tables and figures based on the FADN data for German  #
# crop farms presented in the appendix in the subfolders 'Tables/FADN'         #
# and 'Figures/FADN'. The data are not publicly available but can be           #
# requested from DG-Agri. The script was successfully run with R version 4.2.0.#
#                                                                              #
# ---------------------------------------------------------------------------- #

# Content:
# 0) Preparation (working directory, packages, data, colours)
# 1) Generation of synthetic data
# 2) Descriptive statistics
# 3) Estimation of production function
# 4) Estimation of production frontier
# 5) Data envelopment analysis

# -------------------- #
#### 0. Preparation ####
# -------------------- #

# Set working directory
my.wd <- getwd() # NOTE: getwd() can be replaced with any working directory
setwd(my.wd)

# Create subfolders for tables and figures
dir.create(file.path("Tables"))
dir.create(file.path("Tables", "FADN"))
dir.create(file.path("Figures"))
dir.create(file.path("Figures", "FADN"))

# Install missing packages
list.of.packages <- c("dplyr", "broom", "haven", "frontier", "Benchmarking",
                      "dotwhisker", "ggplot2", "ggpubr", "synthpop", 
                      "stargazer")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library(dplyr) # for data management
library(broom) # data preparation
library(haven) # imports .dta files
library(frontier) # efficiency analysis
library(Benchmarking) # DEA analysis
library(dotwhisker) # coefficient plots
library(ggplot2) # plots
library(ggpubr) # to merge plots
library(synthpop) # creates synthetic data
library(stargazer) # for descriptive statistics

# Load FADN data
load("rOutput/FADN.Rds")
summary(lm(log(output)~log(land)+log(labor)+log(material)+log(capstock), data=dat))
summary(frontier::sfa(log(output)~log(land)+log(labor)+log(material)+log(capstock), data=dat))

# Colours
cbPalette <- c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e")

# ------------------------------------- #
#### 1. Generation of synthetic data ####
# ------------------------------------- #

# Set seed for replication
my.seed <- 1234

# Define original data set
dat_orig <- as.data.frame(dat[c("output", "land", "labor", "material", "capstock")])


# ------------------------------------- #
##### 1.a Synthetic data using CART #####
# ------------------------------------- #

# Generate five synthetic data sets (cart)
dat_synt <- syn(dat_orig, method = "cart", m = 5, visit.sequence = 
                  c("output", "land", "labor", "material", "capstock"),
                seed = my.seed,
                smoothing = list(output = "density", 
                                 land = "density", 
                                 labor = "density", 
                                 material = "density", 
                                 capstock = "density"))

# Store synthetic data in data frame
dat_syn1 <- as.data.frame(dat_synt$syn[[1]])
dat_syn2 <- as.data.frame(dat_synt$syn[[2]])
dat_syn3 <- as.data.frame(dat_synt$syn[[3]])
dat_syn4 <- as.data.frame(dat_synt$syn[[4]])
dat_syn5 <- as.data.frame(dat_synt$syn[[5]])

# ----------------------------------------- #
##### 1.b Synthetic data using NORMRANK #####
# ----------------------------------------- #

# Generate five synthetic data sets (normrank)
dat_synt <- syn(dat_orig, method = "parametric", m = 5, visit.sequence = 
                  c("output", "land", "labor", "material", "capstock"),
                seed = my.seed,
                smoothing = list(output = "density", 
                                 land = "density", 
                                 labor = "density", 
                                 material = "density", 
                                 capstock = "density"))

# Store synthetic data in data frame
dat_syn6 <- as.data.frame(dat_synt$syn[[1]])
dat_syn7 <- as.data.frame(dat_synt$syn[[2]])
dat_syn8 <- as.data.frame(dat_synt$syn[[3]])
dat_syn9 <- as.data.frame(dat_synt$syn[[4]])
dat_syn10 <- as.data.frame(dat_synt$syn[[5]])

# ------------------------------- #
#### 2. Descriptive statistics ####
# ------------------------------- #

# Original data
stargazer::stargazer(dat_orig,
                     digits=2,
                     title="Summary statistics for FADN sample (2018), original",
                     nobs = FALSE,
                     omit.summary.stat = c("p25", "p75"),
                     covariate.labels=c("Output", "Land", "Labor",
                                        "Material", "Capital"),
                     type = "text",
                     style = "ajps",
                     out = "Tables/FADN/TabA5a_descr_orig.txt")

# One synthetic data: CART
stargazer::stargazer(dat_syn1,
                     digits=2,
                     title="Summary statistics for FADN sample (2018), CART",
                     nobs = FALSE,
                     omit.summary.stat = c("p25", "p75"),
                     covariate.labels=c("Output", "Land", "Labor",
                                        "Material", "Capital"),
                     type = "text",
                     style = "ajps",
                     out = "Tables/FADN/TabA5b_descr_cart.txt")

# One synthetic data: NORMRANK
stargazer::stargazer(dat_syn6,
                     digits=2,
                     title="Summary statistics for FADN sample (2018), NORMRANK",
                     nobs = FALSE,
                     omit.summary.stat = c("p25", "p75"),
                     covariate.labels=c("Output", "Land", "Labor",
                                        "Material", "Capital"),
                     type = "text",
                     style = "ajps",
                     out = "Tables/FADN/TabA5c_descr_normrank.txt")


# -------------------------------------------- #
#### 3. Estimation of production function   ####
# -------------------------------------------- #

fn.CD <- log(output) ~ log(land) + log(labor) + 
  log(material) + log(capstock)

fn.CD_logs <- loutput ~ lland + llabor + 
  lmaterial + lcapstock

fn.TL <- log(output) ~ log(land/mean(land)) + log(labor/mean(labor)) + log(material/mean(material)) + log(capstock/mean(capstock)) + 
  I(0.5*log(land/mean(land))*log(land/mean(land))) + I(log(land/mean(land))*log(labor/mean(labor))) + I(log(land/mean(land))*log(material/mean(material))) + 
  + I(log(land/mean(land))*log(capstock/mean(capstock))) + 
  I(0.5*log(labor/mean(labor))*log(labor/mean(labor))) + I(log(labor/mean(labor))*log(material/mean(material))) + 
  + I(log(labor/mean(labor))*log(capstock/mean(capstock))) + 
  I(0.5*log(material/mean(material))*log(material/mean(material))) + I(log(material/mean(material))*log(capstock/mean(capstock))) + 
  I(0.5*log(capstock/mean(capstock))*log(capstock/mean(capstock))) 

# --------------------------------------------#
##### 3a Cobb Douglas Production Function #####
# --------------------------------------------#

# Original data
prod.orig <- tidy(lm(fn.CD, data = dat_orig))
prod.orig <- prod.orig %>% mutate(model = "Original data")

# CART: 
  
  # Synthetic data 1
  prod.syn1 <- tidy(lm(fn.CD, data = dat_syn1))
  prod.syn1 <- prod.syn1 %>% mutate(model = "Synthetic data 1")
  
  # Synthetic data 2
  prod.syn2 <- tidy(lm(fn.CD, data = dat_syn2))
  prod.syn2 <- prod.syn2 %>% mutate(model = "Synthetic data 2")
  
  # Synthetic data 3
  prod.syn3 <- tidy(lm(fn.CD, data = dat_syn3))
  prod.syn3 <- prod.syn3 %>% mutate(model = "Synthetic data 3")
  
  # Synthetic data 4
  prod.syn4 <- tidy(lm(fn.CD, data = dat_syn4))
  prod.syn4 <- prod.syn4 %>% mutate(model = "Synthetic data 4")
  
  # Synthetic data 5
  prod.syn5 <- tidy(lm(fn.CD, data = dat_syn5))
  prod.syn5 <- prod.syn5 %>% mutate(model = "Synthetic data 5")

# NORMRANK: 
  
  # Synthetic data 6
  prod.syn6 <- tidy(lm(fn.CD, data = dat_syn6))
  prod.syn6 <- prod.syn6 %>% mutate(model = "Synthetic data 6")
  
  # Synthetic data 7
  prod.syn7 <- tidy(lm(fn.CD, data = dat_syn7))
  prod.syn7 <- prod.syn7 %>% mutate(model = "Synthetic data 7")
  
  # Synthetic data 8
  prod.syn8 <- tidy(lm(fn.CD, data = dat_syn8))
  prod.syn8 <- prod.syn8 %>% mutate(model = "Synthetic data 8")
  
  # Synthetic data 9
  prod.syn9 <- tidy(lm(fn.CD, data = dat_syn9))
  prod.syn9 <- prod.syn9 %>% mutate(model = "Synthetic data 9")
  
  # Synthetic data 10
  prod.syn10 <- tidy(lm(fn.CD, data = dat_syn10))
  prod.syn10 <- prod.syn10 %>% mutate(model = "Synthetic data 10")

# ----------------------------------------------------------------------- #
# Figure A.8. Parameter estimates of the Cobb-Douglas production function #
# with 95% confidence intervals (95% CIs)                                 #
# ----------------------------------------------------------------------- #
  
  # Combine all models
  
    models_cart <- bind_rows(prod.orig, prod.syn1, prod.syn2, prod.syn3,
                             prod.syn4, prod.syn5)
    models_normrank <- bind_rows(prod.orig, prod.syn6, prod.syn7, prod.syn8,
                            prod.syn9, prod.syn10)

  # Make plot for CART
  elast_CD_prod_cart <- dwplot(models_cart, 
                               show_intercept = FALSE,
                               ci = .95,
                               vars_order = c("log(land)", "log(labor)", "log(material)", "log(capstock)"),
                               model_order = c("Original data", 
                                               "Synthetic data 1", 
                                               "Synthetic data 2", 
                                               "Synthetic data 3", 
                                               "Synthetic data 4", 
                                               "Synthetic data 5")
  ) %>% 
    relabel_predictors(
      c("log(land)" = "Land",
        "log(labor)" = "Labor",
        "log(material)" = "Material",
        "log(capstock)" = "Capital")
    ) +
    theme_light(base_family="Times") + 
    xlab("Estimated elasticities with 95% CIs") + ylab("") +
    ggtitle("CART method") +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank()
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(-0.10,1.00)) +
    scale_colour_manual(values=cbPalette)
  
  # Make plot for NORMRANK
  elast_CD_prod_normrank <- dwplot(models_normrank, 
                              show_intercept = FALSE,
                              ci = .95,
                              vars_order = c("log(land)", "log(labor)", "log(material)", "log(capstock)"),
                              model_order = c("Original data", 
                                              "Synthetic data 6", 
                                              "Synthetic data 7", 
                                              "Synthetic data 8", 
                                              "Synthetic data 9", 
                                              "Synthetic data 10")) %>% 
    relabel_predictors(
      c("log(land)" = "Land",
        "log(labor)" = "Labor",
        "log(material)" = "Material",
        "log(capstock)" = "Capital")
    ) +
    theme_light(base_family="Times") + 
    xlab("Estimated elasticities with 95% CIs") + ylab("") +
    ggtitle("NORMRANK method") +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank()
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(-0.10,1.00)) +
    scale_colour_manual(values=cbPalette)
  
  
  # Save plots
  elast_CD_prod <- ggarrange(elast_CD_prod_cart, elast_CD_prod_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/FADN/FigA8_elast_CD_prod.png", elast_CD_prod, width = 23, height = 11, units = "cm")
  
# ----------------------------------------#
##### 3b Translog Production Function #####
# ----------------------------------------#

# Create functions to estimate production function and calculate elasticities

  prod_tl <- function(data) {
    prod <- lm(fn.TL, data = data)
    coef <- prod$coefficients
    el_land <- coef["log(land/mean(land))"] + 
      coef["I(0.5 * log(land/mean(land)) * log(land/mean(land)))"]*log(data$land/mean(data$land)) + 
      coef["I(log(land/mean(land)) * log(labor/mean(labor)))"]*log(data$labor/mean(data$labor)) +
      coef["I(log(land/mean(land)) * log(material/mean(material)))"]*log(data$material/mean(data$material)) +
      coef["I(log(land/mean(land)) * log(capstock/mean(capstock)))"]*log(data$capstock/mean(data$capstock))
    el_labor <- coef["log(labor/mean(labor))"] + 
      coef["I(log(land/mean(land)) * log(labor/mean(labor)))"]*log(data$land/mean(data$land)) + 
      coef["I(0.5 * log(labor/mean(labor)) * log(labor/mean(labor)))"]*log(data$labor/mean(data$labor)) +
      coef["I(log(labor/mean(labor)) * log(material/mean(material)))"]*log(data$material/mean(data$material)) +
      coef["I(log(labor/mean(labor)) * log(capstock/mean(capstock)))"]*log(data$capstock/mean(data$capstock))
    el_material <- coef["log(material/mean(material))"] + 
      coef["I(log(land/mean(land)) * log(material/mean(material)))"]*log(data$land/mean(data$land)) + 
      coef["I(log(labor/mean(labor)) * log(material/mean(material)))"]*log(data$labor/mean(data$labor)) +
      coef["I(0.5 * log(material/mean(material)) * log(material/mean(material)))"]*log(data$material/mean(data$material)) +
      coef["I(log(material/mean(material)) * log(capstock/mean(capstock)))"]*log(data$capstock/mean(data$capstock))
    el_capital <- coef["log(capstock/mean(capstock))"] + 
      coef["I(log(land/mean(land)) * log(capstock/mean(capstock)))"]*log(data$land/mean(data$land)) + 
      coef["I(log(labor/mean(labor)) * log(capstock/mean(capstock)))"]*log(data$labor/mean(data$labor)) +
      coef["I(log(material/mean(material)) * log(capstock/mean(capstock)))"]*log(data$material/mean(data$material)) +
      coef["I(0.5 * log(capstock/mean(capstock)) * log(capstock/mean(capstock)))"]*log(data$capstock/mean(data$capstock))
    result <- list()
    result$elast <- cbind(el_land, el_labor, el_material, el_capital)
    result$model <- tidy(prod)
    return(result)
  }

# Estimate TL prod and calculate elasticities
  
  # Original data
  prod <- prod_tl(dat_orig)
  prod.orig <- prod$model[c(1:5),] %>% 
    mutate(model = "Original data") # adds model label
  dat_orig <- cbind(dat_orig,prod$elast)
  
  # Synthetic data 1
  prod <- prod_tl(dat_syn1)
  prod.syn1 <- prod$model[c(1:5),] %>% 
    mutate(model = "Synthetic data 1") 
  dat_syn1 <- cbind(dat_syn1,prod$elast)
  
  # Synthetic data 2
  prod <- prod_tl(dat_syn2)
  prod.syn2 <- prod$model[c(1:5),] %>% 
    mutate(model = "Synthetic data 2") 
  dat_syn2 <- cbind(dat_syn2,prod$elast)
  
  # Synthetic data 3
  prod <- prod_tl(dat_syn3)
  prod.syn3 <- prod$model[c(1:5),] %>% 
    mutate(model = "Synthetic data 3") 
  dat_syn3 <- cbind(dat_syn3,prod$elast)
  
  # Synthetic data 4
  prod <- prod_tl(dat_syn4)
  prod.syn4 <- prod$model[c(1:5),] %>% 
    mutate(model = "Synthetic data 4") 
  dat_syn4 <- cbind(dat_syn4,prod$elast)
  
  # Synthetic data 5
  prod <- prod_tl(dat_syn5)
  prod.syn5 <- prod$model[c(1:5),] %>%
    mutate(model = "Synthetic data 5") 
  dat_syn5 <- cbind(dat_syn5,prod$elast)
  
  # Synthetic data 6
  prod <- prod_tl(dat_syn6)
  prod.syn6 <- prod$model[c(1:5),] %>% 
    mutate(model = "Synthetic data 6") 
  dat_syn6 <- cbind(dat_syn6,prod$elast)
  
  # Synthetic data 7
  prod <- prod_tl(dat_syn7)
  prod.syn7 <- prod$model[c(1:5),] %>% 
    mutate(model = "Synthetic data 7") 
  dat_syn7 <- cbind(dat_syn7,prod$elast)
  
  # Synthetic data 8
  prod <- prod_tl(dat_syn8)
  prod.syn8 <- prod$model[c(1:5),] %>% 
    mutate(model = "Synthetic data 8") 
  dat_syn8 <- cbind(dat_syn8,prod$elast)
  
  # Synthetic data 9
  prod <- prod_tl(dat_syn9)
  prod.syn9 <- prod$model[c(1:5),] %>% 
    mutate(model = "Synthetic data 9") 
  dat_syn9 <- cbind(dat_syn9,prod$elast)
  
  # Synthetic data 10
  prod <- prod_tl(dat_syn10)
  prod.syn10 <- prod$model[c(1:5),] %>% 
    mutate(model = "Synthetic data 10") 
  dat_syn10 <- cbind(dat_syn10,prod$elast)

#--------------#
# Figure A.11. #
#--------------#

  # Combine all models
      
    models_cart <- bind_rows(prod.orig, prod.syn1, prod.syn2, prod.syn3,
                             prod.syn4, prod.syn5)
    models_normrank <- bind_rows(prod.orig, prod.syn6, prod.syn7, prod.syn8,
                            prod.syn9, prod.syn10)
  
  # Make plot for CART
  elast_TL_prod_cart <- dwplot(models_cart, 
                               show_intercept = FALSE,
                               ci = .95,
                               vars_order = c("log(land/mean(land))", "log(labor/mean(labor))", "log(material/mean(material))", "log(capstock/mean(capstock))"),
                               model_order = c("Original data", 
                                               "Synthetic data 1", 
                                               "Synthetic data 2", 
                                               "Synthetic data 3", 
                                               "Synthetic data 4", 
                                               "Synthetic data 5")) %>% 
    relabel_predictors(
      c("log(land/mean(land))" = "Land",
        "log(labor/mean(labor))" = "Labor",
        "log(material/mean(material))" = "Material",
        "log(capstock/mean(capstock))" = "Capital")
    ) +
    theme_light(base_family="Times") + 
    xlab("Estimated elasticities with 95% CIs") + ylab("") +
    ggtitle("CART method") +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank()
    ) +
    scale_x_continuous(breaks = seq(0, 0.75, by = 0.25), limits=c(-0.06,0.95)) +
    scale_colour_manual(values=cbPalette)
  
  # Make plot for NORMRANK
  elast_TL_prod_normrank <- dwplot(models_normrank, 
                              show_intercept = FALSE,
                              ci = .95,
                              vars_order = c("log(land/mean(land))", "log(labor/mean(labor))", "log(material/mean(material))", "log(capstock/mean(capstock))"),
                              model_order = c("Original data", 
                                              "Synthetic data 6", 
                                              "Synthetic data 7", 
                                              "Synthetic data 8", 
                                              "Synthetic data 9", 
                                              "Synthetic data 10")) %>% 
    relabel_predictors(
      c("log(land/mean(land))" = "Land",
        "log(labor/mean(labor))" = "Labor",
        "log(material/mean(material))" = "Material",
        "log(capstock/mean(capstock))" = "Capital")
    ) +
    theme_light(base_family="Times") + 
    xlab("Estimated elasticities with 95% CIs") + ylab("") +
    ggtitle("NORMRANK method") +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank()
    ) +
    scale_x_continuous(breaks = seq(0, 0.75, by = 0.25), limits=c(-0.06,0.95)) +
    scale_colour_manual(values=cbPalette) 
  
  # Save plots
  elast_TL_prod <- ggarrange(elast_TL_prod_cart, elast_TL_prod_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/FADN/FigA11_elast_TL_prod.png", elast_TL_prod, width = 23, height = 11, units = "cm")


#--------------#
# Figure A.12. #
#--------------#
  
  # Make data frame with elasticities from all datasets (original and synthetic)
    
    # Land
    nr_obs <- length(dat_orig$el_land)
    el_land <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_land) <- c("el", "data") 
    el_land$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                      rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                      rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                      rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                      rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                      rep("Synthetic data 10", nr_obs))
    el_land$el <- c(dat_orig$el_land,dat_syn1$el_land,dat_syn2$el_land,dat_syn3$el_land,dat_syn4$el_land,dat_syn5$el_land,
                    dat_syn6$el_land,dat_syn7$el_land,dat_syn8$el_land,dat_syn9$el_land,dat_syn10$el_land)
    el_land$inp <- "Land"
    
    # Labor
    el_labor <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_labor) <- c("el", "data") 
    el_labor$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                       rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                       rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                       rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                       rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                       rep("Synthetic data 10", nr_obs))
    el_labor$el <- c(dat_orig$el_labor,dat_syn1$el_labor,dat_syn2$el_labor,dat_syn3$el_labor,dat_syn4$el_labor,dat_syn5$el_labor,
                     dat_syn6$el_labor,dat_syn7$el_labor,dat_syn8$el_labor,dat_syn9$el_labor,dat_syn10$el_labor)
    el_labor$inp <- "Labor"
    
    # Material
    el_material <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_material) <- c("el", "data") 
    el_material$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                          rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                          rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                          rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                          rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                          rep("Synthetic data 10", nr_obs))
    el_material$el <- c(dat_orig$el_material,dat_syn1$el_material,dat_syn2$el_material,dat_syn3$el_material,dat_syn4$el_material,dat_syn5$el_material,
                        dat_syn6$el_material,dat_syn7$el_material,dat_syn8$el_material,dat_syn9$el_material,dat_syn10$el_material)
    el_material$inp <- "Material"
    
    # Capital
    el_capital <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_capital) <- c("el", "data") 
    el_capital$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                         rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                         rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                         rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                         rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                         rep("Synthetic data 10", nr_obs))
    el_capital$el <- c(dat_orig$el_capital,dat_syn1$el_capital,dat_syn2$el_capital,dat_syn3$el_capital,dat_syn4$el_capital,dat_syn5$el_capital,
                       dat_syn6$el_capital,dat_syn7$el_capital,dat_syn8$el_capital,dat_syn9$el_capital,dat_syn10$el_capital)
    el_capital$inp <- "Capital"
    
    # All
    el_prod <- rbind(el_land,el_labor,el_material,el_capital)
    el_prod_cart <- el_prod %>% 
      filter(data == "Original data" | data == "Synthetic data 1" |
               data == "Synthetic data 2" | data == "Synthetic data 3" |
               data == "Synthetic data 4" | data == "Synthetic data 5")
    el_prod_normrank <- el_prod %>% 
      filter(data == "Original data" | data == "Synthetic data 6" |
               data == "Synthetic data 7" | data == "Synthetic data 8" |
               data == "Synthetic data 9" | data == "Synthetic data 10")
    
  # Make sure data are plotted in right order 
    
    el_prod_cart$inp <- factor(el_prod_cart$inp,
                               levels=c("Land","Labor","Material","Capital"))
    el_prod_cart$data <- factor(el_prod_cart$data,
                                levels=c("Synthetic data 10",
                                         "Synthetic data 9",
                                         "Synthetic data 8",
                                         "Synthetic data 7",
                                         "Synthetic data 6",
                                         "Synthetic data 5",
                                         "Synthetic data 4",
                                         "Synthetic data 3",
                                         "Synthetic data 2",
                                         "Synthetic data 1",
                                         "Original data"))
    el_prod_normrank$inp <- factor(el_prod_normrank$inp,
                              levels=c("Land","Labor","Material","Capital"))
    el_prod_normrank$data <- factor(el_prod_normrank$data,
                               levels=c("Synthetic data 10",
                                        "Synthetic data 9",
                                        "Synthetic data 8",
                                        "Synthetic data 7",
                                        "Synthetic data 6",
                                        "Synthetic data 5",
                                        "Synthetic data 4",
                                        "Synthetic data 3",
                                        "Synthetic data 2",
                                        "Synthetic data 1",
                                        "Original data"))

  # Make plot for CART
  elast_TL_fl_prod_cart <- ggplot(el_prod_cart, aes(x=data, y=el, fill=data)) +
    geom_violin(size=0.2) +
    theme_light() +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position="bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      axis.text.y=element_blank() 
    ) +
    coord_flip() + # To switch X and Y axes 
    xlab("") +
    ylab("Elasticity") +
    ggtitle("CART method") +
    facet_wrap(~ inp) +
    scale_y_continuous(breaks = seq(-0.3, 1.2, by = 0.3), limits=c(-0.4,1.3)) +
    scale_fill_manual(values=cbPalette) 
  
  # Make plot for NORMRANK
  elast_TL_fl_prod_normrank <- ggplot(el_prod_normrank, aes(x=data, y=el, fill=data)) +
    geom_violin(size=0.2) +
    theme_light() +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position="bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      axis.text.y=element_blank() 
    ) +
    coord_flip() + # To switch X and Y axes 
    xlab("") +
    ylab("Elasticity") +
    ggtitle("NORMRANK method") +
    facet_wrap(~ inp) +
    scale_y_continuous(breaks = seq(-0.3, 1.2, by = 0.3), limits=c(-0.4,1.3)) +
    scale_fill_manual(values=cbPalette)         
  
  # Save plots
  elast_TL_fl_prod <- ggarrange(elast_TL_fl_prod_cart, elast_TL_fl_prod_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/FADN/FigA12_elast_TL_fl_prod.png", elast_TL_fl_prod, width = 23, height = 11, units = "cm")

# -------------------------------------------- #
#### 4. Estimation of production frontier   ####
# -------------------------------------------- #
  
# --------------------------------------------- #
##### 4.a Cobb-Douglas Production Frontier  #####
# --------------------------------------------- #

  # Create functions to estimate production frontier
  
    front_cd <- function(data) {
      front <- frontier::sfa(fn.CD, data = data)
      result <- list()
      # Store TE scores
      result$te <- frontier::efficiencies(front, asInData=TRUE)
      # prepare for coefficient plot (dwplot)
      front <- summary(front)
      front <- as_tibble(front$mleParam)
      front$term <- c("Intercept", "Land", "Labor", "Material", "Capital", 
                      "Sigma2", "Gamma")
      front <- front %>% 
        select(term, everything()) #re-order variables
      names(front) <- names(prod.orig)[c(1:5)]
      result$estimates <- front
      return(result)
    }
    
    # Original data
    front <- front_cd(dat_orig)
    front.orig <- front$estimates
    front.orig$model <- "Original data"
    te_orig <- front$te
    
    # Synthetic data 1
    front <- front_cd(dat_syn1)
    front.syn1 <- front$estimates
    front.syn1$model <- "Synthetic data 1"
    te_syn1 <- front$te
    
    # Synthetic data 2
    front <- front_cd(dat_syn2)
    front.syn2 <- front$estimates
    front.syn2$model <- "Synthetic data 2"
    te_syn2 <- front$te
    
    # Synthetic data 3
    front <- front_cd(dat_syn3)
    front.syn3 <- front$estimates
    front.syn3$model <- "Synthetic data 3"
    te_syn3 <- front$te
    
    # Synthetic data 4
    front <- front_cd(dat_syn4)
    front.syn4 <- front$estimates
    front.syn4$model <- "Synthetic data 4"
    te_syn4 <- front$te
    
    # Synthetic data 5
    front <- front_cd(dat_syn5)
    front.syn5 <- front$estimates
    front.syn5$model <- "Synthetic data 5"
    te_syn5 <- front$te
    
    # Synthetic data 6
    front <- front_cd(dat_syn6)
    front.syn6 <- front$estimates
    front.syn6$model <- "Synthetic data 6"
    te_syn6 <- front$te
    
    # Synthetic data 7
    front <- front_cd(dat_syn7)
    front.syn7 <- front$estimates
    front.syn7$model <- "Synthetic data 7"
    te_syn7 <- front$te
    
    # Synthetic data 8
    front <- front_cd(dat_syn8)
    front.syn8 <- front$estimates
    front.syn8$model <- "Synthetic data 8"
    te_syn8 <- front$te
    
    # Synthetic data 9
    front <- front_cd(dat_syn9)
    front.syn9 <- front$estimates
    front.syn9$model <- "Synthetic data 9"
    te_syn9 <- front$te
    
    # Synthetic data 10
    front <- front_cd(dat_syn10)
    front.syn10 <- front$estimates
    front.syn10$model <- "Synthetic data 10"
    te_syn10 <- front$te
    
#--------------#
# Figure A.13. #
#--------------#

  # Combine all models
    
    models_cart <- bind_rows(front.orig, front.syn1, front.syn2, front.syn3,
                             front.syn4, front.syn5)
    models_normrank <- bind_rows(front.orig, front.syn6, front.syn7, front.syn8,
                            front.syn9, front.syn10)
  
  # Make plot for CART
  elast_CD_front_cart <- dwplot(models_cart, 
                                show_intercept = FALSE,
                                ci = .95,
                                vars_order = c("Land", "Labor", "Material", "Capital"),
                                model_order = c("Original data", 
                                                "Synthetic data 1", 
                                                "Synthetic data 2", 
                                                "Synthetic data 3", 
                                                "Synthetic data 4", 
                                                "Synthetic data 5")) +
    theme_light(base_family="Times") + 
    xlab("Estimated elasticities with 95% CIs") + 
    ylab("") +
    ggtitle("CART method") +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(-0.12,1)) +
    scale_colour_manual(values=cbPalette)
  
  # Make plot for NORMRANK
  elast_CD_front_normrank <- dwplot(models_normrank, 
                               show_intercept = FALSE,
                               ci = .95,
                               vars_order = c("Land", "Labor", "Material", "Capital"),
                               model_order = c("Original data", 
                                               "Synthetic data 6", 
                                               "Synthetic data 7", 
                                               "Synthetic data 8", 
                                               "Synthetic data 9", 
                                               "Synthetic data 10")) +
    theme_light(base_family="Times") + 
    xlab("Estimated elasticities with 95% CIs") + 
    ylab("") +
    ggtitle("NORMRANK method") +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(-0.12,1)) +
    scale_colour_manual(values=cbPalette)
  
  # Save plots
  elast_CD_front <- ggarrange(elast_CD_front_cart, elast_CD_front_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/FADN/FigA13_elast_CD_front.png", elast_CD_front, width = 23, height = 11, units = "cm")
  
  # ---------- #
  # Table A.6. #
  # ---------- #

  # Combine all TE scores
  descr_te_cd <- cbind(te_orig,te_syn1,te_syn2,te_syn3,te_syn4,te_syn5,
                     te_syn6,te_syn7,te_syn8,te_syn9,te_syn10)
  
  # Create table with descriptive statistics
  stargazer(as.data.frame(descr_te_cd),
            digits=2,
            summary = TRUE,
            iqr = TRUE,
            title="Summary statistics: Technical efficiency scores, Cobb-Douglas, half-normal distribution",
            nobs = FALSE,
            covariate.labels=c("Original data", "Synthetic data 1", 
                               "Synthetic data 2", "Synthetic data 3",
                               "Synthetic data 4","Synthetic data 5",
                               "Synthetic data 6","Synthetic data 7",
                               "Synthetic data 8","Synthetic data 9",
                               "Synthetic data 10","Synthetic data 11",
                               "Synthetic data 12","Synthetic data 13",
                               "Synthetic data 14","Synthetic data 15"),
            type = "text",
            style = "ajps",
            out = "Tables/FADN/TabA6_descr_te_cd.txt")

#-------------#
# Figure A.9. #
#-------------#

  # Make data frame with technical efficiency scores from all datasets
  
  
  # Original and CART data
  nr_obs <- length(te_orig)
  efficiencies_cart <- data.frame(matrix(ncol = 2, nrow = 6*nr_obs))
  colnames(efficiencies_cart) <- c("TE", "data") 
  efficiencies_cart$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                              rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                              rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs))
  efficiencies_cart$TE <- c(te_orig,te_syn1,te_syn2,te_syn3,te_syn4,te_syn5)
  
  # Original and NORMRANK data
  efficiencies_normrank <- data.frame(matrix(ncol = 2, nrow = 6*nr_obs))
  colnames(efficiencies_normrank) <- c("TE", "data") 
  efficiencies_normrank$data <- c(rep("Original data", nr_obs),rep("Synthetic data 6", nr_obs),
                             rep("Synthetic data 7", nr_obs),rep("Synthetic data 8", nr_obs),
                             rep("Synthetic data 9", nr_obs),rep("Synthetic data 10", nr_obs))
  efficiencies_normrank$TE <- c(te_orig,te_syn6,te_syn7,te_syn8,te_syn9,te_syn10)
  
  # Make sure data are plotted in right order 
  
    efficiencies_cart$data <- factor(efficiencies_cart$data,
                                     levels=c("Synthetic data 5",
                                              "Synthetic data 4",
                                              "Synthetic data 3",
                                              "Synthetic data 2",
                                              "Synthetic data 1",
                                              "Original data"))
    efficiencies_normrank$data <- factor(efficiencies_normrank$data,
                                    levels=c("Synthetic data 10",
                                             "Synthetic data 9",
                                             "Synthetic data 8",
                                             "Synthetic data 7",
                                             "Synthetic data 6",
                                             "Original data"))

  # Make plot for CART
  te_cd_cart <-     
    efficiencies_cart %>%
    ggplot(aes(x = TE , y = data,  fill = data)) +
    ggridges::geom_density_ridges(scale = 1, 
                                  stat = "binline", 
                                  binwidth = 0.01) +
    theme_light(base_family="Times") + 
    xlab("Technical efficiency") + ylab("") +
    ggtitle("CART method") +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      plot.title = element_text(size=10, hjust = 0.5),
      axis.text.y=element_blank() 
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(0,1)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 1", 
                                 "Synthetic data 2", 
                                 "Synthetic data 3", 
                                 "Synthetic data 4", 
                                 "Synthetic data 5"))
    
  # Make plot for NORMRANK
  te_cd_normrank <-     
    efficiencies_normrank %>%
    ggplot(aes(x = TE , y = data,  fill = data)) +
    ggridges::geom_density_ridges(scale = 1, 
                                  stat = "binline", 
                                  binwidth = 0.01) +
    theme_light(base_family="Times") + 
    xlab("Technical efficiency") + ylab("") +
    ggtitle("NORMRANK method") +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      plot.title = element_text(size=10, hjust = 0.5),
      axis.text.y=element_blank() 
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(0,1)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 6", 
                                 "Synthetic data 7", 
                                 "Synthetic data 8", 
                                 "Synthetic data 9", 
                                 "Synthetic data 10")) 


  # Save plots
  te_cd <- ggarrange(te_cd_cart, te_cd_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/FADN/FigA9_te_cd.png", te_cd, width = 23, height = 11, units = "cm")

# ----------------------------------------- #
##### 4.b Translog Production Frontier  #####
# ----------------------------------------- # 

# Create functions to estimate production frontiers and calculate elasticities

  front_tl <- function(data) {
    front <- frontier::sfa(fn.TL, data = data)
    result <- list()
    # Store TE scores
    result$te <- frontier::efficiencies(front, asInData=TRUE) 
    #calculate elasticities
    coef <- front$mleParam
    el_land_front <- coef["log(land/mean(land))"] + 
      coef["I(0.5 * log(land/mean(land)) * log(land/mean(land)))"]*log(data$land/mean(data$land)) + 
      coef["I(log(land/mean(land)) * log(labor/mean(labor)))"]*log(data$labor/mean(data$labor)) +
      coef["I(log(land/mean(land)) * log(material/mean(material)))"]*log(data$material/mean(data$material)) +
      coef["I(log(land/mean(land)) * log(capstock/mean(capstock)))"]*log(data$capstock/mean(data$capstock))
    el_labor_front <- coef["log(labor/mean(labor))"] + 
      coef["I(log(land/mean(land)) * log(labor/mean(labor)))"]*log(data$land/mean(data$land)) + 
      coef["I(0.5 * log(labor/mean(labor)) * log(labor/mean(labor)))"]*log(data$labor/mean(data$labor)) +
      coef["I(log(labor/mean(labor)) * log(material/mean(material)))"]*log(data$material/mean(data$material)) +
      coef["I(log(labor/mean(labor)) * log(capstock/mean(capstock)))"]*log(data$capstock/mean(data$capstock))
    el_material_front <- coef["log(material/mean(material))"] + 
      coef["I(log(land/mean(land)) * log(material/mean(material)))"]*log(data$land/mean(data$land)) + 
      coef["I(log(labor/mean(labor)) * log(material/mean(material)))"]*log(data$labor/mean(data$labor)) +
      coef["I(0.5 * log(material/mean(material)) * log(material/mean(material)))"]*log(data$material/mean(data$material)) +
      coef["I(log(material/mean(material)) * log(capstock/mean(capstock)))"]*log(data$capstock/mean(data$capstock))
    el_capital_front <- coef["log(capstock/mean(capstock))"] + 
      coef["I(log(land/mean(land)) * log(capstock/mean(capstock)))"]*log(data$land/mean(data$land)) + 
      coef["I(log(labor/mean(labor)) * log(capstock/mean(capstock)))"]*log(data$labor/mean(data$labor)) +
      coef["I(log(material/mean(material)) * log(capstock/mean(capstock)))"]*log(data$material/mean(data$material)) +
      coef["I(0.5 * log(capstock/mean(capstock)) * log(capstock/mean(capstock)))"]*log(data$capstock/mean(data$capstock))
    # prepare for coefficient plot (dwplot)
    front <- summary(front)
    front <- as_tibble(front$mleParam)
    front <- front[c(1:5),]
    front$term <- c("Intercept", "Land", "Labor", "Material", "Capital")
    front <- front %>% 
      select(term, everything()) #re-order variables
    names(front) <- names(prod.orig)[c(1:5)]
    result$estimates <- front
    result$elast <- cbind(el_land_front, el_labor_front, el_material_front, 
                          el_capital_front)
    return(result)
  }
  
  # Estimate TL prod and calculate elasticities

    # Original data
    front <- front_tl(dat_orig)
    front.orig <- front$estimates %>% 
      mutate(model = "Original data") # adds model label
    te_orig <- front$te
    dat_orig <- cbind(dat_orig,front$elast)
    
    # Synthetic data 1
    front <- front_tl(dat_syn1)
    front.syn1 <- front$estimates %>% 
      mutate(model = "Synthetic data 1") 
    te_syn1 <- front$te
    dat_syn1 <- cbind(dat_syn1,front$elast)  
  
    # Synthetic data 2
    front <- front_tl(dat_syn2)
    front.syn2 <- front$estimates %>% 
      mutate(model = "Synthetic data 2") 
    te_syn2 <- front$te
    dat_syn2 <- cbind(dat_syn2,front$elast)
    
    # Synthetic data 3
    front <- front_tl(dat_syn3)
    front.syn3 <- front$estimates %>% 
      mutate(model = "Synthetic data 3") 
    te_syn3 <- front$te
    dat_syn3 <- cbind(dat_syn3,front$elast)
    
    # Synthetic data 4
    front <- front_tl(dat_syn4)
    front.syn4 <- front$estimates %>% 
      mutate(model = "Synthetic data 4") 
    te_syn4 <- front$te
    dat_syn4 <- cbind(dat_syn4,front$elast)
    
    # Synthetic data 5
    front <- front_tl(dat_syn5)
    front.syn5 <- front$estimates %>% 
      mutate(model = "Synthetic data 5") 
    te_syn5 <- front$te
    dat_syn5 <- cbind(dat_syn5,front$elast)
    
    # Synthetic data 6
    front <- front_tl(dat_syn6)
    front.syn6 <- front$estimates %>% 
      mutate(model = "Synthetic data 6") 
    te_syn6 <- front$te
    dat_syn6 <- cbind(dat_syn6,front$elast)
    
    # Synthetic data 7
    front <- front_tl(dat_syn7)
    front.syn7 <- front$estimates %>% 
      mutate(model = "Synthetic data 7") 
    te_syn7 <- front$te
    dat_syn7 <- cbind(dat_syn7,front$elast)
    
    # Synthetic data 8
    front <- front_tl(dat_syn8)
    front.syn8 <- front$estimates %>% 
      mutate(model = "Synthetic data 8") 
    te_syn8 <- front$te
    dat_syn8 <- cbind(dat_syn8,front$elast)
    
    # Synthetic data 9
    front <- front_tl(dat_syn9)
    front.syn9 <- front$estimates %>% 
      mutate(model = "Synthetic data 9") 
    te_syn9 <- front$te
    dat_syn9 <- cbind(dat_syn9,front$elast)
    
    # Synthetic data 10
    front <- front_tl(dat_syn10)
    front.syn10 <- front$estimates %>% 
      mutate(model = "Synthetic data 10") 
    te_syn10 <- front$te
    dat_syn10 <- cbind(dat_syn10,front$elast)

#--------------#
# Figure A.14. #
#--------------#

  # Combine all models
  models_cart <- bind_rows(front.orig, front.syn1, front.syn2, front.syn3,
                         front.syn4, front.syn5)

  models_normrank <- bind_rows(front.orig, front.syn6, front.syn7, front.syn8,
                        front.syn9, front.syn10)

  # Make plot for CART
  elast_TL_front_cart <- dwplot(models_cart, 
                                show_intercept = FALSE,
                                ci = .95,
                                vars_order = c("Land", "Labor", "Material", "Capital"),
                                model_order = c("Original data", 
                                                "Synthetic data 1", 
                                                "Synthetic data 2", 
                                                "Synthetic data 3", 
                                                "Synthetic data 4", 
                                                "Synthetic data 5")
  ) +
    theme_light(base_family="Times") + 
    xlab("Estimated elasticities with 95% CIs") + ylab("") +
    ggtitle("CART method") +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank() 
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(-0.07,1)) +
    scale_colour_manual(values=cbPalette)
  
  # Make plot for NORMRANK
  elast_TL_front_normrank <- dwplot(models_normrank, 
                               show_intercept = FALSE,
                               ci = .95,
                               vars_order = c("Land", "Labor", "Material", "Capital"),
                               model_order = c("Original data", 
                                               "Synthetic data 6", 
                                               "Synthetic data 7", 
                                               "Synthetic data 8", 
                                               "Synthetic data 9", 
                                               "Synthetic data 10")
  ) +
    theme_light(base_family="Times") + 
    xlab("Estimated elasticities with 95% CIs") + ylab("") +
    ggtitle("NORMRANK method") +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank() 
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(-0.07,1)) +
    scale_colour_manual(values=cbPalette)

  # Save plots
  elast_TL_front <- ggarrange(elast_TL_front_cart, elast_TL_front_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/FADN/FigA14_elast_TL_front.png", elast_TL_front, width = 23, height = 11, units = "cm")
  

#--------------#
# Figure A.15. #
#--------------#

  # Make data frame with elasticities from all datasets (original and synthetic)

    # Land
    nr_obs <- length(dat_orig$el_land)
    el_land <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_land) <- c("el", "data") 
    el_land$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                      rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                      rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                      rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                      rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                      rep("Synthetic data 10", nr_obs))
    el_land$el <- c(dat_orig$el_land,dat_syn1$el_land,dat_syn2$el_land,dat_syn3$el_land,dat_syn4$el_land,dat_syn5$el_land,
                    dat_syn6$el_land,dat_syn7$el_land,dat_syn8$el_land,dat_syn9$el_land,dat_syn10$el_land)
    el_land$inp <- "Land"
    
    # Labor
    el_labor <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_labor) <- c("el", "data") 
    el_labor$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                       rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                       rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                       rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                       rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                       rep("Synthetic data 10", nr_obs))
    el_labor$el <- c(dat_orig$el_labor,dat_syn1$el_labor,dat_syn2$el_labor,dat_syn3$el_labor,dat_syn4$el_labor,dat_syn5$el_labor,
                     dat_syn6$el_labor,dat_syn7$el_labor,dat_syn8$el_labor,dat_syn9$el_labor,dat_syn10$el_labor)
    el_labor$inp <- "Labor"
    
    # Material
    el_material <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_material) <- c("el", "data") 
    el_material$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                          rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                          rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                          rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                          rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                          rep("Synthetic data 10", nr_obs))
    el_material$el <- c(dat_orig$el_material,dat_syn1$el_material,dat_syn2$el_material,dat_syn3$el_material,dat_syn4$el_material,dat_syn5$el_material,
                        dat_syn6$el_material,dat_syn7$el_material,dat_syn8$el_material,dat_syn9$el_material,dat_syn10$el_material)
    el_material$inp <- "Material"
    
    # Capital
    el_capital <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_capital) <- c("el", "data") 
    el_capital$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                         rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                         rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                         rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                         rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                         rep("Synthetic data 10", nr_obs))
    el_capital$el <- c(dat_orig$el_capital,dat_syn1$el_capital,dat_syn2$el_capital,dat_syn3$el_capital,dat_syn4$el_capital,dat_syn5$el_capital,
                       dat_syn6$el_capital,dat_syn7$el_capital,dat_syn8$el_capital,dat_syn9$el_capital,dat_syn10$el_capital)
    el_capital$inp <- "Capital"
    
    # All
    el_front <- rbind(el_land,el_labor,el_material,el_capital)
    el_front_cart <- el_front %>% 
      filter(data == "Original data" | data == "Synthetic data 1" |
               data == "Synthetic data 2" | data == "Synthetic data 3" |
               data == "Synthetic data 4" | data == "Synthetic data 5")
    el_front_normrank <- el_front %>% 
      filter(data == "Original data" | data == "Synthetic data 6" |
               data == "Synthetic data 7" | data == "Synthetic data 8" |
               data == "Synthetic data 9" | data == "Synthetic data 10")

  # Make sure data are plotted in right order 
    
    el_front_cart$inp <- factor(el_front_cart$inp,
                                levels=c("Land","Labor","Material","Capital"))
    el_front_cart$data <- factor(el_front_cart$data,
                                 levels=c("Synthetic data 5",
                                          "Synthetic data 4",
                                          "Synthetic data 3",
                                          "Synthetic data 2",
                                          "Synthetic data 1",
                                          "Original data"))
    el_front_normrank$inp <- factor(el_front_normrank$inp,
                               levels=c("Land","Labor","Material","Capital"))
    el_front_normrank$data <- factor(el_front_normrank$data,
                                levels=c("Synthetic data 10",
                                         "Synthetic data 9",
                                         "Synthetic data 8",
                                         "Synthetic data 7",
                                         "Synthetic data 6",
                                         "Original data"))
    
  # Make plot for CART
  elast_TL_fl_front_cart <- ggplot(el_front_cart, aes(x=data, y=el, fill=data)) +
    geom_violin(size=0.2) +
    theme_light() +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position="bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      axis.text.y=element_blank() 
    ) +
    coord_flip() + # To switch X and Y axes 
    xlab("") +
    ylab("Elasticity") +
    ggtitle("CART method") +
    facet_wrap(~ inp) +
    scale_y_continuous(breaks = seq(-0.3, 1.2, by = 0.3), limits=c(-0.4,1.3)) +
    scale_fill_manual(values=cbPalette) 
  
  # Make plot for NORMRANK
  elast_TL_fl_front_normrank <- ggplot(el_front_normrank, aes(x=data, y=el, fill=data)) +
    geom_violin(size=0.2) +
    theme_light() +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      legend.position="bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      axis.text.y=element_blank() 
    ) +
    coord_flip() + # To switch X and Y axes 
    xlab("") +
    ylab("Elasticity") +
    ggtitle("NORMRANK method") +
    facet_wrap(~ inp) +
    scale_y_continuous(breaks = seq(-0.3, 1.2, by = 0.3), limits=c(-0.4,1.3)) +
    scale_fill_manual(values=cbPalette)         
  
  # Save plots
  elast_TL_fl_front <- ggarrange(elast_TL_fl_front_cart, elast_TL_fl_front_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/FADN/FigA15_elast_TL_fl_front.png", elast_TL_fl_front, width = 23, height = 11, units = "cm")

  # ---------- #
  # Table A.7. #
  # ---------- # 
  
  # Combine all TE scores
  descr_te_tl <- cbind(te_orig,te_syn1,te_syn2,te_syn3,te_syn4,te_syn5,
                      te_syn6,te_syn7,te_syn8,te_syn9,te_syn10)

  # Create table with descriptive statistics
  stargazer(as.data.frame(descr_te_tl),
            digits=2,
            summary = TRUE,
            iqr = TRUE,
            title="Summary statistics: Technical efficiency scores, Translog",
            nobs = FALSE,
            covariate.labels=c("Original data", "Synthetic data 1", 
                               "Synthetic data 2", "Synthetic data 3",
                               "Synthetic data 4","Synthetic data 5",
                               "Synthetic data 6","Synthetic data 7",
                               "Synthetic data 8","Synthetic data 9",
                               "Synthetic data 10"),
            type = "text",
            style = "ajps",
            out = "Tables/FADN/TabA7_descr_te_tl.txt")

  
#--------------#
# Figure A.16. #
#--------------#
  
  # Make data frame with technical efficiency scores from all datasets
  
    # Original and CART data
    nr_obs <- length(te_orig)
    efficiencies_cart <- data.frame(matrix(ncol = 2, nrow = 6*nr_obs))
    colnames(efficiencies_cart) <- c("TE", "data") 
    efficiencies_cart$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                                rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                                rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs))
    efficiencies_cart$TE <- c(te_orig,te_syn1,te_syn2,te_syn3,te_syn4,te_syn5)
    
    # Original and NORMRANK data
    efficiencies_normrank <- data.frame(matrix(ncol = 2, nrow = 6*nr_obs))
    colnames(efficiencies_normrank) <- c("TE", "data") 
    efficiencies_normrank$data <- c(rep("Original data", nr_obs),rep("Synthetic data 6", nr_obs),
                               rep("Synthetic data 7", nr_obs),rep("Synthetic data 8", nr_obs),
                               rep("Synthetic data 9", nr_obs),rep("Synthetic data 10", nr_obs))
    efficiencies_normrank$TE <- c(te_orig,te_syn6,te_syn7,te_syn8,te_syn9,te_syn10)

  # Make sure data are plotted in right order
    
    efficiencies_cart$data <- factor(efficiencies_cart$data,
                                     levels=c("Synthetic data 5",
                                              "Synthetic data 4",
                                              "Synthetic data 3",
                                              "Synthetic data 2",
                                              "Synthetic data 1",
                                              "Original data"))
    efficiencies_normrank$data <- factor(efficiencies_normrank$data,
                                    levels=c("Synthetic data 10",
                                             "Synthetic data 9",
                                             "Synthetic data 8",
                                             "Synthetic data 7",
                                             "Synthetic data 6",
                                             "Original data"))

  # Make plot for CART
  te_tl_cart <-     
    efficiencies_cart %>%
    ggplot(aes(x = TE , y = data,  fill = data)) +
    ggridges::geom_density_ridges(scale = 1, 
                                  stat = "binline", 
                                  binwidth = 0.01) +
    theme_light(base_family="Times") + 
    xlab("Technical efficiency") + ylab("") +
    ggtitle("CART method") +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      plot.title = element_text(size=10, hjust = 0.5),
      axis.text.y=element_blank() 
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(0,1)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 1", 
                                 "Synthetic data 2", 
                                 "Synthetic data 3", 
                                 "Synthetic data 4", 
                                 "Synthetic data 5"))
  
  # Make plot for NORMRANK
  te_tl_normrank <-     
    efficiencies_normrank %>%
    ggplot(aes(x = TE , y = data,  fill = data)) +
    ggridges::geom_density_ridges(scale = 1, 
                                  stat = "binline", 
                                  binwidth = 0.01) +
    theme_light(base_family="Times") + 
    xlab("Technical efficiency") + ylab("") +
    ggtitle("NORMRANK method") +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      plot.title = element_text(size=10, hjust = 0.5),
      axis.text.y=element_blank() 
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(0,1)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 6", 
                                 "Synthetic data 7", 
                                 "Synthetic data 8", 
                                 "Synthetic data 9", 
                                 "Synthetic data 10")) 

  # Save plots
  te_tl <- ggarrange(te_tl_cart, te_tl_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/FADN/FigA16_te_tl.png", te_tl, width = 23, height = 11, units = "cm")


# ------------------------------------ #
#### 5. Data envelopment analysis   ####
# ------------------------------------ #

  # Define inputs and outputs

    inputs_orig <- data.frame(dat[c("land","labor","material","capstock")])
    outputs_orig <- data.frame(dat[c("output")])
    inputs_syn1 <- data.frame(dat_syn1[c("land","labor","material","capstock")])
    outputs_syn1 <- data.frame(dat_syn1[c("output")]) 
    inputs_syn2 <- data.frame(dat_syn2[c("land","labor","material","capstock")])
    outputs_syn2 <- data.frame(dat_syn2[c("output")]) 
    inputs_syn3 <- data.frame(dat_syn3[c("land","labor","material","capstock")])
    outputs_syn3 <- data.frame(dat_syn3[c("output")]) 
    inputs_syn4 <- data.frame(dat_syn4[c("land","labor","material","capstock")])
    outputs_syn4 <- data.frame(dat_syn4[c("output")]) 
    inputs_syn5 <- data.frame(dat_syn5[c("land","labor","material","capstock")])
    outputs_syn5 <- data.frame(dat_syn5[c("output")]) 
    inputs_syn6 <- data.frame(dat_syn6[c("land","labor","material","capstock")])
    outputs_syn6 <- data.frame(dat_syn6[c("output")]) 
    inputs_syn7 <- data.frame(dat_syn7[c("land","labor","material","capstock")])
    outputs_syn7 <- data.frame(dat_syn7[c("output")]) 
    inputs_syn8 <- data.frame(dat_syn8[c("land","labor","material","capstock")])
    outputs_syn8 <- data.frame(dat_syn8[c("output")]) 
    inputs_syn9 <- data.frame(dat_syn9[c("land","labor","material","capstock")])
    outputs_syn9 <- data.frame(dat_syn9[c("output")]) 
    inputs_syn10 <- data.frame(dat_syn10[c("land","labor","material","capstock")])
    outputs_syn10 <- data.frame(dat_syn10[c("output")]) 
    
  # Estimate output-oriented distance functions
  
    # Original data
    odf_orig <- dea(inputs_orig, outputs_orig, RTS = "vrs", ORIENTATION = "out")
    te_dea_orig <- 1/eff(odf_orig, asInData=TRUE) #stores the TE scores
    summary(te_dea_orig)
    odf_orig <- NULL

    # Synthetic data 1
    odf_syn1 <- dea(inputs_syn1, outputs_syn1, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn1 <- 1/eff(odf_syn1, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn1)
    odf_syn1 <- NULL
    
    # Synthetic data 2
    odf_syn2 <- dea(inputs_syn2, outputs_syn2, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn2 <- 1/eff(odf_syn2, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn2)
    odf_syn2 <- NULL
    
    # Synthetic data 3
    odf_syn3 <- dea(inputs_syn3, outputs_syn3, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn3 <- 1/eff(odf_syn3, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn3)
    odf_syn3 <- NULL
    
    # Synthetic data 4
    odf_syn4 <- dea(inputs_syn4, outputs_syn4, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn4 <- 1/eff(odf_syn4, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn4)
    odf_syn4 <- NULL
    
    # Synthetic data 5
    odf_syn5 <- dea(inputs_syn5, outputs_syn5, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn5 <- 1/eff(odf_syn5, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn5)
    odf_syn5 <- NULL
    
    # Synthetic data 6
    odf_syn6 <- dea(inputs_syn6, outputs_syn6, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn6 <- 1/eff(odf_syn6, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn6)
    odf_syn6 <- NULL
    
    # Synthetic data 7
    odf_syn7 <- dea(inputs_syn7, outputs_syn7, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn7 <- 1/eff(odf_syn7, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn7)
    odf_syn7 <- NULL
    
    # Synthetic data 8
    odf_syn8 <- dea(inputs_syn8, outputs_syn8, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn8 <- 1/eff(odf_syn8, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn8)
    odf_syn8 <- NULL
    
    # Synthetic data 9
    odf_syn9 <- dea(inputs_syn9, outputs_syn9, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn9 <- 1/eff(odf_syn9, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn9)
    odf_syn9 <- NULL
    
    # Synthetic data 10
    odf_syn10 <- dea(inputs_syn10, outputs_syn10, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn10 <- 1/eff(odf_syn10, asInData=TRUE) #stores the TE scores
    summary(te_dea_syn10)
    odf_syn10 <- NULL

#------------#
# Table A.8. #
#------------#
    
  # Combine all TE scores from DEA
  descr_dea_out <- cbind(te_dea_orig,te_dea_syn1,te_dea_syn2,te_dea_syn3,te_dea_syn4,te_dea_syn5,
                         te_dea_syn6,te_dea_syn7,te_dea_syn8,te_dea_syn9,te_dea_syn10)

  # Create table with descriptive statistics
  stargazer(as.data.frame(descr_dea_out),
            digits=2,
            summary = TRUE,
            iqr = TRUE,
            title="Output-oriented technical efficiency scores from DEA",
            nobs = FALSE,
            covariate.labels=c("Original data", "Synthetic data 1", 
                               "Synthetic data 2", "Synthetic data 3",
                               "Synthetic data 4","Synthetic data 5",
                               "Synthetic data 6","Synthetic data 7",
                               "Synthetic data 8","Synthetic data 9",
                               "Synthetic data 10"),
            type = "text",
            style = "ajps",
            out = "Tables/FADN/TabA8_descr_dea_out.txt")

  
#--------------#
# Figure A.10. #
#--------------#
  
  # Make data frame with technical efficiency scores from all datasets
  
    # Original and CART data
    nr_obs <- length(te_dea_orig)
    efficiencies_cart <- data.frame(matrix(ncol = 2, nrow = 6*nr_obs))
    colnames(efficiencies_cart) <- c("TE", "data") 
    efficiencies_cart$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                                rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                                rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs))
    efficiencies_cart$TE <- c(te_dea_orig,te_dea_syn1,te_dea_syn2,te_dea_syn3,te_dea_syn4,te_dea_syn5)
    
    # Original and NORMRANK data
    efficiencies_normrank <- data.frame(matrix(ncol = 2, nrow = 6*nr_obs))
    colnames(efficiencies_normrank) <- c("TE", "data") 
    efficiencies_normrank$data <- c(rep("Original data", nr_obs),rep("Synthetic data 6", nr_obs),
                               rep("Synthetic data 7", nr_obs),rep("Synthetic data 8", nr_obs),
                               rep("Synthetic data 9", nr_obs),rep("Synthetic data 10", nr_obs))
    efficiencies_normrank$TE <- c(te_dea_orig,te_dea_syn6,te_dea_syn7,te_dea_syn8,te_dea_syn9,te_dea_syn10)

  # Make sure data are plotted in right order 
        
    efficiencies_cart$data <- factor(efficiencies_cart$data,
                                     levels=c("Synthetic data 5",
                                              "Synthetic data 4",
                                              "Synthetic data 3",
                                              "Synthetic data 2",
                                              "Synthetic data 1",
                                              "Original data"))
    efficiencies_normrank$data <- factor(efficiencies_normrank$data,
                                    levels=c("Synthetic data 10",
                                             "Synthetic data 9",
                                             "Synthetic data 8",
                                             "Synthetic data 7",
                                             "Synthetic data 6",
                                             "Original data"))
    
  # Make plot for CART
  te_dea_out_cart <-     
    efficiencies_cart %>%
    ggplot(aes(x = TE , y = data,  fill = data)) +
    ggridges::geom_density_ridges(scale = 1, 
                                  stat = "binline", 
                                  binwidth = 0.01) +
    theme_light(base_family="Times") + 
    xlab("Technical efficiency") + ylab("") +
    ggtitle("CART method") +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      plot.title = element_text(size=10, hjust = 0.5),
      axis.text.y=element_blank() 
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(0,1.01)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 1", 
                                 "Synthetic data 2", 
                                 "Synthetic data 3", 
                                 "Synthetic data 4", 
                                 "Synthetic data 5"))
  
  # Make plot for NORMRANK
  te_dea_out_normrank <-     
    efficiencies_normrank %>%
    ggplot(aes(x = TE , y = data,  fill = data)) +
    ggridges::geom_density_ridges(scale = 1, 
                                  stat = "binline", 
                                  binwidth = 0.01) +
    theme_light(base_family="Times") + 
    xlab("Technical efficiency") + ylab("") +
    ggtitle("NORMRANK method") +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(colour = "white"),
      legend.title = element_blank(),
      plot.title = element_text(size=10, hjust = 0.5),
      axis.text.y=element_blank() 
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits=c(0,1.01)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 6", 
                                 "Synthetic data 7", 
                                 "Synthetic data 8", 
                                 "Synthetic data 9", 
                                 "Synthetic data 10")) 
  
  # Save plots
  te_dea_out <- ggarrange(te_dea_out_cart, te_dea_out_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/FADN/FigA10_te_dea_out.png", te_dea_out, width = 23, height = 11, units = "cm")

# --------------------------------------- #
#            End of the script            #
# --------------------------------------- #

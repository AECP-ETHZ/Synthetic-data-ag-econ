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
# The program stores all tables and figures based on the Norwegian dairy       #
# production data presented in the article in the subfolders 'Tables/Norway'   #
# and 'Figures/Norway'. The data are automatically downloaded from the         #
# internet. The script was successfully run with R version 4.2.0.              #
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
dir.create(file.path("Tables", "Norway"))
dir.create(file.path("Figures"))
dir.create(file.path("Figures", "Norway"))

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

# Load Norwegian dairy data
temp <- tempfile()
download.file("https://sites.google.com/site/sfbook2014/home/for-stata-v12-v13-v14/sfbook_data_12.zip",temp)
dat <- read_dta(unz(temp, "norway.dta"))
head(dat)

# Define colour palettes
cbPalette <- c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e")
bwPalette <- c("#252525", "#636363", "#969696", "#bdbdbd", "#d9d9d9", "#f7f7f7")
  
# ------------------------------------- #
#### 1. Generation of synthetic data ####
# ------------------------------------- #

# Set seed for replication
my.seed <- 1234

# Define original data set
dat_orig <- as.data.frame(dat[c("y1", "x1", "x2", "x3", "x4", "x5", "x6")])

# ------------------------------------- #
##### 1.a Synthetic data using CART #####
# ------------------------------------- #

# Generate five synthetic data sets (cart)
dat_synt <- syn(dat_orig, method = "cart", m = 5, visit.sequence = 
                  c("y1", "x1", "x2", "x3", "x4",
                    "x5", "x6"),
                seed = my.seed,
                smoothing = list(y1 = "density", 
                                 x1 = "density", 
                                 x2 = "density", 
                                 x3 = "density", 
                                 x4 = "density", 
                                 x5 = "density", 
                                 x6 = "density"))

# Store synthetic data in data frames
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
                  c("y1", "x1", "x2", "x3", "x4",
                    "x5", "x6"),
                seed = my.seed,
                smoothing = list(y1 = "density", 
                                 x1 = "density", 
                                 x2 = "density", 
                                 x3 = "density", 
                                 x4 = "density", 
                                 x5 = "density", 
                                 x6 = "density"))

# Store synthetic data in data frames
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
          title="Descriptive statistics for the Norwegian sample (original)",
          nobs = FALSE,
          omit.summary.stat = c("p25", "p75"),
          covariate.labels=c("Milk", "Land", "Labor",
                             "Feed","Material",
                             "Cattle","Capital"),
          type = "text",
          style = "ajps",
          out = "Tables/Norway/TabA1a.descr_orig.txt")

# One synthetic data: CART
stargazer::stargazer(dat_syn1,
                digits=2,
                title="Descriptive statistics for the Norwegian sample (CART)",
                nobs = FALSE,
                omit.summary.stat = c("p25", "p75"),
                covariate.labels=c("Milk", "Land", "Labor",
                                  "Feed","Material",
                                  "Cattle","Capital"),
                type = "text",
                style = "ajps",
                out = "Tables/Norway/TabA1b.descr_cart.txt")

# One synthetic data: NORMRANK
stargazer::stargazer(dat_syn6,
            digits=2,
            title="Descriptive statistics for the Norwegian sample (NORMRANK)",
            nobs = FALSE,
            omit.summary.stat = c("p25", "p75"),
            covariate.labels=c("Milk", "Land", "Labor",
                               "Feed","Material",
                               "Cattle","Capital"),
            type = "text",
            style = "ajps",
            out = "Tables/Norway/TabA1c.descr_normrank.txt")

# -------------------------------------------- #
#### 3. Estimation of production function   ####
# -------------------------------------------- #

# Define Cobb-Douglas functional form
fn.CD <- log(y1) ~ log(x1) + log(x2) + 
  log(x3) + log(x4) + log(x5) + log(x6)

# Define (de-meaned) translog functional form
fn.TL <- log(y1) ~ log(x1/mean(x1)) + log(x2/mean(x2)) + log(x3/mean(x3)) + 
  log(x4/mean(x4)) + log(x5/mean(x5)) + log(x6/mean(x6)) + 
  I(0.5*log(x1/mean(x1))*log(x1/mean(x1))) + I(log(x1/mean(x1))*log(x2/mean(x2))) + 
  I(log(x1/mean(x1))*log(x3/mean(x3))) + I(log(x1/mean(x1))*log(x4/mean(x4))) + 
  I(log(x1/mean(x1))*log(x5/mean(x5))) + I(log(x1/mean(x1))*log(x6/mean(x6))) + 
  I(0.5*log(x2/mean(x2))*log(x2/mean(x2))) + I(log(x2/mean(x2))*log(x3/mean(x3))) + 
  I(log(x2/mean(x2))*log(x4/mean(x4))) + I(log(x2/mean(x2))*log(x5/mean(x5))) + 
  I(log(x2/mean(x2))*log(x6/mean(x6))) + I(0.5*log(x3/mean(x3))*log(x3/mean(x3))) + 
  I(log(x3/mean(x3))*log(x4/mean(x4))) + I(log(x3/mean(x3))*log(x5/mean(x5))) + 
  I(log(x3/mean(x3))*log(x6/mean(x6))) + I(0.5*log(x4/mean(x4))*log(x4/mean(x4))) + 
  I(log(x4/mean(x4))*log(x5/mean(x5))) + I(log(x4/mean(x4))*log(x6/mean(x6))) +
  I(0.5*log(x5/mean(x5))*log(x5/mean(x5))) + I(log(x5/mean(x5))*log(x6/mean(x6))) + 
  I(0.5*log(x6/mean(x6))*log(x6/mean(x6)))

# ---------------------------------------------#
##### 3a. Cobb Douglas Production Function #####
# ---------------------------------------------#

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
  
# --------------------------------------------------------------------- #
# Figure 1. Parameter estimates of the Cobb-Douglas production function #
#  with 95% confidence intervals                                        #
# --------------------------------------------------------------------- #
  
  # Combine all models
    
    models_cart <- bind_rows(prod.orig, prod.syn1, prod.syn2, prod.syn3,
                             prod.syn4, prod.syn5)
    models_normrank <- bind_rows(prod.orig, prod.syn6, prod.syn7, prod.syn8,
                                 prod.syn9, prod.syn10)
  
  # Make plot for CART
  elast_CD_prod_cart <- dwplot(models_cart, 
                          show_intercept = FALSE,
                          ci = .95,
                          vars_order = c("log(x1)", "log(x2)", "log(x3)", 
                                         "log(x4)", "log(x5)", "log(x6)"),
                          model_order = c("Original data", 
                                         "Synthetic data 1", 
                                         "Synthetic data 2", 
                                         "Synthetic data 3", 
                                         "Synthetic data 4", 
                                         "Synthetic data 5")) %>% 
                          relabel_predictors(
                                  c("log(x1)" = "Land",
                                    "log(x2)" = "Labor",
                                    "log(x3)" = "Feed",
                                    "log(x4)" = "Material",
                                    "log(x5)" = "Cattle",
                                    "log(x6)" = "Capital")) +
                          theme_light(base_family="Times") + 
                          xlab("Estimated elasticities with 95% CIs") + 
                          ylab("") +
                          ggtitle("CART method") +
                          theme(
                            plot.title = element_text(size=10, hjust = 0.5),
                            legend.position = "bottom",
                            legend.background = element_rect(colour = "white"),
                            legend.title = element_blank()) +
                          scale_x_continuous(breaks = seq(0, 0.3, by = 0.1), limits=c(-0.0001,0.35)) +
                          scale_colour_manual(values=cbPalette)
        
  # Make plot for NORMRANK
  elast_CD_prod_normrank <- dwplot(models_normrank, 
                        show_intercept = FALSE,
                        ci = .95,
                        vars_order = c("log(x1)", "log(x2)", "log(x3)", 
                                       "log(x4)", "log(x5)", "log(x6)"),
                        model_order = c("Original data", 
                                        "Synthetic data 6", 
                                        "Synthetic data 7", 
                                        "Synthetic data 8", 
                                        "Synthetic data 9", 
                                        "Synthetic data 10")) %>% 
                        relabel_predictors(
                            c("log(x1)" = "Land",
                              "log(x2)" = "Labor",
                              "log(x3)" = "Feed",
                              "log(x4)" = "Material",
                              "log(x5)" = "Cattle",
                              "log(x6)" = "Capital")
                        ) +
                        theme_light(base_family="Times") + 
                        xlab("Estimated elasticities with 95% CIs") + ylab("") +
                        ggtitle("NORMRANK method") +
                        theme(
                          plot.title = element_text(size=10, hjust = 0.5),
                          legend.position = "bottom",
                          legend.background = element_rect(colour = "white"),
                          legend.title = element_blank()) +
                        scale_x_continuous(breaks = seq(0, 0.3, by = 0.1), limits=c(-0.0001,0.35)) +
                        scale_colour_manual(values=cbPalette)

  # Save plots
  elast_CD_prod <- ggarrange(elast_CD_prod_cart, elast_CD_prod_normrank, 
                             ncol = 2, nrow = 1)
  ggsave("Figures/Norway/Fig1.elast_CD_prod.png", elast_CD_prod, 
         width = 24, height = 11, units = "cm")

  
  # Same figures in black and white for journal
  
    # Make plot for CART
    elast_CD_prod_cart <- dwplot(models_cart, 
                                 dot_args = list(aes(shape = model)),
                                 dodge_size = .6, 
                                 show_intercept = FALSE,
                                 ci = .95,
                                 vars_order = c("log(x1)", "log(x2)", "log(x3)", 
                                                "log(x4)", "log(x5)", "log(x6)"),
                                 model_order = c("Original data", 
                                                 "Synthetic data 1", 
                                                 "Synthetic data 2", 
                                                 "Synthetic data 3", 
                                                 "Synthetic data 4", 
                                                 "Synthetic data 5")) %>% 
                                relabel_predictors( c("log(x1)" = "Land",
                                                      "log(x2)" = "Labor",
                                                      "log(x3)" = "Feed",
                                                      "log(x4)" = "Material",
                                                      "log(x5)" = "Cattle",
                                                      "log(x6)" = "Capital")) +
                                theme_bw() + 
                                xlab("Estimated elasticities with 95% CIs") + 
                                ylab("") +
                                ggtitle("CART method") +
                                theme(
                                    plot.title = element_text(size=10, hjust = 0.5),
                                    legend.position = "bottom",
                                    legend.background = element_rect(colour = "white"),
                                    legend.title = element_blank(),
                                    ## remove the vertical grid lines
                                    panel.grid.minor.y = element_blank() ,
                                    panel.grid.major.y = element_blank() ,
                                    # set the horizontal lines
                                    panel.grid.minor.x = element_line( size=.05, color="lightgrey" ),
                                    panel.grid.major.x = element_line( size=.05, color="lightgrey" )) +
                                scale_x_continuous(breaks = seq(0, 0.3, by = 0.1), limits=c(-0.0001,0.35)) +
                                scale_colour_grey(name = "Model", 
                                                  start = .3, end = .3, # Same colour for all models
                                                  breaks = c("Original data", 
                                                             "Synthetic data 1", 
                                                             "Synthetic data 2", 
                                                             "Synthetic data 3", 
                                                             "Synthetic data 4", 
                                                             "Synthetic data 5"),
                                                  labels = c("Original data", 
                                                             "Synthetic data 1", 
                                                             "Synthetic data 2", 
                                                             "Synthetic data 3", 
                                                             "Synthetic data 4", 
                                                             "Synthetic data 5")) +
                                scale_shape_manual(name = "Model", # Different shapes for models
                                                   values = c(16, 0, 1, 2, 5, 6),
                                                   breaks = c("Original data", 
                                                              "Synthetic data 1", 
                                                              "Synthetic data 2", 
                                                              "Synthetic data 3", 
                                                              "Synthetic data 4", 
                                                              "Synthetic data 5"),
                                                   labels = c("Original data", 
                                                              "Synthetic data 1", 
                                                              "Synthetic data 2", 
                                                              "Synthetic data 3", 
                                                              "Synthetic data 4", 
                                                              "Synthetic data 5")) +
                                guides(
                                  shape = guide_legend("Model"), 
                                  colour = guide_legend("Model")
                                ) # Combine the legends for shape and colour
    
    # Make plot for NORMRANK
    elast_CD_prod_normrank <- dwplot(models_normrank, 
                                     dot_args = list(aes(shape = model)),
                                     dodge_size = .6, 
                                     show_intercept = FALSE,
                                     ci = .95,
                                     vars_order = c("log(x1)", "log(x2)", "log(x3)", 
                                                    "log(x4)", "log(x5)", "log(x6)"),
                                     model_order = c("Original data", 
                                                     "Synthetic data 6", 
                                                     "Synthetic data 7", 
                                                     "Synthetic data 8", 
                                                     "Synthetic data 9", 
                                                     "Synthetic data 10")) %>% 
                                    relabel_predictors( c("log(x1)" = "Land",
                                                          "log(x2)" = "Labor",
                                                          "log(x3)" = "Feed",
                                                          "log(x4)" = "Material",
                                                          "log(x5)" = "Cattle",
                                                          "log(x6)" = "Capital")) +
                                    theme_bw() + 
                                    xlab("Estimated elasticities with 95% CIs") + 
                                    ylab("") +
                                    ggtitle("NORMRANK method") +
                                    theme(
                                        plot.title = element_text(size=10, hjust = 0.5),
                                        legend.position = "bottom",
                                        legend.background = element_rect(colour = "white"),
                                        legend.title = element_blank(),
                                        # remove the vertical grid lines
                                        panel.grid.minor.y = element_blank() ,
                                        panel.grid.major.y = element_blank() ,
                                        # set the horizontal lines
                                        panel.grid.minor.x = element_line( size=.05, color="lightgrey" ),
                                        panel.grid.major.x = element_line( size=.05, color="lightgrey" )) +
                                    scale_x_continuous(breaks = seq(0, 0.3, by = 0.1), limits=c(-0.0001,0.35)) +
                                    scale_colour_grey(name = "Model", 
                                                      start = .3, end = .3, # Same colour for all models
                                                      breaks = c("Original data", 
                                                                 "Synthetic data 6", 
                                                                 "Synthetic data 7", 
                                                                 "Synthetic data 8", 
                                                                 "Synthetic data 9", 
                                                                 "Synthetic data 10"),
                                                      labels = c("Original data", 
                                                                 "Synthetic data 6", 
                                                                 "Synthetic data 7", 
                                                                 "Synthetic data 8", 
                                                                 "Synthetic data 9", 
                                                                 "Synthetic data 10")) +
                                    scale_shape_manual(name = "Model", # Different shapes for models
                                                       values = c(16, 0, 1, 2, 5, 6),
                                                       breaks = c("Original data", 
                                                                  "Synthetic data 6", 
                                                                  "Synthetic data 7", 
                                                                  "Synthetic data 8", 
                                                                  "Synthetic data 9", 
                                                                  "Synthetic data 10"),
                                                       labels = c("Original data", 
                                                                  "Synthetic data 6", 
                                                                  "Synthetic data 7", 
                                                                  "Synthetic data 8", 
                                                                  "Synthetic data 9", 
                                                                  "Synthetic data 10")) +
                                    guides(
                                      shape = guide_legend("Model"), 
                                      colour = guide_legend("Model")
                                    ) # Combine the legends for shape and colour
    
    # Save plots
    elast_CD_prod <- ggarrange(elast_CD_prod_cart, elast_CD_prod_normrank, 
                               ncol = 2, nrow = 1)
    ggsave("Figures/Norway/Fig1.elast_CD_prod_bw.png", elast_CD_prod, 
           width = 24, height = 11, units = "cm")
  
    
# -----------------------------------------#
##### 3b. Translog Production Function #####
# -----------------------------------------#
  
# Create functions to estimate production function and calculate elasticities

  prod_tl <- function(data) {
  prod <- lm(fn.TL, data = data)
  coef <- prod$coefficients
  el_land <- coef["log(x1/mean(x1))"] + 
    coef["I(0.5 * log(x1/mean(x1)) * log(x1/mean(x1)))"]*log(data$x1/mean(data$x1)) + 
    coef["I(log(x1/mean(x1)) * log(x2/mean(x2)))"]*log(data$x2/mean(data$x2)) +
    coef["I(log(x1/mean(x1)) * log(x3/mean(x3)))"]*log(data$x3/mean(data$x3)) +
    coef["I(log(x1/mean(x1)) * log(x4/mean(x4)))"]*log(data$x4/mean(data$x4)) +
    coef["I(log(x1/mean(x1)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
    coef["I(log(x1/mean(x1)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
  el_labor <- coef["log(x2/mean(x2))"] + 
    coef["I(log(x1/mean(x1)) * log(x2/mean(x2)))"]*log(data$x1/mean(data$x1)) + 
    coef["I(0.5 * log(x2/mean(x2)) * log(x2/mean(x2)))"]*log(data$x2/mean(data$x2)) +
    coef["I(log(x2/mean(x2)) * log(x3/mean(x3)))"]*log(data$x3/mean(data$x3)) +
    coef["I(log(x2/mean(x2)) * log(x4/mean(x4)))"]*log(data$x4/mean(data$x4)) +
    coef["I(log(x2/mean(x2)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
    coef["I(log(x2/mean(x2)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
  el_feed <- coef["log(x3/mean(x3))"] + 
    coef["I(log(x1/mean(x1)) * log(x3/mean(x3)))"]*log(data$x1/mean(data$x1)) + 
    coef["I(log(x2/mean(x2)) * log(x3/mean(x3)))"]*log(data$x2/mean(data$x2)) +
    coef["I(0.5 * log(x3/mean(x3)) * log(x3/mean(x3)))"]*log(data$x3/mean(data$x3)) +
    coef["I(log(x3/mean(x3)) * log(x4/mean(x4)))"]*log(data$x4/mean(data$x4)) +
    coef["I(log(x3/mean(x3)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
    coef["I(log(x3/mean(x3)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
  el_material <- coef["log(x4/mean(x4))"] + 
    coef["I(log(x1/mean(x1)) * log(x4/mean(x4)))"]*log(data$x1/mean(data$x1)) + 
    coef["I(log(x2/mean(x2)) * log(x4/mean(x4)))"]*log(data$x2/mean(data$x2)) +
    coef["I(log(x3/mean(x3)) * log(x4/mean(x4)))"]*log(data$x3/mean(data$x3)) +
    coef["I(0.5 * log(x4/mean(x4)) * log(x4/mean(x4)))"]*log(data$x4/mean(data$x4)) +
    coef["I(log(x4/mean(x4)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
    coef["I(log(x4/mean(x4)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
  el_cattle <- coef["log(x5/mean(x5))"] + 
    coef["I(log(x1/mean(x1)) * log(x5/mean(x5)))"]*log(data$x1/mean(data$x1)) + 
    coef["I(log(x2/mean(x2)) * log(x5/mean(x5)))"]*log(data$x2/mean(data$x2)) +
    coef["I(log(x3/mean(x3)) * log(x5/mean(x5)))"]*log(data$x3/mean(data$x3)) +
    coef["I(log(x4/mean(x4)) * log(x5/mean(x5)))"]*log(data$x4/mean(data$x4)) +
    coef["I(0.5 * log(x5/mean(x5)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
    coef["I(log(x5/mean(x5)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
  el_capital <- coef["log(x6/mean(x6))"] + 
    coef["I(log(x1/mean(x1)) * log(x6/mean(x6)))"]*log(data$x1/mean(data$x1)) + 
    coef["I(log(x2/mean(x2)) * log(x6/mean(x6)))"]*log(data$x2/mean(data$x2)) +
    coef["I(log(x3/mean(x3)) * log(x6/mean(x6)))"]*log(data$x3/mean(data$x3)) +
    coef["I(log(x4/mean(x4)) * log(x6/mean(x6)))"]*log(data$x4/mean(data$x4)) +
    coef["I(log(x5/mean(x5)) * log(x6/mean(x6)))"]*log(data$x5/mean(data$x5)) +
    coef["I(0.5 * log(x6/mean(x6)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
  result <- list()
  result$elast <- cbind(el_land, el_labor, el_feed, el_material, el_cattle, el_capital)
  result$model <- tidy(prod)
  return(result)
  }
  
# Estimate TL prod and calculate elasticities

  # Original data
  prod <- prod_tl(dat_orig)
  prod.orig <- prod$model[c(1:7),] %>% 
    mutate(model = "Original data") # adds model label
  dat_orig <- cbind(dat_orig,prod$elast)

  # Synthetic data 1
  prod <- prod_tl(dat_syn1)
  prod.syn1 <- prod$model[c(1:7),] %>% 
    mutate(model = "Synthetic data 1") 
  dat_syn1 <- cbind(dat_syn1,prod$elast)
  
  # Synthetic data 2
  prod <- prod_tl(dat_syn2)
  prod.syn2 <- prod$model[c(1:7),] %>% 
    mutate(model = "Synthetic data 2") 
  dat_syn2 <- cbind(dat_syn2,prod$elast)
  
  # Synthetic data 3
  prod <- prod_tl(dat_syn3)
  prod.syn3 <- prod$model[c(1:7),] %>% 
    mutate(model = "Synthetic data 3") 
  dat_syn3 <- cbind(dat_syn3,prod$elast)
    
  # Synthetic data 4
  prod <- prod_tl(dat_syn4)
  prod.syn4 <- prod$model[c(1:7),] %>% 
    mutate(model = "Synthetic data 4") 
  dat_syn4 <- cbind(dat_syn4,prod$elast)
    
  # Synthetic data 5
  prod <- prod_tl(dat_syn5)
  prod.syn5 <- prod$model[c(1:7),] %>%
    mutate(model = "Synthetic data 5") 
  dat_syn5 <- cbind(dat_syn5,prod$elast)
    
  # Synthetic data 6
  prod <- prod_tl(dat_syn6)
  prod.syn6 <- prod$model[c(1:7),] %>% 
    mutate(model = "Synthetic data 6") 
  dat_syn6 <- cbind(dat_syn6,prod$elast)
    
  # Synthetic data 7
  prod <- prod_tl(dat_syn7)
  prod.syn7 <- prod$model[c(1:7),] %>% 
    mutate(model = "Synthetic data 7") 
  dat_syn7 <- cbind(dat_syn7,prod$elast)
    
  # Synthetic data 8
  prod <- prod_tl(dat_syn8)
  prod.syn8 <- prod$model[c(1:7),] %>% 
    mutate(model = "Synthetic data 8") 
  dat_syn8 <- cbind(dat_syn8,prod$elast)
    
  # Synthetic data 9
  prod <- prod_tl(dat_syn9)
  prod.syn9 <- prod$model[c(1:7),] %>% 
    mutate(model = "Synthetic data 9") 
  dat_syn9 <- cbind(dat_syn9,prod$elast)
    
  # Synthetic data 10
  prod <- prod_tl(dat_syn10)
  prod.syn10 <- prod$model[c(1:7),] %>% 
    mutate(model = "Synthetic data 10") 
  dat_syn10 <- cbind(dat_syn10,prod$elast)
    
# ---------------------------------------------------------------------- #
# Figure A.2. First-order parameter estimates of the Translog production #
# function with 95% confidence intervals                                 #
# ---------------------------------------------------------------------- #
  
  # Combine all models
    
    models_cart <- bind_rows(prod.orig, prod.syn1, prod.syn2, prod.syn3,
                             prod.syn4, prod.syn5)
    models_normrank <- bind_rows(prod.orig, prod.syn6, prod.syn7, prod.syn8,
                                 prod.syn9, prod.syn10)
    
  # Make plot for CART
  elast_TL_prod_cart <- dwplot(models_cart, 
                          show_intercept = FALSE,
                          ci = .95,
                          vars_order = c("log(x1/mean(x1))", "log(x2/mean(x2))", 
                                        "log(x3/mean(x3))", "log(x4/mean(x4))", 
                                        "log(x5/mean(x5))", "log(x6/mean(x6))"),
                          model_order = c("Original data", 
                                         "Synthetic data 1", 
                                         "Synthetic data 2", 
                                         "Synthetic data 3", 
                                         "Synthetic data 4", 
                                         "Synthetic data 5")) %>% 
                          relabel_predictors(
                            c("log(x1/mean(x1))" = "Land",
                              "log(x2/mean(x2))" = "Labor",
                              "log(x3/mean(x3))" = "Feed",
                              "log(x4/mean(x4))" = "Material",
                              "log(x5/mean(x5))" = "Cattle",
                              "log(x6/mean(x6))" = "Capital")) +
                          theme_light(base_family="Times") + 
                          xlab("Estimated elasticities with 95% CIs") + 
                          ylab("") +
                          ggtitle("CART method") +
                          theme(
                            plot.title = element_text(size=10, hjust = 0.5),
                            legend.position = "bottom",
                            legend.background = element_rect(colour = "white"),
                            legend.title = element_blank() ) +
                          scale_x_continuous(breaks = seq(0, 0.4, by = 0.1), limits=c(-0.025,0.4)) +
                          scale_colour_manual(values=cbPalette)
    
  # Make plot for NORMRANK
  elast_TL_prod_normrank <- dwplot(models_normrank, 
                          show_intercept = FALSE,
                          ci = .95,
                          vars_order = c("log(x1/mean(x1))", "log(x2/mean(x2))", 
                                        "log(x3/mean(x3))", "log(x4/mean(x4))", 
                                        "log(x5/mean(x5))", "log(x6/mean(x6))"),
                          model_order = c("Original data", 
                                         "Synthetic data 6", 
                                         "Synthetic data 7", 
                                         "Synthetic data 8", 
                                         "Synthetic data 9", 
                                         "Synthetic data 10")) %>% 
                          relabel_predictors(
                            c("log(x1/mean(x1))" = "Land",
                              "log(x2/mean(x2))" = "Labor",
                              "log(x3/mean(x3))" = "Feed",
                              "log(x4/mean(x4))" = "Material",
                              "log(x5/mean(x5))" = "Cattle",
                              "log(x6/mean(x6))" = "Capital")) +
                          theme_light(base_family="Times") + 
                          xlab("Estimated elasticities with 95% CIs") + 
                          ylab("") +
                          ggtitle("NORMRANK method") +
                          theme(
                            plot.title = element_text(size=10, hjust = 0.5),
                            legend.position = "bottom",
                            legend.background = element_rect(colour = "white"),
                            legend.title = element_blank()) + 
                          scale_x_continuous(breaks = seq(0, 0.4, by = 0.1), limits=c(-0.025,0.4)) +
                          scale_colour_manual(values=cbPalette) 
  
  # Save plots
  elast_TL_prod <- ggarrange(elast_TL_prod_cart, elast_TL_prod_normrank, 
                             ncol = 2, nrow = 1)
  ggsave("Figures/Norway/FigA2.elast_TL_prod.png", elast_TL_prod, 
         width = 23, height = 11, units = "cm")
  
# --------------------------------------------------------------------- #
# Figure A.3. Distribution of farm-level elasticities obtained from the #
#             Translog production function                              #
# --------------------------------------------------------------------- #
  
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
      el_land$el <- c(dat_orig$el_land,dat_syn1$el_land,dat_syn2$el_land,
                      dat_syn3$el_land,dat_syn4$el_land,dat_syn5$el_land,
                      dat_syn6$el_land,dat_syn7$el_land,dat_syn8$el_land,
                      dat_syn9$el_land,dat_syn10$el_land)
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
      el_labor$el <- c(dat_orig$el_labor,dat_syn1$el_labor,dat_syn2$el_labor,
                       dat_syn3$el_labor,dat_syn4$el_labor,dat_syn5$el_labor,
                       dat_syn6$el_labor,dat_syn7$el_labor,dat_syn8$el_labor,
                       dat_syn9$el_labor,dat_syn10$el_labor)
      el_labor$inp <- "Labor"
      
      # Feed
      el_feed <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
      colnames(el_feed) <- c("el", "data") 
      el_feed$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                        rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                        rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                        rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                        rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                        rep("Synthetic data 10", nr_obs))
      el_feed$el <- c(dat_orig$el_feed,dat_syn1$el_feed,dat_syn2$el_feed,
                      dat_syn3$el_feed,dat_syn4$el_feed,dat_syn5$el_feed,
                      dat_syn6$el_feed,dat_syn7$el_feed,dat_syn8$el_feed,
                      dat_syn9$el_feed,dat_syn10$el_feed)
      el_feed$inp <- "Feed"
      
      # Material
      el_material <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
      colnames(el_material) <- c("el", "data") 
      el_material$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                            rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                            rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                            rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                            rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                            rep("Synthetic data 10", nr_obs))
      el_material$el <- c(dat_orig$el_material,dat_syn1$el_material,dat_syn2$el_material,
                          dat_syn3$el_material,dat_syn4$el_material,dat_syn5$el_material,
                          dat_syn6$el_material,dat_syn7$el_material,dat_syn8$el_material,
                          dat_syn9$el_material,dat_syn10$el_material)
      el_material$inp <- "Material"
      
      # Cattle
      el_cattle <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
      colnames(el_cattle) <- c("el", "data") 
      el_cattle$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                          rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                          rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                          rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                          rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                          rep("Synthetic data 10", nr_obs))
      el_cattle$el <- c(dat_orig$el_cattle,dat_syn1$el_cattle,dat_syn2$el_cattle,
                        dat_syn3$el_cattle,dat_syn4$el_cattle,dat_syn5$el_cattle,
                        dat_syn6$el_cattle,dat_syn7$el_cattle,dat_syn8$el_cattle,
                        dat_syn9$el_cattle,dat_syn10$el_cattle)
      el_cattle$inp <- "Cattle"
      
      # Capital
      el_capital <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
      colnames(el_capital) <- c("el", "data") 
      el_capital$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                           rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                           rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                           rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                           rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                           rep("Synthetic data 10", nr_obs))
      el_capital$el <- c(dat_orig$el_capital,dat_syn1$el_capital,dat_syn2$el_capital,
                         dat_syn3$el_capital,dat_syn4$el_capital,dat_syn5$el_capital,
                         dat_syn6$el_capital,dat_syn7$el_capital,dat_syn8$el_capital,
                         dat_syn9$el_capital,dat_syn10$el_capital)
      el_capital$inp <- "Capital"
      
      # All
      el_prod <- rbind(el_land,el_labor,el_feed,el_material,el_cattle,el_capital)
      el_prod_cart <- el_prod %>% # for CART plot
        filter(data == "Original data" | data == "Synthetic data 1" |
                 data == "Synthetic data 2" | data == "Synthetic data 3" |
                 data == "Synthetic data 4" | data == "Synthetic data 5")
      el_prod_normrank <- el_prod %>% # for NORMRANK plot
        filter(data == "Original data" | data == "Synthetic data 6" |
                 data == "Synthetic data 7" | data == "Synthetic data 8" |
                 data == "Synthetic data 9" | data == "Synthetic data 10")
      
  # Make sure data are plotted in right order 
      
    el_prod_cart$inp <- factor(el_prod_cart$inp,
                            levels=c("Land","Labor","Feed","Material",
                                     "Cattle","Capital"))
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
                                 levels=c("Land","Labor","Feed","Material",
                                          "Cattle","Capital"))
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
      scale_y_continuous(breaks = seq(-0.3, 0.6, by = 0.3), limits=c(-0.5,0.6)) +
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
      scale_y_continuous(breaks = seq(-0.3, 0.6, by = 0.3), limits=c(-0.5,0.6)) +
      scale_fill_manual(values=cbPalette) 
  

    # Save plots
    elast_TL_fl_prod <- ggarrange(elast_TL_fl_prod_cart, 
                                  elast_TL_fl_prod_normrank, ncol = 2, nrow = 1)
    ggsave("Figures/Norway/FigA3.elast_TL_fl_prod.png", elast_TL_fl_prod, 
           width = 23, height = 11, units = "cm")
    
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
      # prepare for coefficient plot
      front <- summary(front)
      front <- as_tibble(front$mleParam)
      front$term <- c("Intercept", "Land", "Labor", "Feed", "Material", 
                      "Cattle", "Capital", "Sigma2", "Gamma")
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
  
# ----------------------------------------------------------------------- #
# Figure A.4. Parameter estimates of the Cobb-Douglas production frontier #
#  with 95% confidence intervals                                          #
# ----------------------------------------------------------------------- #
  
  # Combine all models

    models_cart <- bind_rows(front.orig, front.syn1, front.syn2, front.syn3,
                             front.syn4, front.syn5)
    models_normrank <- bind_rows(front.orig, front.syn6, front.syn7, front.syn8,
                                 front.syn9, front.syn10)
    
  # Make plot for CART
  elast_CD_front_cart <- dwplot(models_cart, 
                          show_intercept = FALSE,
                          ci = .95,
                          vars_order = c("Land", "Labor", "Feed", 
                                         "Material", "Cattle", "Capital"),
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
                            legend.title = element_blank()
                          ) +
                          scale_x_continuous(breaks = seq(0, 0.3, by = 0.1), limits=c(-0.015,0.32)) +
                          scale_colour_manual(values=cbPalette)
  
  # Make plot for NORMRANK
  elast_CD_front_normrank <- dwplot(models_normrank, 
                                show_intercept = FALSE,
                                ci = .95,
                                vars_order = c("Land", "Labor", "Feed", 
                                               "Material", "Cattle", "Capital"),
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
                        scale_x_continuous(breaks = seq(0, 0.3, by = 0.1), limits=c(-0.015,0.32)) +
                        scale_colour_manual(values=cbPalette)
  
  # Save plots
  elast_CD_front <- ggarrange(elast_CD_front_cart, elast_CD_front_normrank, 
                              ncol = 2, nrow = 1)
  ggsave("Figures/Norway/FigA4.elast_CD_front.png", elast_CD_front, 
         width = 23, height = 11, units = "cm")
  
# -------------------------------------------------------------------- #
# Table A.2. Technical efficiency scores, Cobb-Douglas functional form #
# -------------------------------------------------------------------- #
  
  # Combine all TE scores
  descr_te_cd <- cbind(te_orig,te_syn1,te_syn2,te_syn3,te_syn4,te_syn5,
                       te_syn6,te_syn7,te_syn8,te_syn9,te_syn10)
  
  # Create table with descriptive statistics
  stargazer(as.data.frame(descr_te_cd),
            digits=2,
            summary = TRUE,
            iqr = TRUE,
            title="Summary statistics: Technical efficiency scores, Cobb-Douglas",
            nobs = FALSE,
            covariate.labels=c("Original data", "Synthetic data 1", 
                               "Synthetic data 2", "Synthetic data 3",
                               "Synthetic data 4","Synthetic data 5",
                               "Synthetic data 6","Synthetic data 7",
                               "Synthetic data 8","Synthetic data 9",
                               "Synthetic data 10"),
            type = "text",
            style = "ajps",
            out = "Tables/Norway/TabA2.descr_te_cd.txt")

# -------------------------------------------------------------- #
# Figure 2. Distribution of technical efficiency scores based on # 
# the Cobb-Douglas production frontier.                          #
# -------------------------------------------------------------- #
  
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
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.15,1)) + 
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
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.15,1)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 6", 
                                 "Synthetic data 7", 
                                 "Synthetic data 8", 
                                 "Synthetic data 9", 
                                 "Synthetic data 10")) 

  # Save plots
  te_cd <- ggarrange(te_cd_cart, te_cd_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/Norway/Fig2.te_cd.png", te_cd, 
         width = 23, height = 11, units = "cm")
  
  # Same figures in black and white for journal
  
    # Make plot for CART
    te_cd_cart <-     
      efficiencies_cart %>%
      ggplot(aes(x = TE , y = data,  fill = data)) +
      ggridges::geom_density_ridges(scale = 1, 
                                    stat = "binline", 
                                    binwidth = 0.01) +
      theme_bw(base_family="Times") + 
      xlab("Technical efficiency") + ylab("") +
      ggtitle("CART method") +
      theme(
        legend.position = "bottom",
        legend.background = element_rect(colour = "white"),
        legend.title = element_blank(),
        plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y=element_blank() 
      ) +
      scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.15,1)) + 
      scale_fill_manual(values=bwPalette,
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
      theme_bw(base_family="Times") + 
      xlab("Technical efficiency") + ylab("") +
      ggtitle("NORMRANK method") +
      theme(
        legend.position = "bottom",
        legend.background = element_rect(colour = "white"),
        legend.title = element_blank(),
        plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y=element_blank() 
      ) +
      scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.15,1)) + 
      scale_fill_manual(values=bwPalette,
                        limits = c("Original data", 
                                   "Synthetic data 6", 
                                   "Synthetic data 7", 
                                   "Synthetic data 8", 
                                   "Synthetic data 9", 
                                   "Synthetic data 10")) 
    
    # Save plots
    te_cd <- ggarrange(te_cd_cart, te_cd_normrank, ncol = 2, nrow = 1)
    ggsave("Figures/Norway/Fig2.te_cd_bw.png", te_cd, 
           width = 23, height = 11, units = "cm")
  
# ----------------------------------------- #
##### 4.b Translog Production Frontier  #####
# ----------------------------------------- # 
  
# Create functions to estimate production frontiers and calculate elasticities
  
  front_tl <- function(data) {
    front <- frontier::sfa(fn.TL, data = data)
    result <- list()
    # Store TE scores
    result$te <- frontier::efficiencies(front, asInData=TRUE) 
    # Calculate elasticities
    coef <- front$mleParam
    el_land_front <- coef["log(x1/mean(x1))"] + 
      coef["I(0.5 * log(x1/mean(x1)) * log(x1/mean(x1)))"]*log(data$x1/mean(data$x1)) + 
      coef["I(log(x1/mean(x1)) * log(x2/mean(x2)))"]*log(data$x2/mean(data$x2)) +
      coef["I(log(x1/mean(x1)) * log(x3/mean(x3)))"]*log(data$x3/mean(data$x3)) +
      coef["I(log(x1/mean(x1)) * log(x4/mean(x4)))"]*log(data$x4/mean(data$x4)) +
      coef["I(log(x1/mean(x1)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
      coef["I(log(x1/mean(x1)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
    el_labor_front <- coef["log(x2/mean(x2))"] + 
      coef["I(log(x1/mean(x1)) * log(x2/mean(x2)))"]*log(data$x1/mean(data$x1)) + 
      coef["I(0.5 * log(x2/mean(x2)) * log(x2/mean(x2)))"]*log(data$x2/mean(data$x2)) +
      coef["I(log(x2/mean(x2)) * log(x3/mean(x3)))"]*log(data$x3/mean(data$x3)) +
      coef["I(log(x2/mean(x2)) * log(x4/mean(x4)))"]*log(data$x4/mean(data$x4)) +
      coef["I(log(x2/mean(x2)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
      coef["I(log(x2/mean(x2)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
    el_feed_front <- coef["log(x3/mean(x3))"] + 
      coef["I(log(x1/mean(x1)) * log(x3/mean(x3)))"]*log(data$x1/mean(data$x1)) + 
      coef["I(log(x2/mean(x2)) * log(x3/mean(x3)))"]*log(data$x2/mean(data$x2)) +
      coef["I(0.5 * log(x3/mean(x3)) * log(x3/mean(x3)))"]*log(data$x3/mean(data$x3)) +
      coef["I(log(x3/mean(x3)) * log(x4/mean(x4)))"]*log(data$x4/mean(data$x4)) +
      coef["I(log(x3/mean(x3)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
      coef["I(log(x3/mean(x3)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
    el_material_front <- coef["log(x4/mean(x4))"] + 
      coef["I(log(x1/mean(x1)) * log(x4/mean(x4)))"]*log(data$x1/mean(data$x1)) + 
      coef["I(log(x2/mean(x2)) * log(x4/mean(x4)))"]*log(data$x2/mean(data$x2)) +
      coef["I(log(x3/mean(x3)) * log(x4/mean(x4)))"]*log(data$x3/mean(data$x3)) +
      coef["I(0.5 * log(x4/mean(x4)) * log(x4/mean(x4)))"]*log(data$x4/mean(data$x4)) +
      coef["I(log(x4/mean(x4)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
      coef["I(log(x4/mean(x4)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
    el_cattle_front <- coef["log(x5/mean(x5))"] + 
      coef["I(log(x1/mean(x1)) * log(x5/mean(x5)))"]*log(data$x1/mean(data$x1)) + 
      coef["I(log(x2/mean(x2)) * log(x5/mean(x5)))"]*log(data$x2/mean(data$x2)) +
      coef["I(log(x3/mean(x3)) * log(x5/mean(x5)))"]*log(data$x3/mean(data$x3)) +
      coef["I(log(x4/mean(x4)) * log(x5/mean(x5)))"]*log(data$x4/mean(data$x4)) +
      coef["I(0.5 * log(x5/mean(x5)) * log(x5/mean(x5)))"]*log(data$x5/mean(data$x5)) +
      coef["I(log(x5/mean(x5)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
    el_capital_front <- coef["log(x6/mean(x6))"] + 
      coef["I(log(x1/mean(x1)) * log(x6/mean(x6)))"]*log(data$x1/mean(data$x1)) + 
      coef["I(log(x2/mean(x2)) * log(x6/mean(x6)))"]*log(data$x2/mean(data$x2)) +
      coef["I(log(x3/mean(x3)) * log(x6/mean(x6)))"]*log(data$x3/mean(data$x3)) +
      coef["I(log(x4/mean(x4)) * log(x6/mean(x6)))"]*log(data$x4/mean(data$x4)) +
      coef["I(log(x5/mean(x5)) * log(x6/mean(x6)))"]*log(data$x5/mean(data$x5)) +
      coef["I(0.5 * log(x6/mean(x6)) * log(x6/mean(x6)))"]*log(data$x6/mean(data$x6))
    # prepare for coefficient plot
    front <- summary(front)
    front <- as_tibble(front$mleParam)
    front <- front[c(1:7),]
    front$term <- c("Intercept", "Land", "Labor", "Feed", 
                    "Material", "Cattle", "Capital")
    front <- front %>% 
      select(term, everything()) #re-order variables
    names(front) <- names(prod.orig)[c(1:5)]
    result$estimates <- front
    result$elast <- cbind(el_land_front, el_labor_front, el_feed_front, 
                          el_material_front, el_cattle_front, el_capital_front)
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
  
# ----------------------------------------------------------- #
# Figure A.5. First-order parameter estimates of the Translog #
# production frontier with 95% confidence intervals           #
# ----------------------------------------------------------- #

  # Combine all models
  models_cart <- bind_rows(front.orig, front.syn1, front.syn2, front.syn3,
                          front.syn4, front.syn5)
  
  models_normrank <- bind_rows(front.orig, front.syn6, front.syn7, front.syn8,
                           front.syn9, front.syn10)
  
  # Make plot for CART
  elast_TL_front_cart <- dwplot(models_cart, 
                          show_intercept = FALSE,
                          ci = .95,
                          vars_order = c("Land", "Labor", "Feed", 
                                         "Material", "Cattle", "Capital"),
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
    scale_x_continuous(breaks = seq(0, 0.4, by = 0.1), limits=c(-0.02,0.4)) +
    scale_colour_manual(values=cbPalette)
  
  # Make plot for NORMRANK
  elast_TL_front_normrank <- dwplot(models_normrank, 
                                show_intercept = FALSE,
                                ci = .95,
                                vars_order = c("Land", "Labor", "Feed", 
                                               "Material", "Cattle", "Capital"),
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
    scale_x_continuous(breaks = seq(0, 0.4, by = 0.1), limits=c(-0.02,0.4)) +
    scale_colour_manual(values=cbPalette)
  
  # Save plots
  elast_TL_front <- ggarrange(elast_TL_front_cart, elast_TL_front_normrank, 
                              ncol = 2, nrow = 1)
  ggsave("Figures/Norway/FigA5.elast_TL_front.png", elast_TL_front, 
         width = 23, height = 11, units = "cm")

# ----------------------------------------------------------------- #
# Figure A.6. Distribution of farm-level elasticities obtained from #
#  the Translog production frontier.                                #
# ----------------------------------------------------------------- #
  
  # Make data frame with elasticities from all datasets (original and synthetic)
  
    # Land
    nr_obs <- length(dat_orig$el_land_front)
    el_land <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_land) <- c("el", "data") 
    el_land$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                      rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                      rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                      rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                      rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                      rep("Synthetic data 10", nr_obs))
    el_land$el <- c(dat_orig$el_land_front,dat_syn1$el_land_front,dat_syn2$el_land_front,
                    dat_syn3$el_land_front,dat_syn4$el_land_front,dat_syn5$el_land_front,
                    dat_syn6$el_land_front,dat_syn7$el_land_front,dat_syn8$el_land_front,
                    dat_syn9$el_land_front,dat_syn10$el_land_front)
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
    el_labor$el <- c(dat_orig$el_labor_front,dat_syn1$el_labor_front,dat_syn2$el_labor_front,
                     dat_syn3$el_labor_front,dat_syn4$el_labor_front,dat_syn5$el_labor_front,
                     dat_syn6$el_labor_front,dat_syn7$el_labor_front,dat_syn8$el_labor_front,
                     dat_syn9$el_labor_front,dat_syn10$el_labor_front)
    el_labor$inp <- "Labor"
    
    # Feed
    el_feed <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_feed) <- c("el", "data") 
    el_feed$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                      rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                      rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                      rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                      rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                      rep("Synthetic data 10", nr_obs))
    el_feed$el <- c(dat_orig$el_feed_front,dat_syn1$el_feed_front,dat_syn2$el_feed_front,
                    dat_syn3$el_feed_front,dat_syn4$el_feed_front,dat_syn5$el_feed_front,
                    dat_syn6$el_feed_front,dat_syn7$el_feed_front,dat_syn8$el_feed_front,
                    dat_syn9$el_feed_front,dat_syn10$el_feed_front)
    el_feed$inp <- "Feed"
    
    # Material
    el_material <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_material) <- c("el", "data") 
    el_material$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                          rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                          rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                          rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                          rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                          rep("Synthetic data 10", nr_obs))
    el_material$el <- c(dat_orig$el_material_front,dat_syn1$el_material_front,dat_syn2$el_material_front,
                        dat_syn3$el_material_front,dat_syn4$el_material_front,dat_syn5$el_material_front,
                        dat_syn6$el_material_front,dat_syn7$el_material_front,dat_syn8$el_material_front,
                        dat_syn9$el_material_front,dat_syn10$el_material_front)
    el_material$inp <- "Material"
    
    # Cattle
    el_cattle <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_cattle) <- c("el", "data") 
    el_cattle$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                        rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                        rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                        rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                        rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                        rep("Synthetic data 10", nr_obs))
    el_cattle$el <- c(dat_orig$el_cattle_front,dat_syn1$el_cattle_front,dat_syn2$el_cattle_front,
                      dat_syn3$el_cattle_front,dat_syn4$el_cattle_front,dat_syn5$el_cattle_front,
                      dat_syn6$el_cattle_front,dat_syn7$el_cattle_front,dat_syn8$el_cattle_front,
                      dat_syn9$el_cattle_front,dat_syn10$el_cattle_front)
    el_cattle$inp <- "Cattle"
    
    # Capital
    el_capital <- data.frame(matrix(ncol = 2, nrow = 11*nr_obs))
    colnames(el_capital) <- c("el", "data") 
    el_capital$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                         rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                         rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs),
                         rep("Synthetic data 6", nr_obs),rep("Synthetic data 7", nr_obs),
                         rep("Synthetic data 8", nr_obs),rep("Synthetic data 9", nr_obs),
                         rep("Synthetic data 10", nr_obs))
    el_capital$el <- c(dat_orig$el_capital_front,dat_syn1$el_capital_front,dat_syn2$el_capital_front,
                       dat_syn3$el_capital_front,dat_syn4$el_capital_front,dat_syn5$el_capital_front,
                       dat_syn6$el_capital_front,dat_syn7$el_capital_front,dat_syn8$el_capital_front,
                       dat_syn9$el_capital_front,dat_syn10$el_capital_front)
    el_capital$inp <- "Capital"
    
    # All
    el_front <- rbind(el_land,el_labor,el_feed,el_material,el_cattle,el_capital)
    el_front_cart <- el_front %>% # for CART plot
      filter(data == "Original data" | data == "Synthetic data 1" |
               data == "Synthetic data 2" | data == "Synthetic data 3" |
               data == "Synthetic data 4" | data == "Synthetic data 5")
    el_front_normrank <- el_front %>% # for NORMRANK plot
      filter(data == "Original data" | data == "Synthetic data 6" |
               data == "Synthetic data 7" | data == "Synthetic data 8" |
               data == "Synthetic data 9" | data == "Synthetic data 10")
  
  # Make sure data are plotted in right order 
    
    el_front_cart$inp <- factor(el_front_cart$inp,
                               levels=c("Land","Labor","Feed",
                                        "Material","Cattle","Capital"))
    el_front_cart$data <- factor(el_front_cart$data,
                                levels=c("Synthetic data 5",
                                         "Synthetic data 4",
                                         "Synthetic data 3",
                                         "Synthetic data 2",
                                         "Synthetic data 1",
                                         "Original data"))
    el_front_normrank$inp <- factor(el_front_normrank$inp,
                              levels=c("Land","Labor","Feed",
                                       "Material","Cattle","Capital"))
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
    scale_y_continuous(breaks = seq(-0.3, 0.6, by = 0.3), limits=c(-0.5,0.6)) +
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
    scale_y_continuous(breaks = seq(-0.3, 0.6, by = 0.3), limits=c(-0.5,0.6)) +
    scale_fill_manual(values=cbPalette)         
  
  # Save plots
  elast_TL_fl_front <- ggarrange(elast_TL_fl_front_cart, elast_TL_fl_front_normrank, 
                                 ncol = 2, nrow = 1)
  ggsave("Figures/Norway/FigA6.elast_TL_fl_front.png", elast_TL_fl_front, 
         width = 23, height = 11, units = "cm")
  
# ---------------------------------------------------------------- #
# Table A.3. Technical efficiency scores, Translog functional form #
# ---------------------------------------------------------------- # 
  
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
            out = "Tables/Norway/TabA3.descr_te_tl.txt")

# ---------------------------------------------------------------- #
# Figure A.7. Distribution of technical efficiency scores based on # 
# the Translog production frontier.                                #
# ---------------------------------------------------------------- #
  
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
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.15,1)) + 
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
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.15,1)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 6", 
                                 "Synthetic data 7", 
                                 "Synthetic data 8", 
                                 "Synthetic data 9", 
                                 "Synthetic data 10")) 
  
  # Save plots
  te_tl <- ggarrange(te_tl_cart, te_tl_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/Norway/FigA7.te_tl.png", te_tl, 
         width = 23, height = 11, units = "cm")
  
# ------------------------------------ #
#### 5. Data envelopment analysis   ####
# ------------------------------------ #
  
  # Define inputs and outputs
    
    inputs_orig <- data.frame(dat_orig[c("x1","x2","x3","x4","x5","x6")])
    outputs_orig <- data.frame(dat_orig[c("y1")]) 
    inputs_syn1 <- data.frame(dat_syn1[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn1 <- data.frame(dat_syn1[c("y1")]) 
    inputs_syn2 <- data.frame(dat_syn2[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn2 <- data.frame(dat_syn2[c("y1")]) 
    inputs_syn3 <- data.frame(dat_syn3[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn3 <- data.frame(dat_syn3[c("y1")]) 
    inputs_syn4 <- data.frame(dat_syn4[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn4 <- data.frame(dat_syn4[c("y1")]) 
    inputs_syn5 <- data.frame(dat_syn5[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn5 <- data.frame(dat_syn5[c("y1")]) 
    inputs_syn6 <- data.frame(dat_syn6[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn6 <- data.frame(dat_syn6[c("y1")]) 
    inputs_syn7 <- data.frame(dat_syn7[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn7 <- data.frame(dat_syn7[c("y1")]) 
    inputs_syn8 <- data.frame(dat_syn8[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn8 <- data.frame(dat_syn8[c("y1")]) 
    inputs_syn9 <- data.frame(dat_syn9[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn9 <- data.frame(dat_syn9[c("y1")]) 
    inputs_syn10 <- data.frame(dat_syn10[c("x1","x2","x3","x4","x5","x6")])
    outputs_syn10 <- data.frame(dat_syn10[c("y1")]) 
  
  # Estimate output-oriented distance functions
  
    # Original data
    odf_orig <- dea(inputs_orig, outputs_orig, RTS = "vrs", ORIENTATION = "out")
    te_dea_orig <- 1/eff(odf_orig, asInData=TRUE) # For Shephard TE
    summary(te_dea_orig)
    odf_orig <- NULL
    
    # Synthetic data 1
    odf_syn1 <- dea(inputs_syn1, outputs_syn1, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn1 <- 1/eff(odf_syn1, asInData=TRUE) 
    summary(te_dea_syn1)
    odf_syn1 <- NULL
    
    # Synthetic data 2
    odf_syn2 <- dea(inputs_syn2, outputs_syn2, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn2 <- 1/eff(odf_syn2, asInData=TRUE) 
    summary(te_dea_syn2)
    odf_syn2 <- NULL
    
    # Synthetic data 3
    odf_syn3 <- dea(inputs_syn3, outputs_syn3, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn3 <- 1/eff(odf_syn3, asInData=TRUE) 
    summary(te_dea_syn3)
    odf_syn3 <- NULL
    
    # Synthetic data 4
    odf_syn4 <- dea(inputs_syn4, outputs_syn4, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn4 <- 1/eff(odf_syn4, asInData=TRUE) 
    summary(te_dea_syn4)
    odf_syn4 <- NULL
    
    # Synthetic data 5
    odf_syn5 <- dea(inputs_syn5, outputs_syn5, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn5 <- 1/eff(odf_syn5, asInData=TRUE) 
    summary(te_dea_syn5)
    odf_syn5 <- NULL
    
    # Synthetic data 6
    odf_syn6 <- dea(inputs_syn6, outputs_syn6, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn6 <- 1/eff(odf_syn6, asInData=TRUE) 
    summary(te_dea_syn6)
    odf_syn6 <- NULL
    
    # Synthetic data 7
    odf_syn7 <- dea(inputs_syn7, outputs_syn7, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn7 <- 1/eff(odf_syn7, asInData=TRUE) 
    summary(te_dea_syn7)
    odf_syn7 <- NULL
    
    # Synthetic data 8
    odf_syn8 <- dea(inputs_syn8, outputs_syn8, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn8 <- 1/eff(odf_syn8, asInData=TRUE) 
    summary(te_dea_syn8)
    odf_syn8 <- NULL
    
    # Synthetic data 9
    odf_syn9 <- dea(inputs_syn9, outputs_syn9, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn9 <- 1/eff(odf_syn9, asInData=TRUE) 
    summary(te_dea_syn9)
    odf_syn9 <- NULL
    
    # Synthetic data 10
    odf_syn10 <- dea(inputs_syn10, outputs_syn10, RTS = "vrs", ORIENTATION = "out")
    te_dea_syn10 <- 1/eff(odf_syn10, asInData=TRUE) 
    summary(te_dea_syn10)
    odf_syn10 <- NULL
    
# ----------------------------------------------------------- #
# Table A.4. Output-oriented technical efficiency scores, DEA #
# ----------------------------------------------------------- #
    
  # Combine all TE scores from DEA
  descr_dea <- cbind(te_dea_orig,te_dea_syn1,te_dea_syn2,te_dea_syn3,
                     te_dea_syn4,te_dea_syn5,te_dea_syn6,te_dea_syn7,
                     te_dea_syn8,te_dea_syn9,te_dea_syn10)
  
  # Create table with descriptive statistics
  stargazer(as.data.frame(descr_dea),
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
            out = "Tables/Norway/TabA4.descr_TE_dea.txt")
  
# -------------------------------------------------------------- #
# Figure 3. Distribution of output-oriented technical efficiency # 
#  scores based on data envelopment analysis (DEA)               #
# -------------------------------------------------------------- #
  
  # Make data frame with technical efficiency scores from all datasets
  
    # Original and CART data
    nr_obs <- length(te_dea_orig)
    efficiencies_cart <- data.frame(matrix(ncol = 2, nrow = 6*nr_obs))
    colnames(efficiencies_cart) <- c("TE", "data") 
    efficiencies_cart$data <- c(rep("Original data", nr_obs),rep("Synthetic data 1", nr_obs),
                              rep("Synthetic data 2", nr_obs),rep("Synthetic data 3", nr_obs),
                              rep("Synthetic data 4", nr_obs),rep("Synthetic data 5", nr_obs))
    efficiencies_cart$TE <- c(te_dea_orig,te_dea_syn1,te_dea_syn2,te_dea_syn3,
                              te_dea_syn4,te_dea_syn5)
    
    # Original and NORMRANK data
    efficiencies_normrank <- data.frame(matrix(ncol = 2, nrow = 6*nr_obs))
    colnames(efficiencies_normrank) <- c("TE", "data") 
    efficiencies_normrank$data <- c(rep("Original data", nr_obs),rep("Synthetic data 6", nr_obs),
                               rep("Synthetic data 7", nr_obs),rep("Synthetic data 8", nr_obs),
                               rep("Synthetic data 9", nr_obs),rep("Synthetic data 10", nr_obs))
    efficiencies_normrank$TE <- c(te_dea_orig,te_dea_syn6,te_dea_syn7,te_dea_syn8,
                                  te_dea_syn9,te_dea_syn10)
    
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
  te_dea_cart <-     
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
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.15,1)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 1", 
                                 "Synthetic data 2", 
                                 "Synthetic data 3", 
                                 "Synthetic data 4", 
                                 "Synthetic data 5"))
  
  # Make plot for NORMRANK
  te_dea_normrank <-     
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
    scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.10,1.01)) + 
    scale_fill_manual(values=cbPalette,
                      limits = c("Original data", 
                                 "Synthetic data 6", 
                                 "Synthetic data 7", 
                                 "Synthetic data 8", 
                                 "Synthetic data 9", 
                                 "Synthetic data 10")) 
  
  # Save plots
  te_dea <- ggarrange(te_dea_cart, te_dea_normrank, ncol = 2, nrow = 1)
  ggsave("Figures/Norway/Fig3.te_dea.png", te_dea, 
         width = 23, height = 11, units = "cm")
  
  # Same figures in black and white for journal
    
    # Make plot for CART
    te_dea_cart <-     
      efficiencies_cart %>%
      ggplot(aes(x = TE , y = data,  fill = data)) +
      ggridges::geom_density_ridges(scale = 1, 
                                    stat = "binline", 
                                    binwidth = 0.01) +
      theme_bw(base_family="Times") + 
      xlab("Technical efficiency") + ylab("") +
      ggtitle("CART method") +
      theme(
        legend.position = "bottom",
        legend.background = element_rect(colour = "white"),
        legend.title = element_blank(),
        plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y=element_blank() 
      ) +
      scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.10,1.01)) + 
      scale_fill_manual(values=bwPalette,
                        limits = c("Original data", 
                                   "Synthetic data 1", 
                                   "Synthetic data 2", 
                                   "Synthetic data 3", 
                                   "Synthetic data 4", 
                                   "Synthetic data 5"))
    
    # Make plot for NORMRANK
    te_dea_normrank <-     
      efficiencies_normrank %>%
      ggplot(aes(x = TE , y = data,  fill = data)) +
      ggridges::geom_density_ridges(scale = 1, 
                                    stat = "binline", 
                                    binwidth = 0.01) +
      theme_bw(base_family="Times") + 
      xlab("Technical efficiency") + ylab("") +
      ggtitle("NORMRANK method") +
      theme(
        legend.position = "bottom",
        legend.background = element_rect(colour = "white"),
        legend.title = element_blank(),
        plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y=element_blank() 
      ) +
      scale_x_continuous(breaks = seq(0.25, 1, by = 0.25), limits=c(0.15,1.01)) + 
      scale_fill_manual(values=bwPalette,
                        limits = c("Original data", 
                                   "Synthetic data 6", 
                                   "Synthetic data 7", 
                                   "Synthetic data 8", 
                                   "Synthetic data 9", 
                                   "Synthetic data 10")) 
    
    # Save plots
    te_dea <- ggarrange(te_dea_cart, te_dea_normrank, ncol = 2, nrow = 1)
    ggsave("Figures/Norway/Fig3.te_dea_bw.png", te_dea, 
           width = 23, height = 11, units = "cm")  

# --------------------------------------- #
#            End of the script            #
# --------------------------------------- #
    
    

# Load packages
library(DemoDecomp)
library(HMDHFDplus)
library(here)
library(colorspace)
library(tidyverse)

# Set smart working directory
i_am('TFR_from_parity.R')

# Set Human Fertility Database credentials.
# DO NOT LEAVE THEM HERE IF YOU SHARE THIS FILE
myHFDusername <- 'myHFDusername'
myHFDpassword <- 'myHFDpassword'

# Download conditional age-specific fertility rates by birth order for Sweden
# in 2000 and 2015
ASFRmi <- readHFDweb(
  'SWE',
  'mi',
  username = myHFDusername,
  password = myHFDpassword
  ) %>%
  filter(Year == 2000 | Year==2015)

# Split data into values for 2000 (period 1) and 2015 (period 2)
fxp1 <- ASFRmi %>%
  filter(Year == 2000) %>%
  pivot_longer(m1x:m5px,names_to = 'parity',values_to = 'fxp') %>%
  arrange(parity,Age) %>%
  pull(fxp)

fxp2 <- ASFRmi %>%
  filter(Year == 2015) %>%
  pivot_longer(m1x:m5px,names_to = 'parity',values_to = 'fxp') %>%
  arrange(parity,Age) %>%
  pull(fxp)

# Define function that takes as input a long vector of parity- and age-specific
# fertility rates and computes TFR. This is the tricky part but also what allows
# us to use general decomposition techniques such as the line-integral method,
# the stepwise-replacement method, and the life-table response experiment method
# as implemented in the DemoDecomp package to decompose TFR differences into 
# parity- and age-specific components.

# The flexible implementation below allows for any parity and any age range
compute_TFR_from_parity_mat <- function(fxp,dims) {
  
  fxp_mat <- matrix(fxp,nrow=dims[2],ncol=dims[1])
  prob_mat <- apply(fxp_mat,M=2,function(x) x/(1+0.5*x))
  prob_mat <- t(apply(prob_mat,M=1,function(x) if_else(is.na(x),0,x)))
  
  births <- matrix(0,nrow=dims[2],ncol=dims[1])
  pops <- matrix(0,nrow=dims[2],ncol=dims[1])
  
  pops[1,1] <- 1000
  
  # calculating the births and population by age / parity 
  # the 0.5 is because we assume the births happen halfway through the year on avg.
  for (a in 1:(dims[2]-1)) {
    
    births[a,1] <- pops[a,1] * prob_mat[a,1]
    pops[a+1,1] <- pops[a,1] - births[a,1]  
    
    for (p in 2:(dims[1]-1)) {
      
      births[a,p] <- sum(pops[a,p],births[a,p-1]*0.5, na.rm=T) * prob_mat[a,p]
      births[a,p][is.na(births[a,p])] <- 0
      pops[a+1,p] <- pops[a,p] + births[a,p-1] - births[a,p]
      
    }
    
    births[a,dims[1]] <- sum(pops[a,dims[1]],births[a,dims[1]-1]*0.5, na.rm=T) * prob_mat[a,dims[1]]
    births[a,dims[1]][is.na(births[a,dims[1]])] <- 0
    pops[a+1,dims[1]] <- pops[a,dims[1]] + births[a,dims[1]-1]
    
  }
  
  # calculating the TFR as a sum of the births by age & parity
  TFR <- sum(births) / 1000
  
  TFR
}


# Run decomposition (stepwise)
decomp <- tibble(
  Age = rep(unique(ASFRmi$Age),5),
  parity = rep(1:5,each=44),
  pars1 = fxp1,
  pars2 = fxp2,
  C = stepwise_replacement(
    func = compute_TFR_from_parity_mat,
    pars1 = pars1,
    pars2 = pars2,
    dims = c(5,44)
  )
)

# Visualize the age- and parity- specific contributions to TFR change
decomp_plot <- decomp %>%
  mutate(age=fct_reorder(as.character(Age),Age)) %>%
  ggplot() +
  geom_bar(aes(x=age, y=C, fill=as.factor(parity)),stat="identity") +
  scale_fill_discrete_sequential(palette='SunsetDark') +
  labs(
    x='Age',
    y='Contribution to TFR Change',
    title='Decomposition of TFR Change 2000-2015 Sweden by Age and Parity',
    fill='Parity'
    ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
    )

decomp_plot

# Save graph
svg(here('decomp_results.svg'),width=9,height = 6)
decomp_plot
dev.off()

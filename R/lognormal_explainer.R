# Sample data
library(tidyverse)
theme_set(theme_classic())


host <- c(1.1, 1.87, 2.3, 15, 25, 50.4)  # converted to days
#soils <- c(0.44, 50*216, 0.0903, 7.5) # double check this
soil <- c(0.44, 1*24, 160*24,6*24)
sw <- c(1.4*24, 231*24)
subsurf <- 24*365*c(3, 1000)


#
# creates a lognormal distribution from a vector of x (independant) values, mu, and sigma
my_lognormal <- function(x, mu, sigma) {
  (1/(x*sigma*sqrt(2*pi)))*exp(-1*(log(x)-mu)^2/(2*sigma^2))
}

# fits a lognormal distribution to a vector of raw data
fit_lognormal <- function(raw.data, grid) {
  ln_fit <- MASS::fitdistr(raw.data, "log-normal")
  mu <- ln_fit$estimate["meanlog"]
  sigma <- ln_fit$estimate["sdlog"]
  
  density.values <- my_lognormal(x=grid, mu=mu, sigma=sigma)
  density_df <- tibble(x=grid,
                       density = density.values)
  density_df
}

# Note that I want the probability that any given point is at doubling_time=x 
# AND abdundance = y
# p(a,b) = p(a)*p(b)
# fit_lognormal gives me probability densities functions for one variable
# I think I just need to calculate a prob density for each variable (doubling_time and abundance)
# and multiply them together

# host_doubling_density <- fit_lognormal(raw.data = host, 
#                               grid = seq(from=0.1, to=100, by=0.1)) |> 
#   rename(c(doubling.time="x", doubling.density="density"))
# host_abundance_density <- fit_lognormal(raw.data = c(0.01, 0.011, 0.009), # this is bullshit made up numbers
#                                         grid = seq(from=0.001, to=1, by=0.001)) |> 
#   rename(c(rel.abund="x", rel.abund.density="density"))
# joint_density <- expand_grid(host_doubling_density, 
#                              host_abundance_density, .name_repair = "unique") |> 
#   mutate(joint.density = doubling.density * rel.abund.density, 
#          norm.joint.density = joint.density/max(joint.density, na.rm=TRUE)) # this works!!!



calc_joint_density <- function(doubling.data, 
                               doubling.grid = 10^(seq(from=log10(0.1), to=log10(87600000), length.out=100)), 
                               rel.abund.data, 
                               rel.abund.grid = seq(from=0.01, to=1, by=0.01),
                               envt = NULL) {
  # calc lognormal distributions for doubling times and relative abundances
  doubling_density <- fit_lognormal(raw.data = doubling.data, 
                                    grid = doubling.grid) |> 
    rename(c(doubling.time="x", doubling.density="density"))
  abundance_density <- fit_lognormal(raw.data = rel.abund.data, # this is bullshit made up numbers
                                     grid = rel.abund.grid) |> 
    rename(c(rel.abund="x", rel.abund.density="density"))
  
  # Calculate a joint density distribution from each distribution
  joint_density <- expand_grid(doubling_density, 
                               abundance_density) |> 
    mutate(joint.density = doubling.density * rel.abund.density, 
           norm.joint.density = joint.density/max(joint.density, na.rm=TRUE)) # this works!!!
  
  # Add a column for environment if need be
  if(!is.null(envt)) {
    joint_density <- joint_density |> 
      mutate(envt = as.character(envt))
  }
  
  joint_density
}



# What if I define a universal grid for all situations?
# rel abundance is easy - it shoudl got from 0 to 1, by 0.01 or 0.001, linearly
# Doubling time should be log-spaced - so it should go from 10^-1 hours to 
#    10,000 year * 365 days/year * 24 hours/day = 87,600,000
#    But it should be log spaced, so it should go evenly from log10(0.1) to log10(87,600,000)
#doubling.grid.univ <- 10^(seq(from=log10(0.1), to=log10(87600000), length.out=100))
#rel.abund.grid.univ <- seq(from=0.01, to=1, by=0.01)

joint_density_host <- calc_joint_density(doubling.data = host,
                                         rel.abund.data = c(0.01, 0.02, 0.03),
                                         envt="host")
joint_density_soil <- calc_joint_density(doubling.data = soil,
                                         rel.abund.data = c(0.1, 0.2, 0.4),
                                         envt="soil")
joint_density_sw <- calc_joint_density(doubling.data = sw,
                                       rel.abund.data = c(0.5, 0.6, 0.7),
                                       envt="sw")

densities <- rbind(joint_density_host, 
                   joint_density_soil,
                   joint_density_sw)



# What if I just sample from the lognormal distribution
n.points <- 10000
set.seed(2112)
soil.double.samps <- 10^rnorm(n=n.points, mean=mean(log10(soil)), sd=sd(log10(soil)/3))
sw.double.samps <- 10^rnorm(n=n.points, mean=mean(log10(sw)), sd=sd(log10(sw)/3))
host.double.samps <- 10^rnorm(n=n.points, mean=mean(log10(host)), sd=sd(log10(host)/3))
subsurf.double.samps <- 10^rnorm(n=n.points, mean=mean(log10(subsurf)), sd=sd(log10(subsurf)/3))

soil.abund.samps <- 10^rnorm(n=n.points, mean=log10(7.5), sd=log10(7.5)/10) # fourfold
sw.abund.samps <- 10^rnorm(n=n.points, mean=log10(1.6), sd=log10(1.6)/10) # twofold uncertainty
host.abund.samps <- 10^rnorm(n=n.points, mean=log10(0.01), sd=abs(log10(0.01)/10))
subsurf.abund.samps <- 10^rnorm(n=n.points, mean=log10(74), sd=log10(1.5)) # terrestrial subsurface is 20fold, marine susburface is 8fold

# OK so now I've got samples of all the environment doubling times and abundances
soils_df <- tibble(double=soil.double.samps, abund=soil.abund.samps, envt="soil")
sw_df <- tibble(double=sw.double.samps, abund=sw.abund.samps, envt="seawater")
host_df <- tibble(double=host.double.samps, abund=host.abund.samps, envt="host")
subsurf_df <- tibble(double=subsurf.double.samps, abund=subsurf.abund.samps, envt="subsurface" )
all <- rbind(soils_df, sw_df, host_df, subsurf_df) |> 
  slice_sample(prop = 1)

# Make labels for ggrepel; it's repetitive but I don't care
blob_labels <- all |> 
  group_by(envt) |> 
  summarise(mean.dt = mean(double, na.rm=TRUE),
            mean.abund = mean(abund, na.rm=TRUE))


# setup for plot
breaks.vec <- c(1, 24, 24*31, 24*365.25, 24*365.25*10, 24*365*100, 24*365*1000, 24*365*10000)
breaks.names <- c("1 hr", "1 day", "1 month", "1 yr", "10 yr", "100 yr", "1000 yr", "10,000 yr")
colors <- c("#908C13", "#2B5597", "#990000","#674422")

# the plot
ggplot(all, aes(x=double, y=abund)) + 
  stat_density_2d(
    aes(fill = envt, alpha = after_stat(level)), 
    geom = "polygon",
    color = NA,
    bins = 500
  ) + 
  scale_x_log10(name = "doubling/turnover time", 
                breaks = breaks.vec,
                labels = breaks.names) + 
  ggrepel::geom_text_repel(data=blob_labels, 
                            aes(x=mean.dt, y=mean.abund, label=envt)) +
  scale_y_log10(name = "biomass, Gt C") +
  scale_fill_manual(values = colors) + 
  theme(text = element_text(size = 10),
        legend.position = "none",
        axis.text.x = element_text(angle=-45, hjust=0),
        plot.margin = margin(0.1,0.25,0.1,0.1, "in"))
ggsave("plots/fig_1_v2.png", height = 2.5, width=3.5, units = "in", dpi=300)

########
# All belwo is BS I think
########

  


ggplot(joint_density_host ,
       , aes(x=doubling.time, y=rel.abund)) +
  geom_tile(aes(alpha=norm.joint.density), fill="black", color=NA) +
  scale_x_log10() + 
  theme_classic( )
 ggsave("~/Downloads/eatme.png", height = 30, width = 30, dpi = 100)   
  


ggplot() + 
  geom_line(data=test_df, aes(x=x, y=dens)) + 
  geom_rug(data = raw_data, aes(x=raw.data))

# Calculate parameters

mean <- mean(raw.data) 

sd <- sd(raw.data)



# Accessing the lognormal density function

density.values <- dlnorm(x = 1:100, 
                         meanlog = meanlog, sdlog = sdlog) 


d_test <- tibble(density.values)
ggplot(d_test, aes(x=density.values)) + 
  geom_density() + 
  geom_rug(data = raw_data, aes(x=raw.data)) +
  scale_x_log10() 


# Plot the lognormal distribution

plot(1:100, 
     density.values, type = "l", xlab = "Value", ylab = "Density") 



### Trying it with normal distributuion
WOuldn't it just be easier to take the log of the data, calculate mean and median, 
# plot the resulting distributions, and don't log-transform the graph?
  # I think so
  
  log.raw.data <- log10(raw.data)
log.mean.raw <- mean(log.raw.data)
log.sd.raw <- sd(log.raw.data)

my_normal <- function(x, mu, sigma) {
  (1/(x*sigma*sqrt(2*pi)))*exp(-1*(x-mu)^2/(2*sigma^2))
}

fit_normal <- function(raw.data, grid) {
  mu <- mean(raw.data)
  sigma <- sd(raw.data)
  norm_fit <- MASS::fitdistr(x=raw.data, densfun="normal", mu=mu, sigma=sigma)
  
  
}
fit_normal(host, grid=seq(0.1, 100, by=0.1)) 

fit_normal(raw.data)

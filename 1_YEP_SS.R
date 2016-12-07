require(R2jags)
require(doBy)
# rm(list=ls())


################
# NOTE: IF USED- MUST STANDARDIZE COVARIATES BEFORE PUTTING THEM IN A MATRIX
# THE SCALE FUNCTION PERFORMS COLUMN-WISE STANDARDIZATION TO A MATRIX - NOT WHAT WE WANT




########################################################
# Format data
dat <- read.csv("MN_GN_CPUE_all_lakes.csv") 
head(dat)
length(unique(dat$DOW))

temp <- read.csv("temperature_data_for_ty.csv") 
head(temp)

wqdat <- read.csv("wq_data_for_ty.csv")
head(wqdat)

## Clean data
# Select out lakes with n >= X
tbl1 <- table(dat$DOW)
dat <- dat[dat$DOW %in% names(tbl1)[tbl1 >=8],]
length(unique(dat$DOW))
head(dat)
range(dat$Year)

# Select temp data for lakes in dat
temp <- temp[temp$DOW %in% dat$DOW,]
length(unique(temp$DOW))

# Convert ice-off date to decimal day of year
temp$ice_off_date <- as.POSIXct(temp$ice_off_date, format="%m/%d/%Y")
temp$ice_off_julian <- as.numeric(strftime(temp$ice_off_date, format = "%j"))
head(temp)

# Not all lakes in dat have temp, so remove lakes in dat that are not in temp
dat <- dat[dat$DOW %in% temp$DOW,]
length(unique(dat$DOW))

# Grab water quality data for lakes in dat
wqdat <- wqdat[wqdat$DOW %in% dat$DOW,]
length(unique(dat$DOW))

# Look at years contained in all datasets - match them up
range(dat$Year)
range(temp$Year)
range(wqdat$SAMPLE_YEAR)

# Remove sample year 2016 from dat - other data sets only go to 2015
# and start all data sets at 1987 (latest date for cpe data)
dat <- dat[dat$Year < 2016, ]
range(dat$Year)

temp <- temp[temp$Year > 1986, ]
range(temp$Year)

wqdat <- wqdat[wqdat$SAMPLE_YEAR > 1986, ]
range(wqdat$SAMPLE_YEAR)

# Sort data frames by DOW and Year for simplicity
dat <- dat[order(dat$DOW, dat$Year), ]
temp <- temp[order(temp$DOW, temp$Year), ]
wqdat <- wqdat[order(wqdat$DOW, wqdat$SAMPLE_YEAR), ]

# max(table(dat$DOW, dat$Year))

# Fill in missing years for all lakes with NA in dat
datY <- expand.grid(DOW = unique(dat$DOW), Year = min(dat$Year):max(dat$Year))
head(datY)

# Grab YEP data
yep <- dat[c("DOW", "Year", "YEP")]

# Merge yep with datY that has missing data for years with no data for each lake
yep <- merge(yep,datY,all=TRUE)
dim(yep)
head(yep)
length(unique(yep$Year))

# Grab WAE data
wae <- dat[c("DOW", "Year", "WAE")]
wae <- merge(wae,datY,all=TRUE)
dim(wae)
head(wae)
length(unique(wae$Year))

# Convert to matrices [t, i]
nSites <- length(unique(yep$DOW))
nYears <- length(unique(yep$Year))
# Matrices of CPE data
# Response variable
yepM <- matrix(yep$YEP, nrow=nYears, ncol=nSites, byrow=F)
# Log-tranform CPE
yepM <- log(yepM + 0.1)
# Predictor
waeM <- matrix(wae$WAE, nrow=nYears, ncol=nSites, byrow=F)
# Log-tranform CPE
waeM <- log(waeM + 0.1)
waeM[is.na(waeM)] <- 0 # assuming mean values for unsampled years

# Create matrices of temperature and water quality data
head(temp)
dd <- matrix(as.numeric(scale(temp$gdd_wtr_5c)), nrow=nYears, ncol=nSites, byrow=F)
sum(is.na(dd)) # Check for any missing values (NAs)
dim(dd)

iceoff <- matrix(as.numeric(scale(temp$ice_off_julian)), nrow=nYears, ncol=nSites, byrow=F)
sum(is.na(iceoff))

iceduration <- matrix(as.numeric(scale(temp$ice_duration_days)), nrow=nYears, ncol=nSites, byrow=F)
sum(is.na(iceduration))

coefvar <- matrix(as.numeric(scale(temp$coef_var_30.60)), nrow=nYears, ncol=nSites, byrow=F)
sum(is.na(coefvar))

cor(as.numeric(scale(temp$gdd_wtr_5c)), as.numeric(scale(temp$coef_var_30.60)))

# Create lake-level covariate vectors [i]: elevation, basin size, gradient
area <- as.numeric(by(dat$LAKE_AREA_GIS_ACRES, dat$DOW, mean))
area <- log(area)
area <- as.numeric(scale(area))
depth <- as.numeric(by(dat$MAX_DEPTH_FEET, dat$DOW, mean))
depth <- log(depth)
depth <- as.numeric(scale(depth))


########################################################
sink("SSModel1.txt")
cat("
  model {
    # Likelihood
    # State process

    for (t in 1:(T-1)){       #loop over year
      for (i in 1:S){         #loop over site
        r.proc[t,i] ~ dnorm(0 , tau.proc)
          r[t,i] <- int.r + r.proc[t,i] + beta[1] * dd[t,i] + beta[2] * coefvar[t,i] + 
                    beta[3] * wae[t, i] + beta[4] * depth[i] + beta[5] * dd[t,i] * depth[i]

        # Expoenetial growth rate model (N_t+1 = N_t * e^r)
        # Discrete-time growth rate model (log-linear model)
          logN[(t+1),i] <- logN[t,i] + r[t,i]
      }
    }
    # Observation process
    for (t in 1:T) {
      for (i in 1:S){
        y[t,i] ~ dnorm(logN[t,i], tau.obs)
      }
    }
    # Population sizes on real scale
    for (t in 1:T){
      for (i in 1:S){
        N[t,i] <- exp(logN[t,i])
      }
    }
    # Priors and constraints
    for (i in 1:S){
      logN[1,i] ~ dnorm(0,0.01)    # Prior for initial population size at every site                                 
    }
    int.r ~ dnorm(0, 0.1)                # Prior for mean growth rate
    
    for (a in 1:nBeta){
      beta[a] ~ dnorm(0, 0.01)    #Prior for climate covs
    }
    
    sigma.proc ~ dunif(0, 5)             # Prior for sd of state process
    tau.proc <- pow(sigma.proc, -2)      #Precision process
    
    sigma.obs ~ dunif(0, 5)              # Prior for sd of observation process
    tau.obs <- pow(sigma.obs, -2)        #Precision observation

    

}
    ",fill = TRUE)
sink()


########################################################
# Compile data
data <- list(y = yepM, # YEP CPE [t, i]
             dd = dd, # DD [t, i]
             coefvar = coefvar,
             wae = waeM,
             T = nrow(yepM),
             S = ncol(yepM),
             depth = depth,
             nBeta=5) 

########################################################
# Initial values
inits <- function() {
  list(sigma.proc = runif(1),
       sigma.obs =runif(1)
       )
}


# Parameters monitored
parameters <- c("int.r", "sigma.obs", "sigma.proc", "beta")

# MCMC settings
ni <- 4000
nt <- 2
nb <- 2000
nc <- 3

start.time = Sys.time()         # Set timer (2.3 mins)
# Call JAGS from R 

out <- jags(data, inits, parameters, "SSModel1.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Find which parameters, if any, have Rhat > 1.1
which(out$BUGSoutput$summary[, c("Rhat")] > 1.1)

# See what max Rhat value is
# max(out$BUGSoutput$summary[, c("Rhat")])


# Summarize posteriors
print(out, dig = 3)

hist(out$BUGSoutput$sims.list$beta[,2], col='gray')
mean(out$BUGSoutput$sims.list$beta[,2])
quantile(out$BUGSoutput$sims.list$beta[,2], c(0.025, 0.975))


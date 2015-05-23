# use your own directory, if not already set
# setwd('/Users/fabiodelia/Dropbox/Studium CSE/Statistische Physik und Computer Simulation/mdatom/sources')

dat = read.csv('resultsFor.txt', sep='', header = FALSE)
# dat$V1 is step
# dat$V9 is average energy at this step

# V9 is read as factor... convert
dat[,9] <- as.numeric(as.character(dat[,9]))
# remove cases with errors (NaN, inf etc)
dat = dat[dat$V9>1,]

plot(dat$V1, dat$V9, type='l', xlab='Steps', ylab='av. Energy', col='blue')
title('Maxwellian thermalisation')
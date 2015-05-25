# use your own directory, if not already set
setwd('/Users/fabiodelia/Dropbox/Studium CSE/Statistische Physik und Computer Simulation/mdatom/sources/results')

dat1 = read.csv('resultsFor1.txt', sep='', header = FALSE)
dat3 = read.csv('resultsFor3.txt', sep='', header = FALSE)
# dat$V1 is step
# dat$V9 is average energy at this step

# V9 is read as factor... convert
dat1[,9] <- as.numeric(as.character(dat1[,9]))
dat3[,9] <- as.numeric(as.character(dat3[,9]))

# remove cases with errors (NaN, inf etc)
dat1 <- do.call(data.frame,lapply(dat1, function(x) replace(x, is.infinite(x),0)))
dat3 <- do.call(data.frame,lapply(dat3, function(x) replace(x, is.infinite(x),0)))
dat1 = dat1[dat1$V9>0,]
dat3 = dat3[dat3$V9>0,]

xrange = range(0,10000)
yrange = range(min(dat1$V9, dat3$V9), max(dat1$V9, dat3$V9))

plot(xrange, yrange, type='n', xlab='Steps', ylab='E_tot')
lines(dat1$V1, dat1$V9, type='l', col='blue')
lines(dat3$V1, dat3$V9, type='l', col='green')
legend('topright', inset=.05, title='coupling', c('collisions', 'pressure scaling'), fill=c('blue', 'green'))

title(main='Coupling to a temperature or pressure bath, 2 atoms')


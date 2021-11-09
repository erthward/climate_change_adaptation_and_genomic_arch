library(movMF) # von Mises mix dists
#library(rotations) # has von Mises CDF fn
#library(circular) # has fns for working with circular data

# params to control behavior
plot.circ = F

# TODO:

        # figure out circ KDE plot

        # decide if/how to do KS testing...

        # figure out how to call R function from py script on savio, get params back
        # or if to just do after all sims are done instead

par(mfrow=c(3,2))

# read data from some of my simulation output
data <- read.csv('../output/output/output_PID-182509_DIR_short.csv')

# fn to read the data for a given set of values of the columns
get.subdf.data <- function(df, genicity, linkage, neutrality, nullness, it){
    subdf.data <- data[data$genicity == genicity &
                 data$linkage == linkage &
                 data$neutrality == neutrality &
                 data$nullness == nullness &
                 data$it == it, ]$dir
    # drop NAs
    subdf.data <- subdf.data[!is.na(subdf.data)]
    return(subdf.data)
}



non <- get.subdf.data(df, 100, 'independent', 'nonneut', 'non-null', 0)
null <- get.subdf.data(df, 100, 'independent', 'nonneut', 'null', 0)

hist(null, breaks=200, main='null')
hist(non, breaks=200, main='non')

# express angs (in degs) as x,y coords on the unit circle (for movMF)
ang.2.unit.circ.coord <- function(ang){
    x <- cos(ang/360*2*pi)
    y <- sin(ang/360*2*pi)
    return(c(x,y))
}
angs.2.unit.circ.coords <- function(angs){
    coords <- t(apply(as.matrix(angs), MARGIN=1, FUN=ang.2.unit.circ.coord))
    return(coords)
}
non.coords <- angs.2.unit.circ.coords(non)
null.coords <- angs.2.unit.circ.coords(null)


# fn to rotate [-180,180] data to [0,360] or vice versa
rotate.angs <- function(angs, dir='R'){
    rotated = c()
    if (dir == 'R'){ # rotate 'right'
        for (ang in angs){
            if (ang < 0){
                ang = ang + 360
            }
            rotated = c(rotated, ang)
        }
    } else if (dir == 'L'){ # rotate 'left'
        for (ang in angs){
            if (ang >180){
                ang = ang - 360
            }
            rotated = c(rotated, ang)
        }
    }
    return(rotated)
}


# fn to convert unit-circ coords back to angs (in degs)
unit.circ.coords.2.angs <- function(coords){
    angs <- atan2(coords[,2], coords[,1])*360/2/pi
    angs = rotate.angs(angs, 'R')
    return(angs)
}


# fn to plot data as a circular KDE
# taken from: https://www.r-bloggers.com/2011/03/circular-or-spherical-data-and-density-estimation/
# great food for thought, as I could either overplot and color by linkage/genicity
# or overplot and color areas of under/overrepresentation against my simulated null
plot.circular.kde <- function(angs){
    rad.angs = circular(rotate.angs(angs, 'L')/360*2*pi, type='angle', units='radians', rotation='counter')
    circ.dens = density(rad.angs+3*pi/2,bw=20)
    plot(rad.angs, stack=TRUE, shrink=.35, cex=0, sep=0.0,
         axes=FALSE,tol=.8,zero=c(0),bins=24,
         xlim=c(-2,2),ylim=c(-2,2), ticks=TRUE, tcl=.075)
    lines(circ.dens, col="darkgrey", lwd=3)
    text(0,0.8,"N", cex=2); text(0,-0.8,"S",cex=2); 
    text(0.8,0,"E",cex=2); text(-0.8,0,"W",cex=2)
}

# fit 4-part von Mises mixture dists
# for null
d <- movMF(null.coords, 4)
mu <- atan2(d$theta[,2], d$theta[,1])
kappa <- sqrt(rowSums(d$theta^2))
theta <- d$theta
alpha <- d$alpha

# draw random variates from the fitted dist and convert back to angs in degrees
null.rand <- rmovMF(10000, theta, alpha)
hist(unit.circ.coords.2.angs(null.rand), breaks=200, main='null, fitted')

# and for non-null
d <- movMF(non.coords, 4)
mu <- atan2(d$theta[,2], d$theta[,1])
kappa <- sqrt(rowSums(d$theta^2))
theta <- d$theta
alpha <- d$alpha

# draw random variates from the fitted dist and convert back to angs in degrees
non.rand <- rmovMF(10000, theta, alpha)
hist(unit.circ.coords.2.angs(non.rand), breaks=200, main='non, fitted')


# plot circular KDEs of the random variates
if (plot.circ){
    plot.circular.kde(null.rand)
    plot.circular.kde(non.rand)
}



##############################
##############################

# a main function, to be called from Python, return params of 4-dist vM mixture dist

fit.vM.mix.dist <- function(angs, nullness, n.mix=4, plot.it=F, plot.circ=F){
    # drop NAs
    angs <- angs[!is.na(angs)]

    # convert angs (in degrees) to unit-circle coords
    coords <- angs.2.unit.circ.coords(angs)

    # fit n.mix-part von Mises mixture dists
    d <- movMF(coords, n.mix)
    mu <- atan2(d$theta[,2], d$theta[,1])
    kappa <- sqrt(rowSums(d$theta^2))
    theta <- d$theta
    alpha <- d$alpha

    # plot fitted dist, if requested
    if (plot.it){ 
        if (nullness=='null'){
            color = '#79c2d9' # blue
        } else {
            color = '#c4626e' # red
        }
        # draw random variates from the fitted dist and convert back to angs in degrees
        rand <- rmovMF(10000, theta, alpha)
        hist(unit.circ.coords.2.angs(rand), breaks=200, main='null, fitted')
        if (plot.circ){
            # plot circular KDE of the random variates
            plot.circular.kde(rand)
        }
    }

    # return params
    params = list(mu=mu, kappa=kappa, alpha=alpha)
    return(params)
}



# get list of unique values for each of the non-data columns in the df
uniques <- apply(data[,!(colnames(data) %in% c("dir"))], 2, unique)

# set up an output list to store results
output = list('genicity'=c(), 'linkage'=c(), 'nullness'=c(), 'neutrality'=c(), 'it'=c(),
              'mu.1'=c(), 'mu.2'=c(), 'mu.3'=c(), 'mu.4'=c(),
              'kappa.1'=c(), 'kappa.2'=c(), 'kappa.3'=c(), 'kappa.4'=c(),
              'alpha.1'=c(), 'alpha.2'=c(), 'alpha.3'=c(), 'alpha.4'=c())

# loop over all scenarios
for (genicity in uniques[['genicity']]){
    for (linkage in uniques[['linkage']]){
        for (nullness in uniques[['nullness']]){
            for (neutrality in uniques[['neutrality']]){
                for (it in uniques[['it']]){
                    angs <- get.subdf.data(data, genicity, linkage, neutrality, nullness, it)
                    # only analyze if there are multiple rows' worth of data
                    if (length(angs) > 1){
                        params <- fit.vM.mix.dist(angs, nullness, n.mix=4)
                        # save all the data
                        for (param in names(params)){
                            for (num in seq(1,4)){
                                item.name = paste(param, num, sep='.')
                                output[[item.name]] <- c(output[[item.name]], params[[param]][num])
                            }
                        }
                        output[['genicity']] <- c(output[['genicity']], genicity) 
                        output[['linkage']] <- c(output[['linkage']], linkage) 
                        output[['nullness']] <- c(output[['nullness']], nullness) 
                        output[['neutrality']] <- c(output[['neutrality']], neutrality) 
                        output[['it']] <- c(output[['it']], it) 
                    }
                }
            }
        }
    }
}


# convert output to df
output.df = data.frame(output)

# write output.df to disk
write.csv(output.df, 'TEST_output.csv', row.names=F)


##############################
##############################



# function to calculate the CDF of a von Mises mixture dist
calc.CDF.movMF <- function(kappa, alpha){
    # vector of quantiles
    qs = seq(0, 100, 0.01)
    cdfs = list()
    for (i in seq(length(mu))){
        cdf = c()
        for (q in qs){
            # TODO: WTF TO DO ABOUT FACT THAT THERE'S NO MU?
            p = pvmises(q, kappa = kappa[i], lower.tail = TRUE)
            cdf = c(cdf, p)
        }
        cdfs[[i]] = cdf
    }
    #TODO: HOW TO SUM ARBITRARY-LENGTH LIST OF EQUAL-LENGTH VECTORS? (just loop it)
    #TODO: RETURN SUM
}


# TODO: FUNCTION TO CALC KS TEST STAT FOR A GIVEN RESULT

# TODO: FN TO RUN KS TESTS FOR 1-4 MIX DISTS, RETURN RESULTS OF BEST FITTING

# TODO: RUN FOR BOTH






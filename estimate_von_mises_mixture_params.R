library(movMF) # von Mises mix dists
library(rotations) # has von Mises CDF fn

# TODO:

        # figure out circ KDE plot

        # decide if/how to do KS testing...

        # figure out how to call R function from py script on savio, get params back
        # or if to just do after all sims are done instead

par(mfrow=c(3,2))

# read data from some of my simulation output
data <- read.csv('../output/output/output_PID-182509_DIR_short.csv')
non <- data[data$genicity == 100 &
           data$linkage == 'independent' &
           data$neutrality == 'nonneut' &
           data$nullness == 'non-null' &
           data$it == 0, ]$dir
null <- data[data$genicity == 100 &
           data$linkage == 'independent' &
           data$neutrality == 'nonneut' &
           data$nullness == 'null' &
           data$it == 0, ]$dir

# drop NAs
non <- non[!is.na(non)]
null <- null[!is.na(null)]

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
rand <- rmovMF(10000, theta, alpha)
hist(unit.circ.coords.2.angs(rand), breaks=200, main='null, fitted')
# plot circular KDE of the random variates
plot.circular.kde(rand)

# and for non-null
d <- movMF(non.coords, 4)
mu <- atan2(d$theta[,2], d$theta[,1])
kappa <- sqrt(rowSums(d$theta^2))
theta <- d$theta
alpha <- d$alpha

# draw random variates from the fitted dist and convert back to angs in degrees
rand <- rmovMF(10000, theta, alpha)
hist(unit.circ.coords.2.angs(rand), breaks=200, main='non, fitted')
# plot circular KDE of the random variates
plot.circular.kde(rand)



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






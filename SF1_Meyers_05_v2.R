library(scales)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)

#model parameters
N = 100 #number of individuals
u = 0.01 #mutation rate
s = 1 #selection coefficient (>0)
k = 0.5 #hedge factor (0 < k < 1)
lambdas = round(10^(seq(0,4,0.2))) #periods of Env oscillation
gens = 2e+4 - 1
En = 0

#genotype fitnesses in EnvA
wAg0 = 1+s
wAg1 = 1+s
wAg2 = 1+s
wAg3 = 1+s*k
wAg4 = 1
#as a vector
wAs = c(wAg0, wAg1, wAg2, wAg3, wAg4)

#genotype fitnesses in EnvB
wBg0 = 1
wBg1 = 1
wBg2 = 1
wBg3 = 1+s*k
wBg4 = 1+s
#as a vector
wBs = c(wBg0, wBg1, wBg2, wBg3, wBg4)
    
#a vector of initial genotype frequencies
g = c(0.2,0.2,0.2,0.2,0.2)
    
#the difference equation for updating genotype freq
gi_prime = function(Env, gis){
    #arguments are Env integer, and current 
    #allele frequency vector
    wis = wAs
    if (Env == 1){
        wis = wBs
    }
    griz = numeric(length=0)
    for (allele in 1:5){
        left = allele - 1
        right = allele + 1
        if (allele == 1){
            left = 5
        }
        if (allele == 5){
            right = 1
        }
        gip = gis[allele]*wis[allele]*(1-u) + (u/2)*(gis[left]*wis[left] + gis[right]*wis[right])
        griz = c(griz, gip)
    }
    return(griz/sum(griz))
}

freqs = matrix(,nrow=0,ncol=7)
colnames(freqs) <- c("g0","g1","g2","g3","g4","Env","Lambda")
for (lambda in lambdas){
    print(lambda)
    gis = g #start with equal frequency for each allele
    for (i in 1:gens){ #for each generation
        if (i %% lambda == 0){ #if divisible by lambda, switch the environment
            En = abs(En-1)
        }
        gp = gi_prime(En, gis) #given the environment and current feqs, update freqs
        freqs = rbind(freqs,c(gp,En,lambda)) #add data to the freqs matrix
        gis = gp #update the current freq vector
    }
}

df <- data.frame(freqs)
X <- df %>% group_by(Lambda, Env)
d <- X %>% summarize(g0=mean(g0), g1=mean(g1), g2=mean(g2), g3=mean(g3), g4=mean(g4))
dA <- d %>% filter(Env==0)
dB <- d %>% filter(Env==1)
dd <- dA %>% gather(genotype, frequency, g0:g4)
p1 = ggplot(dd, aes(x=Lambda, y=frequency, color=genotype))+
    geom_line(alpha=0.6, size=1.5)+
    scale_color_manual(values=c("mediumorchid2", "skyblue4", "seagreen3", "goldenrod2", "cyan4"))+
    scale_x_continuous(trans='log10') +
    ylab("genotype frequency")+
    xlab(expression(lambda))+
    theme_minimal()
print(p1)
    


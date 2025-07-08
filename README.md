# MultiplexSuppressionDriveResistanceSims
Stochastic simulations of multiplex suppression drives with multiple types of resistance alleles (R1, R2, R3)

This code runs a single instance of Wright-Fisher pop-gen and Beverton-Holt population dynamics simulations of  a diecious population for m gRNAs, where each cleavage target site has 5 alleles: W - wild type;  R (R1) - functional resistant; N (R2) — non-functional resistant; P (R3) - partial resistant allele; D - drive, but where a mosaic of drive on a single chromosome is not possible (e.g. gametes WRD is not possible, but only DDD)

Inputs:
m — number of gRNAs (can in principle do simulations for any number, but computation becomes slower)

psi — initial fraction of males in population 

s — female fitness cost of homozygote with deleterious haplotypes on both chromosomes

hDW — dominance coeff/fitness cost of each wild type allele in haplotype with drive on other chromosome 
hDP — dominance coeff/fitness cost of each partial resistance (R3) allele in haplotype with drive on other chromosome
hN — dominance coeff/fitness cost of each non-functional resistance allele (R2) in haplotype with drive on other chromosome

sigmaR — female fitness cost of each occurrence of R (R1) in a genotype (e.g.  w(WWR/WWR) = (1-σR)^2, w(WWW/WWR) = 1-σR, w(WRR/RRR) = (1-σR)^5 )
sigmaP — female fitness cost of each occurrence of P (R3) in a genotype (e.g.  w(WWP/WWP) = (1-σP)^2, w(WWW/WWP) = 1-σP, w(WPP/PPP) = (1-σP)^5 )

Rm — absolute or intrinsic growth rate of population assuming population is all WWW
K — carrying capacity of population assuming population is all WW...W

xinit — initial frequency vector of alleles in males if scalar then assumes it is =xD — initial frequency of drive in males (DD...D)
yinit — initial frequency vector of alleles in females if scalar then assumes it is =yD — initial frequency of drive in females (DD...D)

epsilonW — efficiency of drive cleavage per target site with W allele
epsilonP — efficiency of drive cleavage per target site with P allele

nuW — non-homologous end-joining (NHEJ) rate per generation per individual per target site with W allele
nuP — non-homologous end-joining (NHEJ) rate per generation per individual per target site with P allele

muSite — mutation rate per cleavage target site per generation per individual per target site (e.g. if a target site contains 10bp and bp mutation rate is 1e-9, set mu=1e-8)

betaXY — fraction of type Y alleles from NHEJ mutations from allele X
(i.e.. betaWR are the fraction of NHEJ mutations from W that give R,
        betaWP are the fraction of NHEJ mutations from W that give P,
        betaPR are the fraction of NHEJ mutations from P that give R,
        betaPW are the fraction of NHEJ mutations from P that give W.


xiR — fraction of functional resistant single nucleotide mutations that produce R (R1) alleles from W
xiP — fraction of functional resistant single nucleotide mutations that produce P (R3) alleles from W

preexisting — preexisting =1 runs simulations for a period 1/sigma as a burn-in period to establish a mutation-selection balance 

T — length of simulations in generations (simulations will end earlier if population eliminated) test_threshold — at what frequency should the sum of all resistance alleles reach to terminate simulations (and signify resistance=1)

Gauss — Gauss = 0 uses multinomial distribution; Gauss = 1 uses Poisson-Gaussian hybrid approximation to the multinomial distribution.
%   If Gauss=0 the maximum population size is ~N=1e9 and is very slow — if Gauss=1, the hybrid approx can be used with arbitrarily large N with little speed penalty
 
Cswitch — Cswitch = 1 if together with Gauss=0 uses GSL C version of multinomial random number generator, which is much quicker than the matlab implementation. Cswitch = 0 uses Matlab's multinomial RNG, unless Gauss=1

plotfig — plotfig=1 plots allele frequency time series and population dynamics; plotfig=0 doesn't
plotmale — plots male allele frequencies or not
vis — vis=0 with plotfig=1 plots figures and makes them invisible and saves them (for remote jobs); their visibility state will be 0 and so when opened in Matlab this needs to be changed, e.g. openfig('Figure.fig','visible')

name — name of directory to create and store figures 



Outputs:
resistance — resistance=0 or 1, whether simulation produced resistance (freq(sum all resistance alleles)>test_threshold); if resistance alleles have not passed threshold by T, then resistance=0
Nend — final population size at the end of simulation
tau — the time when resistance occured (freq(R)>test_threshold) or = NaN if population eliminated
N — total population size over time
M — male population size over time
F — female population size over time
x — matrix whose columns are the frequency of each allele in males at subsequent time points
y — matrix whose columns are the frequency of each allele in females at subsequenttime points
    For both of these the order of alleles is W R N P D — so the 3rd row of x is the time series of the freq(N)    
t — corresponding vector of times (generations)


# MR Power Functions 

# Calculate the proportion of variance explained 
## Shim, Heejung, Daniel I. Chasman, Joshua D. Smith, Samia Mora, Paul M. Ridker, Deborah A. Nickerson, Ronald M. Krauss, and Matthew Stephens.  
## “A Multivariate Genome-Wide Association Analysis of 10 LDL Subfractions, and Their Response to Statin Treatment, in 1868 Caucasians.” 
## PloS One 2015; 10(4): e0120758.
snp.pve <- function(eaf, beta, se, n){
  (2*eaf*(1 - eaf)*beta^2) / (2 * beta * eaf * (1-eaf) + se^2 * 2 * n * eaf * (1-eaf))
}

# Calculating F-Statistic
## Burgess S, Thompson SG, CRP CHD Genetics Collaboration. 
## Avoiding bias from weak instruments in Mendelian randomization studies. 
## Int J Epidemiol. 2011;40: 755–764.
f_stat = function(N, K, R){
  f = ((N-K-1) / K) * (R/(1-R))}

# Calculating Power for Binary Outcomes 
## Konstantin Shakhbazov, Peter M Visscher
## Calculating statistical power in Mendelian randomization studies Marie-Jo A Brion
## Int J Epidemiol 2013 42: 1497-1501
## See http://cnsgenomics.com/shiny/mRnd/
results_binary <- function(Nexposure, Noutcome, nsnps, alpha, R2xz, K, OR) {
  threschi <- qchisq(1 - alpha, 1) # threshold chi(1) scale
  f = f_stat(Nexposure, nsnps, R2xz)
  
  b_MR <- K * ( OR/ (1 + K * (OR - 1)) -1)
  
  v_MR <- (K * (1-K) - b_MR^2) / (Noutcome*R2xz)
  NCP <- b_MR^2 / v_MR
  
  # 2-sided test
  power <- 1 - pchisq(threschi, 1, NCP)
  tibble(Parameter = c("Sample Size (outcome)", 'Proportion of Cases', 'Sample Size (exposure)', 'R2', 'F', "Alpha", 'OR', "NCP", "Power"), 
         Value = c(Noutcome, K, Nexposure, R2xz, f, alpha, OR, NCP, power))
}

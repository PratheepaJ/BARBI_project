// copyright to Kris Sankaran https://github.com/krisrs1128/boot_expers/blob/master/ldaSim/inst/extdata/lda.stan.

data {
  int<lower=1> K; // num topics
  int<lower=1> V; // num species
  int<lower=1> D; // num samples
  int<lower=0> n[D, V]; // species counts for each doc

  // hyperparameters
  vector<lower=0>[K] alpha;
  vector<lower=0>[V] gamma;
}

parameters {
  simplex[K] theta[D]; // topic mixtures for each sample
  simplex[V] beta[D, K]; // word dist for k^th topic in d^th sample
}

model {
  for (d in 1:D) {
    theta[d] ~ dirichlet(alpha);
  }

  for (d in 1:D) {
    for (k in 1:K) {
      beta[d, k] ~ dirichlet(gamma);
    }
  }
 
  for (d in 1:D) {
    vector[V] eta[d];
    
    eta[d] = beta[d,1] * theta[d, 1];

    for (k in 2:K) {
      eta[d] = eta[d] + beta[d,k] * theta[d, k];
    }
    
    n[d] ~ multinomial(eta[d]);
    
  }
}

data {
  int<lower=1> nTime;
  int counts[nTime];
  int vaccine[nTime];
  int drop[nTime];
  int nCounts[nTime];
  int nVaccine[nTime];
  int nDrop[nTime];
}

parameters{
  real meanWeek1;
  vector[nTime-1] changes;
  real<lower=0> changeSigma;
  real vaccineChange;
  real dropChange;
}
transformed parameters{
  vector[nTime] means;
  means[1]=meanWeek1;
  for(ii in 2:nTime) means[ii]=means[ii-1]+changes[ii-1]*changeSigma;
}


model {
  meanWeek1 ~ normal(0,10);
  changes[1] ~ normal(0,1);
  for(ii in 2:(nTime-1)) changes[ii]~normal(changes[ii-1],1);
  counts ~ binomial_logit(nCounts,means);
  vaccine ~ binomial_logit(nVaccine,means+vaccineChange);
  drop ~ binomial_logit(nDrop,means+dropChange);
  changeSigma~gamma(1,2);
  vaccineChange~double_exponential(0,1);
  dropChange~double_exponential(0,1);
}

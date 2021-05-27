


data {
    int<lower=1> nLineage;
    int<lower=1> nTime;
    int counts[nLineage,nTime];
}

parameters{
    vector[nLineage-1] meansWeek1;
    matrix[nLineage,nTime-1] changes;
    real<lower=0> changeSigma;
}


transformed parameters{
  matrix[nLineage,nTime] props;
  matrix[nLineage,nTime] means;
  means[1,1]=0.0;
  means[2:nLineage,1]=meansWeek1;
  for(ii in 2:nTime){
    means[,ii]=means[,ii-1]+changes[,ii-1]*changeSigma;
  }
  for(ii in 1:nTime){
    props[,ii]=exp(means[,ii])/sum(exp(means[,ii]));
  }
}


model {
    meansWeek1 ~ normal(0,10);
    changes[,1] ~ normal(0,1);
    for(ii in 2:(nTime-1)){
      changes[,ii]~normal(changes[,ii-1],1);
    }
    for(ii in 1:nTime)counts[,ii] ~ multinomial(props[,ii]);
    changeSigma~gamma(1,1);
}

data {
    int<lower=1> nLineage;
    int<lower=1> nTime;
    int counts[nLineage,nTime];
    int vaccine[nLineage,nTime];
    int drop[nLineage,nTime];
}

parameters{
    vector[nLineage-1] meansWeek1;
    matrix[nLineage-1,nTime-1] changes;
    real<lower=0> changeSigma;
    vector[nLineage] vaccineChange;
    vector[nLineage] dropChange;
}


transformed parameters{
  matrix[nLineage,nTime] props;
  matrix[nLineage,nTime] propsVaccine;
  matrix[nLineage,nTime] propsDrop;
  matrix[nLineage,nTime] means;
  means[1,1]=0.0;
  means[2:nLineage,1]=meansWeek1;
  for(ii in 2:nTime){
    means[1,ii]=0.0;
    means[2:nLineage,ii]=means[2:nLineage,ii-1]+changes[,ii-1]*changeSigma;
  }
  for(ii in 1:nTime){
    props[,ii]=exp(means[,ii])/sum(exp(means[,ii]));
    propsVaccine[,ii]=exp(means[,ii]+vaccineChange)/sum(exp(means[,ii]+vaccineChange));
    propsDrop[,ii]=exp(means[,ii]+dropChange)/sum(exp(means[,ii]+dropChange));
  }
}


model {
    meansWeek1 ~ normal(0,10);
    changes[,1] ~ normal(0,1);
    for(ii in 2:(nTime-1)){
      changes[,ii]~normal(changes[,ii-1],1);
    }
    for(ii in 1:nTime){
      counts[,ii] ~ multinomial(props[,ii]);
      vaccine[,ii] ~ multinomial(propsVaccine[,ii]);
      drop[,ii] ~ multinomial(propsDrop[,ii]);
    }
    changeSigma~gamma(1,2);
    vaccineChange~double_exponential(0,1);
    dropChange~double_exponential(0,1);
}

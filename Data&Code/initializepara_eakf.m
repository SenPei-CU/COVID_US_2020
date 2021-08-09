function [para,paramax,paramin,betamap,alphamap]=initializepara_eakf(dailyincidence,num_ens,parafit)
%initialize ensemble members
%para: Z,D,mu,theta,alpha1,alpha2,...,alpha3142,beta1,...,beta3142
Zlow=2;Zup=5;%latency period
Dlow=2;Dup=5;%infectious period
mulow=0.2;muup=1;%relative transmissibility
thetalow=0;thetaup=0.2;%movement factor
alphalow=0.05;alphaup=1;%reporting rate
betalow=0.2;betaup=2.5;%transmission rate
num_loc=size(dailyincidence,1);
%define alpha
alphamap=4+(1:num_loc)';
%define beta
betamap=4+num_loc+(1:num_loc)';
load popd%population density
PD=log10(popd);
PD_median=median(PD);
%Z,D,mu,theta,alpha1,alpha2,...,alpha3142,beta1,...,beta3142
paramin=[Zlow;Dlow;mulow;thetalow;ones(num_loc,1)*alphalow;ones(num_loc,1)*betalow];
paramax=[Zup;Dup;muup;thetaup;ones(num_loc,1)*alphaup;ones(num_loc,1)*betaup];
para=zeros(size(paramin,1),num_ens);
%parafit:beta;mu;Z;D;alpha;theta
%Z
para(1,:)=median(parafit(3,:))*ones(1,num_ens);
%D
para(2,:)=median(parafit(4,:))*ones(1,num_ens);
%mu
para(3,:)=median(parafit(2,:))*ones(1,num_ens);
%theta
para(4,:)=median(parafit(6,:))*ones(1,num_ens);
for l=1:num_loc
    %alpha
    if sum(dailyincidence(l,1:24))>=20%at least 20 cases by March 15
        para(alphamap(l),:)=parafit(5,ceil(rand(1,num_ens)*size(parafit,2)))*0.6;
    else
        para(alphamap(l),:)=parafit(5,ceil(rand(1,num_ens)*size(parafit,2)));
    end
    %beta
    if sum(dailyincidence(l,1:24))>=20%at least 20 cases by March 15
        scale=1;
        factor=max(0.1,PD(l)/PD_median*scale);
    else
        scale=0.5;
        factor=max(0.1,min(1.5,PD(l)/PD_median*scale));
    end
    para(betamap(l),:)=parafit(1,ceil(rand(1,num_ens)*size(parafit,2)))*factor;
end
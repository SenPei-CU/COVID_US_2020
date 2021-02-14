function infer()
% Run "mex model_eakf.cpp" to compile the cpp program
%%%%%%%%%%%%%%%%%loading data
load delaypara%load delay parameters
load commutedata%load inter-county commuting data
load population%load county population
load countyfips%load county information
load dailyincidence%load county-level daily reported cases
load parafit1%load initial parameters
load MI_inter%load inter-county visitor numbers
%%%%%%%%%%%%%%%%%%%%
rnds=zeros(1e4,size(delaypara,1));%delay for reporting
for t=1:size(delaypara,1)
    Td=delaypara(t,1)*delaypara(t,2)+2.5;
    a=delaypara(t,1);
    b=Td/a;
    rnds(:,t)=ceil(gamrnd(a,b,1e4,1));
end
%%%%%%%%%%%%%%%%%%%%
num_loc=size(part,1)-1;
num_mp=size(nl,1);
num_times=size(dailyincidence,2);%total length of data
obs_case=zeros(size(dailyincidence));
%smooth the data: 7 day moving average
for l=1:num_loc
    for t=1:num_times
        if (t+3)<=num_times
            obs_case(l,t)=(mean(dailyincidence(l,max(1,t-3):min(t+3,num_times))));
        else
            obs_case(l,t)=(mean(dailyincidence(l,max(1,num_times-6):num_times)));
        end
    end
end
T=num_times;%number of days to keep track of reported cases
startday='02/21/2020';
Tstart=datetime(startday);
%set OEV, observation error variance
OEV_case=zeros(size(dailyincidence));
for l=1:num_loc
    for t=1:num_times
        obs_ave=mean(dailyincidence(l,max(1,t-6):t));
        OEV_case(l,t)=max(1e-4,obs_ave^2/100);
    end
end
%%%%%%%%%%%%%%%%%%%%%%
%adjusting inter-county movement
%MI_inter start from March 1st, first column is fips code for each county
%adjusting mobility starting from March 1, day 10
MI_inter=MI_inter(:,2:end);
MI_inter(:,2:end)=min(MI_inter(:,2:end),1);
%relative daily change of inter-county visitor numbers
MI_inter_relative=[ones(num_loc,9),MI_inter(:,2:end)];
for t=10:size(MI_inter,2)
    MI_inter_relative(:,t+9)=MI_inter(:,t)./MI_inter(:,t-1);
    MI_inter_relative(isnan(MI_inter_relative(:,t+9)),t+9)=0;
    MI_inter_relative(isinf(MI_inter_relative(:,t+9)),t+9)=0;
    MI_inter_relative(:,t+9)=min(MI_inter_relative(:,t+9),1.01);
end
MI_inter_relative(:,1:10)=1;
C=C*ones(1,T);%Daily subpopulation size
Cave=Cave*ones(1,T);
for t=10:T%day 10, March 1
    C(:,t)=C(:,t-1);
    Cave(:,t)=Cave(:,t-1);
    for l=1:num_loc
        for j=part(l)+1:part(l+1)-1
            if t<=size(MI_inter_relative,2)
                C(part(l),t)=C(part(l),t)+((1-MI_inter_relative(nl(j),t))*C(j,t));
                Cave(part(l),t)=Cave(part(l),t)+((1-MI_inter_relative(nl(j),t))*Cave(j,t));
                C(j,t)=(MI_inter_relative(nl(j),t)*C(j,t));
                Cave(j,t)=(MI_inter_relative(nl(j),t)*Cave(j,t));
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
num_ens=100;%number of ensemble members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize model
[S,E,Ir,Iu,Seedc]=initialize(nl,part,C(:,1),num_ens,dailyincidence);
obs_temp=zeros(num_loc,num_ens,T);%records of reported cases
%initialize parameters
[para,paramax,paramin,betamap,alphamaps]=initializepara_eakf(dailyincidence,num_ens,parafit);
para_ori=para;
%%%%%%%%%%%%inflate variables and parameters to get larger spread
lambda=2;
para=mean(para,2)*ones(1,num_ens)+lambda*(para-mean(para,2)*ones(1,num_ens));
para=checkbound_para(para,paramax,paramin);
E=mean(E,2)*ones(1,num_ens)+lambda*(E-mean(E,2)*ones(1,num_ens));
Ir=mean(Ir,2)*ones(1,num_ens)+lambda*(Ir-mean(Ir,2)*ones(1,num_ens));
Iu=mean(Iu,2)*ones(1,num_ens)+lambda*(Iu-mean(Iu,2)*ones(1,num_ens));
[S,E,Ir,Iu]=checkbound(S,E,Ir,Iu,C(:,1));
%fix Z, D, mu and theta
para(1,:)=parafit(3,1:num_ens);%Z, latency period
para(2,:)=parafit(4,1:num_ens);%D, infectious period
para(3,:)=parafit(2,1:num_ens);%mu, relative transmissibility
para(4,:)=parafit(6,1:num_ens);%theta, mobility factor
parastd=std(para,0,2);%get ensemble spread of parameters
%%%%%%%%%%%%%%%%%%%%%
dailyIr_post_rec=zeros(num_loc,num_ens,T);%posterior reported infection
dailyIu_post_rec=zeros(num_loc,num_ens,T);%posterior unreported infection
%initialize poseteriors (percentage of population)
S_post=zeros(num_loc,num_times,num_ens);
E_post=zeros(num_loc,num_times,num_ens);
Ir_post=zeros(num_loc,num_times,num_ens);
Iu_post=zeros(num_loc,num_times,num_ens);
%initialize cumulative reported and unreported infections
cumu_dailyIr_post=zeros(num_mp,num_ens);
cumu_dailyIu_post=zeros(num_mp,num_ens);
%%%%%%%%%%%%%%%%%%%%%%%
lambda=1.2;%inflation in EAKF to avoid divergence
num_para=size(para,1);%number of parameters
para_post=zeros(num_para,num_ens,num_times);%posterior parameters
for t=1:num_times%start from Feb 21, stop on Dec 31
    t
    tic
    if t<=41%April 1st
        locs=find(sum(dailyincidence(:,1:24),2)<20);
    else
        locs=(1:num_loc)';
    end
    %increase lower bound of reporting rate from 5% to 15% on 06/01
    paramin(4+(locs))=min(0.05+0.1*(t-1)/100,0.15);
    %fix Z, D and theta
    para(1,:)=parafit(3,1:num_ens);%Z
    para(2,:)=parafit(4,1:num_ens);%D
    para(3,:)=parafit(2,1:num_ens);%mu
    para(4,:)=parafit(6,1:num_ens);%theta
    %seeding
    if t<=size(Seedc,2)
        [S,E,Ir,Iu]=seeding(S,E,Ir,Iu,nl,part,C(:,t),Seedc,t);
    end
    %re-initialize beta
    for l=1:num_loc
        if Seedc(l,t)>0
            para(betamap(l),:)=para_ori(betamap(l),:);
        end
    end
    para=checkbound_para(para,paramax,paramin);
    %%%%%%%%%%%%%%%%%%
    S_temp=S; E_temp=E; Ir_temp=Ir; Iu_temp=Iu;
    %integrate forward one step
    dailyIr_prior=zeros(num_mp,num_ens);
    dailyIu_prior=zeros(num_mp,num_ens);
    for k=1:num_ens%run for each ensemble member
        [S(:,k),E(:,k),Ir(:,k),Iu(:,k)]=adjustmobility(S(:,k),E(:,k),Ir(:,k),Iu(:,k),nl,part,MI_inter_relative,t);
        [S(:,k),E(:,k),Ir(:,k),Iu(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t),Cave(:,t),S(:,k),E(:,k),Ir(:,k),Iu(:,k),para(:,k),betamap,alphamaps);
        dailyIr_prior(:,k)=dailyIr_temp;
        dailyIu_prior(:,k)=dailyIu_temp;
    end
    %%%%%%%%%%%%%%%%%%%%%%
    %integrate forward for 9 days, prepare for observation
    Tproj=9;
    obs_temp1=obs_temp;%
    for t1=t:t+Tproj-1
        for k=1:num_ens%run for each ensemble member
            [S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k)]=adjustmobility(S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k),nl,part,MI_inter_relative,t1);
            [S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t1),Cave(:,t1),S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k),para(:,k),betamap,alphamaps);
            %reporting delay
            monthnum=max(month(Tstart+t1-1)-3,1);
            monthnum=min(monthnum, size(rnds,2));
            for l=1:num_loc
                for j=part(l):part(l+1)-1
                    inci=dailyIr_temp(j);
                    if inci>0
                        rnd=datasample(rnds(:,monthnum),inci);
                        for h=1:length(rnd)
                            if (t1+rnd(h)<=T)
                                obs_temp1(l,k,t1+rnd(h))=obs_temp1(l,k,t1+rnd(h))+1;
                            end
                        end
                    end
                end
            end

        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t1=min(t+Tproj,num_times)
        obs_ens=obs_temp1(:,:,t1);%observation at t1, prior
        %loop through local observations
        for l=1:num_loc
            %%%%%%%%%%%%%%%%%%%case
            if (sum(obs_case(l,1:t1))>0)
                %Get the variance of the ensemble
                obs_var = OEV_case(l,t1);
                prior_var = var(obs_ens(l,:));
                post_var = prior_var*obs_var/(prior_var+obs_var);
                if prior_var==0%if degenerate
                    post_var=1e-3;
                    prior_var=1e-3;
                end
                prior_mean = mean(obs_ens(l,:));
                post_mean = post_var*(prior_mean/prior_var + obs_case(l,t1)/obs_var);
                %%%% Compute alpha and adjust distribution to conform to posterior moments
                alpha = (obs_var/(obs_var+prior_var)).^0.5;
                dy = post_mean + alpha*(obs_ens(l,:)-prior_mean)-obs_ens(l,:);
                %Loop over each state variable (connected to location l)
                %adjust related metapopulation
                neighbors=part(l):part(l+1)-1;%metapopulation live in l
                for h=1:length(neighbors)
                    j=neighbors(h);
                    %E
                    temp=E(j,:);
                    A=cov(temp,obs_ens(l,:));
                    rr=A(2,1)/prior_var;
                    dx=rr*dy;
                    E(j,:)=E(j,:)+dx;
                    %Ir
                    temp=Ir(j,:);
                    A=cov(temp,obs_ens(l,:));
                    rr=A(2,1)/prior_var;
                    dx=rr*dy;
                    Ir(j,:)=Ir(j,:)+dx;
                    %Iu
                    temp=Iu(j,:);
                    A=cov(temp,obs_ens(l,:));
                    rr=A(2,1)/prior_var;
                    dx=rr*dy;
                    Iu(j,:)=Iu(j,:)+dx;
                    %dailyIr
                    temp=dailyIr_prior(j,:);
                    A=cov(temp,obs_ens(l,:));
                    rr=A(2,1)/prior_var;
                    dx=rr*dy;
                    dailyIr_prior(j,:)=round(max(dailyIr_prior(j,:)+dx,0));
                    %dailyIu
                    temp=dailyIu_prior(j,:);
                    A=cov(temp,obs_ens(l,:));
                    rr=A(2,1)/prior_var;
                    dx=rr*dy;
                    dailyIu_prior(j,:)=round(max(dailyIu_prior(j,:)+dx,0));
                end
                %adjust alpha
                temp=para(alphamaps(l),:);
                A=cov(temp,obs_ens(l,:));
                rr=A(2,1)/prior_var;
                dx=rr*dy;
                para(alphamaps(l),:)=para(alphamaps(l),:)+dx;
                %inflation
                if std(para(alphamaps(l),:))<parastd(alphamaps(l))
                    para(alphamaps(l),:)=mean(para(alphamaps(l),:),2)*ones(1,num_ens)+lambda*(para(alphamaps(l),:)-mean(para(alphamaps(l),:),2)*ones(1,num_ens));
                end
                %adjust beta
                temp=para(betamap(l),:);
                A=cov(temp,obs_ens(l,:));
                rr=A(2,1)/prior_var;
                dx=rr*dy;
                para(betamap(l),:)=para(betamap(l),:)+dx;
                %inflation
                if std(para(betamap(l),:))<parastd(betamap(l))
                    para(betamap(l),:)=mean(para(betamap(l),:),2)*ones(1,num_ens)+lambda*(para(betamap(l),:)-mean(para(betamap(l),:),2)*ones(1,num_ens));
                end
            else
                %adjust related metapopulation
                neighbors=part(l):part(l+1)-1;%metapopulation live in l
                for h=1:length(neighbors)
                    j=neighbors(h);
                    %S
                    S(j,:)=C(j,t);
                    %E
                    E(j,:)=0;
                    %Ir
                    Ir(j,:)=0;
                    %Iu
                    Iu(j,:)=0;
                    %dailyIr
                    dailyIr_prior(j,:)=0;
                    %dailyIu
                    dailyIu_prior(j,:)=0;
                end
            end
        end
        para=checkbound_para(para,paramax,paramin);
    end
    %update posterior Ir and Iu
    dailyIr_post=dailyIr_prior;
    dailyIu_post=dailyIu_prior;
    
    cumu_dailyIr_post=cumu_dailyIr_post+dailyIr_post;
    cumu_dailyIu_post=cumu_dailyIu_post+dailyIu_post;
    
    monthnum=max(month(Tstart+t-1)-3,1);
    monthnum=min(monthnum, size(rnds,2));
    for k=1:num_ens
        %update obs_temp
        for l=1:num_loc
            for j=part(l):part(l+1)-1
                inci=dailyIr_post(j,k);
                if inci>0
                    rnd=datasample(rnds(:,monthnum),inci);
                    for h=1:length(rnd)
                        if (t+rnd(h)<=T)
                            obs_temp(l,k,t+rnd(h))=obs_temp(l,k,t+rnd(h))+1;
                        end
                    end
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%update S
    S=C(:,t)*ones(1,num_ens)-E-cumu_dailyIr_post-cumu_dailyIu_post;
    %%%%%%%%%%%%%%%%
    [S,E,Ir,Iu]=checkbound(S,E,Ir,Iu,C(:,t));
    %%%%%%%%save posterior statevariables
    for i=1:num_loc
        for j=1:num_ens
            S_post(i,t,j)=sum(S(part(i):part(i+1)-1,j))./population(i);
            E_post(i,t,j)=sum(E(part(i):part(i+1)-1,j))./population(i);
            Ir_post(i,t,j)=sum(Ir(part(i):part(i+1)-1,j))./population(i);
            Iu_post(i,t,j)= sum(Iu(part(i):part(i+1)-1,j))./population(i);
            dailyIr_post_rec(i,j,t)=sum(dailyIr_post(part(i):part(i+1)-1,j));
            dailyIu_post_rec(i,j,t)=sum(dailyIu_post(part(i):part(i+1)-1,j));
        end
    end
    para_post(:,:,t)=para;
    toc
end
%save outputs
save inference.mat para_post S_post dailyIr_post_rec dailyIu_post_rec obs_temp

function para = checkbound_para(para,paramax,paramin)
for i=1:size(para,1)
    para(i,para(i,:)<paramin(i))=paramin(i)*(1+0.2*rand(sum(para(i,:)<paramin(i)),1));
    para(i,para(i,:)>paramax(i))=paramax(i)*(1-0.2*rand(sum(para(i,:)>paramax(i)),1));
end

function [S,E,Ir,Iu]=checkbound(S,E,Ir,Iu,C)
for k=1:size(S,2)
    S(S(:,k)<0,k)=0; E(E(:,k)<0,k)=0; Ir(Ir(:,k)<0,k)=0; Iu(Iu(:,k)<0,k)=0;
end


function [S,E,Ir,Iu]=adjustmobility(S,E,Ir,Iu,nl,part,MI_inter_relative,t)
num_loc=size(MI_inter_relative,1);
for l=1:num_loc
    for j=part(l)+1:part(l+1)-1
        if t<=size(MI_inter_relative,2)
            S(part(l))=S(part(l))+((1-MI_inter_relative(nl(j),t))*S(j));
            S(j)=(MI_inter_relative(nl(j),t)*S(j));
            E(part(l))=E(part(l))+((1-MI_inter_relative(nl(j),t))*E(j));
            E(j)=(MI_inter_relative(nl(j),t)*E(j));
            Ir(part(l))=Ir(part(l))+((1-MI_inter_relative(nl(j),t))*Ir(j));
            Ir(j)=(MI_inter_relative(nl(j),t)*Ir(j));
            Iu(part(l))=Iu(part(l))+((1-MI_inter_relative(nl(j),t))*Iu(j));
            Iu(j)=(MI_inter_relative(nl(j),t)*Iu(j));
        end
    end
end
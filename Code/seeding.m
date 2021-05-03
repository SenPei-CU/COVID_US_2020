function [S,E,Ir,Iu]=seeding(S,E,Ir,Iu,nl,part,C,Seedc,t)
%set initial infection for all subpopulations
num_loc=size(Seedc,1);
num_ens=size(S,2);
for l=1:num_loc
    seedloc=l;
    seedE=round(rand(1,num_ens)*12*Seedc(l,t));
    seedIu=round(rand(1,num_ens)*10*Seedc(l,t));
    seedIr=round(rand(1,num_ens)*2*Seedc(l,t));
    if Seedc(l,t)>0
        pop=sum(C(part(seedloc):part(seedloc+1)-1));
        for i=part(seedloc):part(seedloc+1)-1
            E(i,:)=round(seedE*C(i)/pop);
            Iu(i,:)=round(seedIu*C(i)/pop);
            Ir(i,:)=round(seedIr*C(i)/pop);
            S(i,:)=S(i,:)-E(i,:)-Ir(i,:)-Iu(i,:);
        end
    end
end
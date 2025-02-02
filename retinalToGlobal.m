

newRes=0.5;
i=1;
j=1;
for alpha=0:newRes:360
    for beta=0:newRes:180
        
        if  beta==0
            beta=180;
        end
        [aziR(j,i),elevR(j,i),~]=rotateNormalToExtraPersonalSpace4('Right',alpha,beta);
        [aziL(j,i),elevL(j,i),~]=rotateNormalToExtraPersonalSpace4('Left',alpha,beta);
        alphalookup(j,i)=alpha;
        betalookup(j,i)=beta;
        j=j+1;
    end
    j=1;
    i=i+1;
end

save('retinal to global lookup table.mat','aziR','elevR','aziL','elevL','alphalookup','betalookup');
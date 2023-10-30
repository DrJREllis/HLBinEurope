function [forceInfect,forceVector,longDistDisp,estab,...
    crypticRate,symptomRate,forceVectorOth,forceInfectOth,...
    lambdaInfectionSum,lambdaVectorSum,lambdaInfectionSumOth,lambdaVectorSumOth,VecTime,InfTime] = ...
    updateRates(N,E,C,I,R,newR,fullInfestation,forceInfect,forceVector,comInd,M,...
    forceVectorOth,forceInfectOth,lambdaInfectionSum,lambdaVectorSum,lambdaInfectionSumOth,lambdaVectorSumOth,...
    longDistDisp,estab,longDistDispOth,extInfCSum,extVCSum,extInfCSumOth,extVCSumOth,estabOth,...
    crypticRate,symptomRate,latentRateOth,symptomRateOth,beta,localDispRate,longDispRate,estRate,lRate,sRate,...
    citArray,adjN,climArray,carrCap,othCC,Noth,Loth,Coth,Ioth,Roth,othArray,VecTime,InfTime,Fi,Kd,Ki,cell,event,t)
% Function to generate new rates for each event after previous event

if isnan(cell)
    return
end


if sum(event==[1,2,3,5,6,8])==1 % VECTOR 0->A EVENT
    
    estab(cell) = estRate;

    forceVector(cell) = 0;
    lambdaVectorSum = sum(forceVector) + sum(longDistDisp) + extVCSum(end) + sum(estab);

    if isnan(VecTime(cell))
      VecTime(cell)=t;
    end

end

if sum(event==9)==1 % VECTOR A->N EVENT

    estab(cell) = 0;

    forceVector(Fi{cell}) = (forceVector(Fi{cell}) + adjN(cell)*localDispRate *Kd{cell}.* climArray(Fi{cell}) .* max(M(Fi{cell}),1-comInd)) .*(carrCap(Fi{cell})>N(Fi{cell}));
    forceVectorOth(Fi{cell}) = (forceVectorOth(Fi{cell}) + adjN(cell)*localDispRate *Kd{cell}.* climArray(Fi{cell}) .* max(M(Fi{cell}),comInd)) .*(othCC(Fi{cell})>Noth(Fi{cell}));
    


    longDistDisp(cell) = N(cell) * adjN(cell) * longDispRate;



    
    if carrCap(cell)<=N(cell)
        forceVector(cell) = 0;
    end            
    lambdaVectorSum = sum(forceVector) + sum(longDistDisp) + extVCSum(end) + sum(estab);
    lambdaVectorSumOth = sum(forceVectorOth) + sum(longDistDispOth) + extVCSumOth(end) + sum(estabOth);

    if I(cell)+C(cell)>0
        S=citArray-E-C-I-R;
        Soth=othArray-Loth-Coth-Ioth-Roth;
        forceInfect(Fi{cell}) = forceInfect(Fi{cell}) + beta * Ki{cell}* (C(cell)+I(cell))*adjN(cell)...
            .* max(0,S(Fi{cell})./(citArray(Fi{cell}).^2));
        forceInfectOth(Fi{cell}) = forceInfectOth(Fi{cell}) + beta * Ki{cell}* (C(cell)+I(cell))*adjN(cell)...
            .* max(0,Soth(Fi{cell})./(othArray(Fi{cell}).^2));
    
        lambdaInfectionSum = sum(forceInfect) + extInfCSum(end) + sum(crypticRate) + sum(symptomRate);
        lambdaInfectionSumOth = sum(forceInfectOth) + extInfCSumOth(end) + sum(latentRateOth) + sum(symptomRateOth);
    end

end


if sum(event==[3,4,6,7,10,11])==1 % CITRUS S->L EVENT
    S=citArray-E-C-I-R;
    forceInfect(cell) = forceInfect(cell) * (S(cell))/ (S(cell)+1) ;
    crypticRate(cell) = E(cell)*lRate;

    if isnan(InfTime(cell))
      InfTime(cell)=t;
    end

    lambdaInfectionSum = sum(forceInfect) + extInfCSum(end) + sum(crypticRate) + sum(symptomRate);
end

if sum(event==12)==1 %CITRUS E->C EVENT
    S=citArray-E-C-I-R;
    Soth=othArray-Loth-Coth-Ioth-Roth;
    forceInfect(Fi{cell}) = forceInfect(Fi{cell}) + beta * Ki{cell}* N(cell)*adjN(cell)...
        .* max(0,S(Fi{cell})./(citArray(Fi{cell}).^2));
    forceInfectOth(Fi{cell}) = forceInfectOth(Fi{cell}) + beta * Kd{cell}* N(cell)*adjN(cell)...
        .* max(0,Soth(Fi{cell})./(othArray(Fi{cell}).^2));

    crypticRate(cell) = E(cell)*lRate;
    symptomRate(cell)= symptomRate(cell)+sRate;

    lambdaInfectionSum = sum(forceInfect) + extInfCSum(end) + sum(crypticRate) + sum(symptomRate);
    lambdaInfectionSumOth = sum(forceInfectOth) + extInfCSumOth(end) + sum(latentRateOth) + sum(symptomRateOth);

    if fullInfestation==1
        longDistDisp(cell) = N(cell) * adjN(cell) * longDispRate * ((I(cell)+C(cell))/citArray(cell));
        lambdaVectorSum = sum(longDistDisp);
    end
end


if sum(event==13)==1 %CITRUS C->I EVENT
    symptomRate(cell)=symptomRate(cell)-sRate;
    lambdaInfectionSum = sum(forceInfect) + extInfCSum(end) + sum(crypticRate) + sum(symptomRate);
end

if sum(event==14)==1 %CITRUS I->R EVENT
    S=citArray-E-C-I-R;
    Soth=othArray-Loth-Coth-Ioth-Roth;


    if fullInfestation==1
        longDistDisp(cell) = N(cell) * (adjN(cell)* newR/(S(cell)+newR)) * longDispRate * ((I(cell)+C(cell))/citArray(cell));
        lambdaVectorSum = sum(longDistDisp);

    else
        forceVector(Fi{cell}) = (forceVector(Fi{cell}) - (adjN(cell)* newR/(S(cell)+newR))...
            *localDispRate *Kd{cell}.* climArray(Fi{cell}) .* max(M(Fi{cell}),1-comInd)) .*(carrCap(Fi{cell})>N(Fi{cell}));
        forceVectorOth(Fi{cell}) = (forceVectorOth(Fi{cell}) - (adjN(cell)* newR/(S(cell)+newR))...
            *localDispRate *Kd{cell}.* climArray(Fi{cell}) .* max(M(Fi{cell}),comInd)) .*(othCC(Fi{cell})>Noth(Fi{cell}));

        lambdaVectorSum = sum(forceVector) + sum(longDistDisp) + extVCSum(end) + sum(estab);
        lambdaVectorSumOth = sum(forceVectorOth) + sum(longDistDispOth) + extVCSumOth(end) + sum(estabOth);
    end

for r=1:newR
    forceInfect(Fi{cell}) = max(0,forceInfect(Fi{cell}) - beta * Ki{cell}* N(cell)*adjN(cell)...
        .* max(0,S(Fi{cell})./(citArray(Fi{cell}).^2)));
    forceInfectOth(Fi{cell}) = forceInfectOth(Fi{cell}) - beta * Kd{cell}* N(cell)*adjN(cell)...
        .* max(0,Soth(Fi{cell})./(othArray(Fi{cell}).^2));
end
    forceInfect(Fi{cell}) = forceInfect(Fi{cell}) - beta * Ki{cell}* (C(cell)+I(cell))*(adjN(cell)* newR/(S(cell)+newR))...
            .* max(0,S(Fi{cell})./(citArray(Fi{cell}).^2));
        forceInfectOth(Fi{cell}) = forceInfectOth(Fi{cell}) - beta * Ki{cell}* (C(cell)+I(cell))*(adjN(cell)* newR/(S(cell)+newR))...
            .* max(0,Soth(Fi{cell})./(othArray(Fi{cell}).^2));

    adjN(cell) = adjN(cell)*(1 - newR/(S(cell)+newR));
    

    lambdaInfectionSum = sum(forceInfect) + extInfCSum(end) + sum(crypticRate) + sum(symptomRate);
    lambdaInfectionSumOth = sum(forceInfectOth) + extInfCSumOth(end) + sum(latentRateOth) + sum(symptomRateOth);

end
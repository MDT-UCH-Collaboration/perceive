function [] = nonAvgTrialCompare(leng, js, trialA, trialB, trialC)

%just looks at raw pspectrum data, no processing 
%by comparing trials we can see if some were good or bad 
%needs four arguments and should be separated by 6 
%ex: trials 1, 6, 13, 19 should be the same contact 
trials = [trialA, trialB, trialC] ;
for i = 1:leng 
    if ismember(i,trials)
        [p,f] = pspectrum(js.LfpMontageTimeDomain(i).TimeDomainData);
        semilogx(f,p);
        semilogx(f,db(p,"power"));
        hold on; 
    end
end 
end 
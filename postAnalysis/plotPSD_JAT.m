





% data1 = data.treatment1.timeDomain;

for i = 1:5
    field = append("treatment", num2str(i));
    [outDAT] = proc1(data.(field).timeDomain); 
    %collect time domain data from given field in Data strcut I made 

    hold on
    plot(outDAT.freq,outDAT.pwr)
end
legend()



function [outDAT] = proc1(dataIN)

hp = highpass(dataIN,0.8,250);
lp = lowpass(hp,59,250);
%dont recall how these numbers were created 
[x,y] = pspectrum(lp,250,'FrequencyLimits',[1 59]);

xl = 10*log10(x);
%plots bode-ish of power from pspectrum
x2 = smoothdata(xl,'gaussian',150);
%smooth the bode 

outDAT.freq = y;
outDAT.pwr = x2;

end
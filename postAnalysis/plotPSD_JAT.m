





% data1 = data.treatment1.timeDomain;

for i = 1:5
    field = append("treatment", num2str(i));
    [outDAT] = proc1(data.(field).timeDomain);

    hold on
    plot(outDAT.freq,outDAT.pwr)
end
legend()









function [outDAT] = proc1(dataIN)

hp = highpass(dataIN,0.8,250);
lp = lowpass(hp,59,250);
[x,y] = pspectrum(lp,250,'FrequencyLimits',[1 59]);

xl = 10*log10(x);

x2 = smoothdata(xl,'gaussian',150);

outDAT.freq = y;
outDAT.pwr = x2;

end
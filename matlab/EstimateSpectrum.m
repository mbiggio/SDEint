clear variables
close all
clc

D = linspace(0,5,10);

for k = 1:numel(D)
    s = system(['./HR_SDE ' num2str(D(k))]);
    disp('Programma Eseguito!')
    data = load('output.dat');
    t = data(:,1);
    v = data(:,2);
    dt = t(2)-t(1);
    fs = 1/dt;
    Hs = spectrum.periodogram;
    P = psd(Hs,v,'Fs',fs,'NFFT',2^nextpow2(numel(t)-1));
    
    figure
    plot(P.Frequencies,P.Data);
    set(gca,'XScale','log','YScale','log');
end

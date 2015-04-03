clear all
close all
clc

data = load('../examples/output_HH_Orio.dat');
isi = diff(data);
[h,x] = hist(isi,4000);
dt = x(2)-x(1);
h = h/sum(h)/dt;
plot(x,h)
set(gca,'yscale','log','xlim',[0 110])
xlabel('ISI length, ms')
ylabel('Probability Density')

print('-dpng', 'ISI_histogram.png')
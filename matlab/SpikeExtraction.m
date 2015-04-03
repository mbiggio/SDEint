
clear variables
close all
clc

diam = 1.0:1.0:6.0;
cv = zeros(size(diam));
rate = zeros(size(diam));
T = 1e4;
transient = T/10;

%%
h = waitbar(0,'Please wait...');
for k = 1:numel(diam)
    waitbar(k / numel(diam));
    data = load(['data_MT/output_HH_diam_' num2str(diam(k),'%2.1f') '.dat']);
    spike_times = data(:,1)/1e3;
    ndx = find(spike_times > transient);
    spike_times = spike_times(ndx);
    
    isi = diff(spike_times);
    cv(k) = std(isi)/mean(isi);
    rate(k) = numel(spike_times)/((T-transient));
end
close(h) 

%%


figure
plot(diam,cv)
set(gca,'YLim',[0.6,1.2])
print('-depsc','CV_2_MT.eps')

figure
plot(diam,rate)
set(gca,'YLim',[0,20])
print('-depsc','rate_2_MT.eps')
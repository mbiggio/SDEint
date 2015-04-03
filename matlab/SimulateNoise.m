clear variables
clear all
close all

data = load('data/output_HH_diam_10.0.dat');
plot(data(:,1),data(:,2));
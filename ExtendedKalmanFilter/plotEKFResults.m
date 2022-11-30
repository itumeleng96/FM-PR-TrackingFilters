clc; clear all; close all;
addpath('../FERS/','../CFAR/','../MeanShiftCluster/');

system("fers ../FERS/Simulation_60_direct.fersxml");
system("fers ../FERS/Simulation_60_echo_2.fersxml");


% h5 Import from FERS simulation
[Ino Qno scale_no] = loadfersHDF5('direct.h5');
[Imov Qmov scale_mov] = loadfersHDF5('echo.h5');


I_Qmov = Imov + j*Qmov;
I_Qmov = I_Qmov.*scale_mov;
I_Qno = Ino + j*Qno;
I_Qno = I_Qno.*scale_no;

I_Qmov=I_Qmov-I_Qno;

fs = 200000;
dopp_bins = 200;
delay = 233e-6;

s1 = I_Qmov;
s2 = I_Qno;

initial=1;
current=fs;                                 %based on samples in transmitted signal
simulation_time = size(I_Qmov,1)/fs ;       %Simulation time: number of data points/sampling frequency


ard = [];

figure('Name','2D image');
%ARD
f=figure(1);
f.Position = [4000 10 1000 800]; 
movegui(f,'northwest');

%CFAR
f2=figure(2);
f2.Position = [4000 10 1000 800]; 
movegui(f2,'northeast');


for i = 1:simulation_time
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    
    %Plot Range-Doppler Map
    [y,ard_] = ardPlot(s1,s2,fs,dopp_bins,delay,i,ard,f);

    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters] = ca_cfarPlot(10*log10(y.'),0.35,fs,dopp_bins,delay,i,f2);                    
     
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids] = meanShiftPlot(targetClusters,10);
    disp(clusterCentroids);
    %Plot tracks from Tracker
    
    ard = ard_;
    %Counting Variables
    initial = current+1;
    current = current + fs;
end
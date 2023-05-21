clc; clear all; close all;
addpath('../FERS/','../CFAR/','../KmeansCentroids','../DPI_Suppression/');


system("fers ../FERS/scenario_4_singleFile.fersxml");


% h5Import from FERS simulation
[Ino Qno scale_no] = loadfersHDF5('direct.h5');
[Imov Qmov scale_mov] = loadfersHDF5('echo.h5');


I_Qmov = Imov + j*Qmov;
I_Qmov = I_Qmov.*scale_mov;
I_Qno = Ino + j*Qno;
I_Qno = I_Qno.*scale_no;

%I_Qmov=I_Qmov-I_Qno;


fs = 200000;
dopp_bins = 200;
delay = 233e-6;
c=299792458;
range_delay = delay*c/2;

proc = struct('cancellationMaxRange_m', range_delay, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m', 12540, ...
              'nSegments', 1, ...
              'nIterations', 20, ...
              'Fs', fs, ...
              'alpha', 0, ...
              'initialAlpha', 0);





s1 = I_Qmov;   %Surv
s2 = I_Qno;    %Ref

initial=1;
current=fs;  %based on samples in transmitted signal
simulation_time = size(I_Qmov,1)/fs ; %Simulation time: number of data points/sampling frequency


figure('Name','2D image');

for i = 1:simulation_time

    s1 = I_Qmov(initial:current); %surv
    s2 = I_Qno(initial:current);  %ref

    s1 = procECA(s2,s1,proc);

    ard_plot(s1,s2,fs,dopp_bins,delay);

    initial = current+1;
    current = current + fs;

    disp(["Time(s):",i]);
    %disp(["Current:",current]);
end
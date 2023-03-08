clc; clear all; close all;
addpath('../FERS/','../CFAR/','../KmeansCentroids','../SuppressionFunction/');

%system("fers ../FERS/Simulation_60_direct.fersxml");
%system("fers ../FERS/Simulation_60_echo_2.fersxml");
%system("fers ../FERS/Simulation_60_Bistatic.fersxml");
system("fers ../FERS/singleFile.fersxml");


% h5Import from FERS simulation
[Ino Qno scale_no] = loadfersHDF5('direct.h5');
[Imov Qmov scale_mov] = loadfersHDF5('echo.h5');


I_Qmov = Imov + j*Qmov;
I_Qmov = I_Qmov.*scale_mov;
I_Qno = Ino + j*Qno;
I_Qno = I_Qno.*scale_no;

% run_ard
fs = 200000;
dopp_bins = 200;
delay = 133e-6;

%I_Qmov=I_Qmov-I_Qno;
%Use the cancellation Function
%Parameters for the cancellation function
proc = struct('cancellationMaxRange_m', 13850, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m', 13734, ...
              'nSegments', 16, ...
              'nIterations', 30, ...
              'Fs', fs, ...
              'alpha', 0, ...
              'initialAlpha', 0);


I_Qmov = procCGLS(I_Qno, I_Qmov, proc);


s1 = I_Qmov;
s2 = I_Qno;

initial=1;
current=fs;  %based on samples in transmitted signal
simulation_time = size(I_Qmov,1)/fs ; %Simulation time: number of data points/sampling frequency
disp(["SimTime:",simulation_time])
%disp(["I_Qmov",size(I_Qmov,1)]);
%disp(["fs",fs]);
%disp(["s1,s2:",size(s1),size(s2)]);

figure('Name','2D image');

for i = 1:simulation_time
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    ard_plot(s1,s2,fs,dopp_bins,delay);

    initial = current+1;
    current = current + fs;

    disp(["Time(s):",i]);
    %disp(["Current:",current]);
end
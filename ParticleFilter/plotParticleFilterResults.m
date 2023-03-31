clc; clear all; close all;

addpath('../FERS/','../CFAR/','../MeanShiftCluster/','../multiTargetTracking/','../DPI_Suppression');

system("fers ../FERS/scenario_1_ref.fersxml");
system("fers ../FERS/scenario_1_surv.fersxml");

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
rdm =[];

figure('Name','2D image');
%ARD
f=figure(1);
f.Position = [4000 10 1000 800]; 
movegui(f,'northwest');

%CFAR
f2=figure(2);
f2.Position = [4000 10 1000 800]; 
movegui(f2,'northeast');

%Multi-Target Tracking 
f3=figure(3);
f3.Position = [4000 10 1000 800]; 
movegui(f3,'southeast');

%Create MTT object
confirmationThreshold=2;
deletionThreshold=4;
gatingThreshold=30;
filterType=3;           %Particle Filter 
multiTargetTracker = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType);

for i = 1:simulation_time
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    
    %Plot Range-Doppler Map
    [y,ard_] = ardPlot(s1,s2,fs,dopp_bins,delay,i,ard,f);

    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters,RDM,rdm_] = ca_cfarPlot(10*log10(y.'),0.1,fs,dopp_bins,delay,i,f2,rdm);                    
    
    
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids] = meanShiftPlot(targetClusters,10,fs,dopp_bins,delay);
    
    %Plot tracks from Tracker - Call Multi-target Tracker
    multiTargetTracker = multiTargetTracker.createNewTracks(clusterCentroids);
    multiTargetTracker = multiTargetTracker.maintainTracks();
    multiTargetTracker = multiTargetTracker.predictionStage();
    multiTargetTracker.plotMultiTargetTracking(fs,dopp_bins,delay,i,f3,RDM)
    multiTargetTracker = multiTargetTracker.updateStage(clusterCentroids);

    ard = ard_;
    rdm= rdm_;
    %Counting Variables
    initial = current+1;
    current = current + fs;
end

%CFAR
f4=figure(4);
f4.Position = [4000 10 1000 800]; 
movegui(f4,'northeast');

%Multi-Target Tracking 
f5=figure(5);
f5.Position = [4000 10 1000 800]; 
movegui(f5,'southeast');

multiTargetTracker.plotRMSE(f4,f5,true,true,simulation_time);


clc; clear all; close all;
addpath('FERS/','CFAR/','KmeansCentroids');

% h5Import from FERS simulation
[Ino Qno scale_no] = loadfersHDF5('direct.h5');
[Imov Qmov scale_mov] = loadfersHDF5('echo.h5');
I_Qmov = Imov + j*Qmov;
I_Qmov = I_Qmov.*scale_mov;
I_Qno = Ino + j*Qno;
I_Qno = I_Qno.*scale_no;
I_Qmov=I_Qmov-I_Qno;

% run_ard
fs = 500000;
dopp_bins = 200;
delay = 233e-6;

s1 = I_Qmov;
s2 = I_Qno;

initial=1;
current=500000;  %based on samples in transmitted signal

%Initialise Particles for Particle Filter
N = 5000; %Number of Particles
initialCentroid = [17,388];
particles = createGaussianParticles(initialCentroid,[10,10],N);
%particles = createUniformParticles([,116],[0,401],N);

weights = ones(N,1)/N;

X_estimated =[];  %Arrary to store kalman estimated values
Centroids = [];
ard = [];
cfar = [];

for i = 1:10
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    [particles_,weights_,X_estimated_,Centroids_,ard_,cfar_]=ardPlotParticle(s1,s2,fs,dopp_bins,delay,particles,weights,i,X_estimated,Centroids,ard,cfar);
    particles = particles_;
    weights = weights_;
    X_estimated = X_estimated_;
    Centroids = Centroids_;
    ard = ard_;
    cfar = cfar_;
    
    initial = current+1;
    current = current + fs;
end
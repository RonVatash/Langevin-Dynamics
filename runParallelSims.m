numOfSims = 1;
ExperimentTime = 60; %sec
Dt = 1e-4;
N % The number of steps to run the simulation for
sampleRate = 30; % fps
R = 1.5e-6; % radius in meter
T  = 300;פפפפפ
eta = 0.001;
lx = 11e-6;
ly = 11e-6;
numOfParticles = 6;
wallShrink = 0;
baseFoldername = 'Parallel Harmonic ';
displayLive = true;

parfor i=1:numOfSims
    currFoldername = strcat(baseFoldername,num2str(i));
    infoChamber(N, Dt, sampleRate, R, T, eta, lx, ly, numOfParticles, wallShrink, currFoldername, displayLive);
end

%% Combining the results from the different simulations into one heat map
load(strcat(baseFoldername,'1','/densityMatrix.mat'));
fullDensity = meanDensity;
for i = 2:1:numOfSims
    currFoldername = strcat(baseFoldername,num2str(i));
    load(strcat(baseFoldername,num2str(i),'/densityMatrix.mat'));
    fullDensity = fullDensity + meanDensity;
end
fullDensity = fullDensity ./ numOfSims;
imagesc(fullDensity);
%% combining the results of the standard deviations to get a the characteristic correlation decay time
load(strcat(baseFoldername,'1','/densitySTDs.mat'));
fullSTDs = densitySTDs;
for i = 2:1:numOfSims
    currFoldername = strcat(baseFoldername,num2str(i));
    load(strcat(baseFoldername,num2str(i),'/densitySTDs.mat'));
    fullSTDs = fullSTDs + densitySTDs;
end
fullSTDs = fullSTDs ./ numOfSims;
t = 0:1/sampleRate:(length(fullSTDs)-1)*(1/sampleRate);
figure
plot(t,fullSTDs);
xlabel('time (sec)');
ylabel('Standard Deviation');
title('Averaged Standard Deviation of the particle densities');
%% Ron's version
numOfSims = 200; 
%numOfSims = 1;
ExperimentTime = 60; %sec
Dt = 1e-4;
N = ExperimentTime/Dt; % The number of steps to run the simulation for
sampleRate = 30; % fps
R = 0.75e-6; % radius in meter
T  = 298;
lx = 2.4e-5;
ly = 2.4e-5;
wallShrink = 0;
eta = 0.001;
numOfParticles = 1;
%wallShrink = 0;
%saveFoldername = 'F:\Experiments\Many body spreading experiment\Simulations'
MainFolder = pwd; %'F:\Experiments\Many body spreading experiment\Simulations\Single particle'
%baseFoldername = 'Parallel Harmonic ';
displayLive = 1;
%particlePositions = ManyBodyDiffusionSim(N,Dt,sampleRate,R,T,eta, lx, ly, numOfParticles, wallShrink,saveFoldername, displayLive)
parfor p=1%:numOfSims
    currFoldername = strcat(MainFolder,'\',num2str(p));
    ManyBodyDiffusionSim(N, Dt, sampleRate, R, T, eta, lx, ly, numOfParticles, wallShrink, currFoldername, displayLive);
end

%% Run one simulation   
numOfSims = 300

0; 
ExperimentTime = 60; %sec
Dt = 1e-4;
N = ExperimentTime/Dt;
sampleRate = 30;
R = 0.75e-6;
eta = 0.001;
lx = 0.75e-4;
ly = lx;
numOfParticles = 1;
wallShrink = 0;
currFoldername = pwd;
displayLive = 0;
T = 300;
%initposX = [1.12854917132807e-06; 2.23527532474401e-06; 1.13321873770091e-06; -1.17441058653142e-06; -2.19776133851434e-06; -1.10358946176809e-06]; % Xtraps locations from previous experiments
%initposY = [-1.99263553440715e-06; -1.23133715505655e-07; 1.96279242978705e-06; -2.03044832971602e-06; -1.19446040672767e-08; 2.05004525027041e-06]; % Ytraps locations
initposX = [1.12854917132807e-06];%; 2.23527532474401e-06];
initposY = [-1.99263553440715e-06];%; -1.23133715505655e-07];
initpos = [initposX, initposY];
saveFoldername = pwd;
hydro = 0;
for i = 1:numOfSims
    currFoldername = strcat(saveFoldername,'\',num2str(i,'%06d'));
    DiffSimRon(N,Dt,sampleRate,R,T,eta, lx, ly, numOfParticles,currFoldername, displayLive, initpos, hydro);
end
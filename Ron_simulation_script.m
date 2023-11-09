%% Run one simulation - brownian particles free diffusion
numOfSims = 7000;

ExperimentTime = 60; %sec
Dt = 1e-4;
%N = 3;
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
initposX = [1.12854917132807e-06]%; 2.23527532474401e-06];
initposY = [-1.99263553440715e-06]%; -1.23133715505655e-07];
initpos = [initposX, initposY];
saveFoldername = pwd;
hydro = 0;

for i = 1:numOfSims
    currFoldername = strcat(saveFoldername,'\',num2str(i,'%06d'));
    DiffSimRon(N,Dt,sampleRate,R,T,eta, lx, ly, numOfParticles,currFoldername, displayLive, initpos, hydro);
end
%% Simuations with or without reset - active and passive particles non instantaneous
clear;clc;
NumOfSims = 1;
numOfFastParticles = 6; %active paricles
numOfSlowParticles = 900; %passive particles
numOfresets = 6; % how many reset events
NumPar2Reset = 1; % particles to reset together (one by one or all together for example)
hologramstepsize = 0.5e-7; % the stepsize - when trap is moving
hologramspeed = 0.02; % the time between every hologram calculating
releasetogether = 0;%%1;
rate = 0.05; % rate of resetting (1/\mu)
ExperimentTime = 10; %sec
fps = 30; % frames per second

%g=exprnd(mu,1,NumRestart);
%g=round(g, 2);
Dt = 1e-4;
sampleRate = 30;
Rslow = 0.75e-6;
Rfast = Rslow;
Meter2Pixles = 15e6; % one micrometer is about 15 pixles on the screen 
Rpixels = Rslow*Meter2Pixles;
eta = 0.001;
%lx = 0.75e-4;
lx = 6e-5;
ly = lx;
multiplication = 1; % if you want one population to have different diffusion coefficient
numOfParticles = numOfSlowParticles+numOfFastParticles;
wallShrink = 0;
%currFoldername = pwd;
displayLive = 1;
mkvideo = 1;
T = 300;
initposXfast = [1.12854917132807e-06; 2.23527532474401e-06; 1.13321873770091e-06; -1.17441058653142e-06; -2.19776133851434e-06; -1.10358946176809e-06]; % Xtraps locations from previous experiments
initposYfast = [-1.99263553440715e-06; -1.23133715505655e-07; 1.96279242978705e-06; -2.03044832971602e-06; -1.19446040672767e-08; 2.05004525027041e-06]; % Ytraps locations
TrapPos = [initposXfast, initposYfast];
[slowparticlesX, slowparticlesY] = randomizePositionsTest2populations([-lx/2, lx/2], [-ly/2, ly/2], numOfSlowParticles, Rslow, TrapPos(1:numOfFastParticles,:), Rfast);
initPositions = zeros(numOfParticles, 2);
initposX = [initposXfast; slowparticlesX]; %; 2.23527532474401e-06];
initposY = [initposYfast; slowparticlesY];%; -1.23133715505655e-07];
initpos = [initposX, initposY];
saveFoldername = pwd;
hydro = 1;
activespeed = 3e-6;%[];%;
AngularVel = 0.5;%[];%1;
trapstiffness = 1e-7 ; % Newton/meter
%FavorDirSteps = 25000; % number of steps in favable direction
    for j = 1:NumOfSims
        SimsFolder = strcat(saveFoldername,'\SimNum_', num2str(j,'%06d'));
        N = ExperimentTime/Dt;
        particlePositions = DiffSimRonActivePassive(N,Dt,sampleRate,Rslow,T,eta, lx, ly, numOfParticles,SimsFolder, displayLive, initpos(1:numOfParticles,:), ...
            hydro, numOfFastParticles, multiplication, Rfast, ...
            activespeed, AngularVel, numOfresets, rate, ...
            trapstiffness, NumPar2Reset, hologramstepsize, hologramspeed, releasetogether, mkvideo);       
    end

save("rate.mat", "rate");
%% To be fixed - 20230926
numOfSims = 7000;

ExperimentTime = 10; %sec
Dt = 1e-4;
%N = 3;
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
initposX = [1.12854917132807e-06]%; 2.23527532474401e-06];
initposY = [-1.99263553440715e-06]%; -1.23133715505655e-07];
initpos = [initposX, initposY];
saveFoldername = pwd;
hydro = 0;

for i = 1:numOfSims
    currFoldername = strcat(saveFoldername,'\',num2str(i,'%06d'));
    DiffSimRon(N,Dt,sampleRate,R,T,eta, lx, ly, numOfParticles,currFoldername, displayLive, initpos, hydro);
end

%% Simuations with reset 
NumRestart=6; mu=20;
NumOfSims = 3000;
fps = 30; % frames per second
g=exprnd(mu,1,NumRestart);
g=round(g, 2);
Dt = 1e-4;
sampleRate = 30;
R = 0.75e-6;
Meter2Pixles = 15e6; % one micrometer is about 15 pixles on the screen 
Rpixels = R*Meter2Pixles;
eta = 0.001;
lx = 0.75e-4;
ly = lx;
numOfParticles = 6;
wallShrink = 0;
%currFoldername = pwd;
displayLive = 0;
T = 300;
initposX = [1.12854917132807e-06; 2.23527532474401e-06; 1.13321873770091e-06; -1.17441058653142e-06; -2.19776133851434e-06; -1.10358946176809e-06]; % Xtraps locations from previous experiments
initposY = [-1.99263553440715e-06; -1.23133715505655e-07; 1.96279242978705e-06; -2.03044832971602e-06; -1.19446040672767e-08; 2.05004525027041e-06]; % Ytraps locations
%initposX = [1.12854917132807e-06]%; 2.23527532474401e-06];
%initposY = [-1.99263553440715e-06]%; -1.23133715505655e-07];
initpos = [initposX, initposY];
saveFoldername = pwd;
hydro = 0;
TrapPos = initpos;
for j = 80:NumOfSims
    SimsFolder = strcat(saveFoldername,'\SimNum_',num2str(j,'%06d'));
    for i = 1:NumRestart
        ExperimentTime = g(i);
        N = ExperimentTime/Dt;
        currFoldername = strcat(SimsFolder,'\WalkNum_',num2str(i,'%06d'));
        %DiffSimRon(N,Dt,sampleRate,R,T,eta, lx, ly, numOfParticles,currFoldername, displayLive, initpos, hydro);
        particlePositions = DiffSimRon(N,Dt,sampleRate,R,T,eta, lx, ly, numOfParticles,currFoldername, displayLive, initpos, hydro);
        finalposX = particlePositions(end,:,1);
        finalposY = particlePositions(end,:,2);
        finalpositions = [finalposX', finalposY'];
        if i ~= NumRestart
        [initpos ResetD] = restartfxnExlodedArea(TrapPos, finalpositions, R);
        ResetD(1,3) = g(i);
        save(strcat(currFoldername,'\RestData_',num2str(i+1,'%06d'),'.mat'),"ResetD");
        end
    end
end


%% prin results
clear;clc;
lx = 0.75e-4;
ly = lx;
NumOfActivePar = 1;
xpos = load("pos_x.csv");
xpos = xpos';
ypos = load("pos_y.csv");
ypos = ypos';
Rslow = 0.75e-6;
Rfast = Rslow;
NumberOfSteps = length(xpos(1,:));
NumberOfParticles = length(xpos(:,1));
NumOfSlowParticles = NumberOfParticles - NumOfActivePar; 
colorFast = 'r';
colorSlow = 'b';
for i = 1:NumberOfSteps
    cla
    xlim([-lx/2 lx/2])
    ylim([-ly/2 ly/2])
    viscircles([xpos(1:NumOfActivePar, i), ypos(1:NumOfActivePar, i)], ...
        ones(NumOfActivePar,1).*Rfast, 'color', colorFast);
    hold on
    viscircles([xpos(1+NumOfActivePar:end, i), ypos(1+NumOfActivePar:end, i)], ...
        ones(NumOfSlowParticles,1).*Rslow, 'color', colorSlow);
    hold off
    pause(0.1)
end
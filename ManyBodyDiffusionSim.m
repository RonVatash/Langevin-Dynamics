function particlePositions = ManyBodyDiffusionSim(N,Dt,sampleRate,R,T,eta, lx, ly, numOfParticles, wallShrink,saveFoldername, Hydro, displayLive)
    tic
    startTime = toc;
    %% Setting the configuration for the simulation
    cfg = simConfig;
    cfg.N = N;
    cfg.Dt = Dt;
    cfg.sampleRate = sampleRate;
    cfg.numOfParticles = numOfParticles;
    cfg.eta = eta;
    cfg.R = ones(1,numOfParticles).*R;
    cfg.T = T;
    cfg.wallPositionsX = [-lx/2, lx/2 - wallShrink];
    cfg.wallPositionsY = [-ly/2, ly/2];
    cfg.xlimits=cfg.wallPositionsX.*1.1;
    cfg.ylimits=cfg.wallPositionsY*1.1;
    cfg.useWalls = false;
    cfg.wallRepulsionType = 'WCA';
    cfg.wallHarmonicK = 2e6*physconst('boltzmann')*cfg.T./1e-6;
%     cfg.wallGaussianA = 1e3*physconst('boltzmann')*cfg.T;
%     cfg.wallGaussianS = 10e-6;
    cfg.useParticleRepulsion = true;
    cfg.WCAEpsilon = 0.2*physconst('boltzmann')*cfg.T; % 0.2 kBT seems to work well
    cfg.useHydro = Hydro;
    cfg.displayLive = displayLive;
    cfg.saveFoldername = saveFoldername;
    cfg.savePeriod = cfg.N ./ 50; %save to disk
    cfg.useTraps = false;
    cfg.initTrapPositions = [0,0]%
%     cfg.A = -10*physconst('boltzmann')*cfg.T;
%     cfg.s = 6e-6;

    %% Creating the additional info object
    wallData = additionalData;
%     maxWallMoves = 100;
%     wallData.wallMoveSteps = nan(1,maxWallMoves);
%     wallData.newWallPositions = nan(1,maxWallMoves);
%     wallData.newWallPositions(1) = cfg.wallPositionsX(2);
%     wallData.wallMoveSteps(1) = 0;
%     wallData.wallMoveInd = 2;
%     wallData.wallCheckInd = 1;
%     wallData.closestParticlePositions = nan(cfg.N,1);
%     wallData.closestParticleDistances = nan(cfg.N,1);
    %% Randomizing particle starting positions
%     rng('shuffle');
    particlesX = [1.12854917132807e-06; 2.23527532474401e-06; 1.13321873770091e-06; -1.17441058653142e-06; -2.19776133851434e-06; -1.10358946176809e-06]; % Xtraps locations from previous experiments
    particlesY = [-1.99263553440715e-06; -1.23133715505655e-07; 1.96279242978705e-06; -2.03044832971602e-06; -1.19446040672767e-08; 2.05004525027041e-06]; % Ytraps locations
    initPositions = zeros(numOfParticles, 2);
    %initPositions(:,1) = particlesX;
    %initPositions(:,2) = particlesY;
    initPositions(1,1) = particlesX(1);
    initPositions(1,2) = particlesY(1);
    cfg.initPositions = initPositions;
%     if cfg.numOfParticles == 1
%         cfg.initPositions = cfg.initPositions';
%     end
    %% Running the simulation
    particlePositions = MDSim(cfg, @computeForces, @(x,y,z) false, @moveWall, @printCurrStep, wallData);
    
    totalTime = toc - startTime;
    %% Running the post-simulation analysis
%     postSimAnalysis(cfg, wallData);
function particlePositions = MDSimTest2populations(cfg, forcesFunc, feedbackCheckFunc, feedbackFunc, printFunc, addedData)
%% Basic definitions 
close all
if ~isempty(cfg.Nfast)
    fastparticales = cfg.Nfast;
    active = 1;
end
phi = rand(fastparticales, 1).*2*pi;

multiplication = cfg.multy;
kB = physconst('Boltzmann'); % Boltzmann constant [J/K]
R = cfg.R(1); % Currently works only with one constant R
gamma = 6*pi*R*cfg.eta; % friction coefficient
D = kB*cfg.T/gamma; % diffusion coefficient
%D = 1.29e-13; % mesured
d = 2; % The dimension of the problem. Currently ONLY WORKS FOR 2d!!
samplePeriod = round(1 / (cfg.Dt*cfg.sampleRate)); % Defines the samplePeriod at which we will sample the simulation. 
particlePositions = zeros(1+ceil(cfg.savePeriod/samplePeriod),cfg.numOfParticles, d); % Preallocate space
particlePositions(1,:,:) = cfg.initPositions; % sets the first step as defined in configuration
ActiveSpeed = cfg.Activespeed; % velocity of active particles;
if ~isempty(cfg.Rfast)
    RF = cfg.Rfast(1); % radius of fiat particles
else
    RF = R;
end
DrF = (kB*cfg.T)/(8*pi*RF^3*cfg.eta); % rotational diffusion coefficient
gammaF = 6*pi*RF*cfg.eta;
if ~isempty(cfg.multy)
    Dfast = kB*cfg.T/gammaF*multiplication;   % multiply by a multiplication factor
else
    Dfast = kB*cfg.T/gammaF;
end
  
%% Checking if the save directory already exists
if ~exist(cfg.saveFoldername, 'dir')
    mkdir(cfg.saveFoldername);
else
    answer = questdlg('Save folder already exists. Overwrite?');
    if strcmp(answer, 'Yes')
        delete(strcat(cfg.saveFoldername, '/pos_x.csv'));
        delete(strcat(cfg.saveFoldername, '/pos_y.csv'));
        delete(strcat(cfg.saveFoldername, '/cfg.m'));
    else
        error('Aborting sim to avoid overwriting the existing saved file');
    end
end
%% Saving the configuration
save(strcat(cfg.saveFoldername, '/cfg.mat'), 'cfg');
%% Setting up the run variables
currStepData = simStepData;
currStepData.particlePositions = cfg.initPositions;
if cfg.useWalls
    currStepData.wallPositionsX = cfg.wallPositionsX;
    currStepData.wallPositionsY = cfg.wallPositionsY;
end

if cfg.useTraps
    %currStepData.A = cfg.A;
    %currStepData.s = cfg.s;
    %trapPositions(1,:,:) = cfg.initTrapPositions';
    %currStepData.trapPositions = (squeeze(trapPositions(1,:,:)));
end
%cfg.useTraps = [];

%% Plotting the initial placements
if cfg.displayLive
    printFunc(cfg, currStepData, addedData);
    xlabel('x [m]');
    ylabel('y [m]');
  
    title('Initial placement');
    
    pause(0.01);
end
%% Checking whether to run hydrodynamic interactions
if ~cfg.useHydro
    Dx_slow = ones(cfg.numOfParticles-fastparticales,1).*D;
    Dx_fast = ones(fastparticales,1).*Dfast;
    Dx = [Dx_fast; Dx_slow];
    Dy = Dx;
end
%% Checking weather to run resetting interuaption
if cfg.resettimes
    v= cfg.hologramstepsize; dt_v= cfg.hologramspeed;
    dtN = dt_v/cfg.Dt;
    currStepData.SpringConst = cfg.SpringConst;
    currStepData.trapPositions = cfg.initTrapPositions;
    ResetTimes = cfg.resettimes;
    ResetTimesStep = round(ResetTimes./cfg.Dt); 
    numOfResets = length(ResetTimes);
    cfg.N = sum(ResetTimesStep) + (3*numOfResets/cfg.Dt);
    ResetTime = ResetTimesStep(1); %ResetTimesStep(1);
    ResetNum = 1;
    currStepData.ResetNum = ResetNum;
end
%% Check whether to make a video 
if cfg.mkvideo
        VideoSim = VideoWriter("Video", "Uncompressed AVI"); % temporary - to be changed
        open(VideoSim)
end
%% Running the simulation
sampleInd = 2;
for i = 2:1:cfg.N
   
    currStepData.stepNum = i;
    %% Printing data
    if mod(i, samplePeriod) == 0
        if cfg.displayLive
            printFunc(cfg, currStepData, addedData);
        end
        particlePositions(sampleInd,:,:) = currStepData.particlePositions;
        sampleInd = sampleInd + 1;
        if cfg.saveframes
           saveas(figure(1), strcat(cfg.saveFoldername, '\stepnum', num2str(i), ".png"));
        end
        if cfg.mkvideo
            FrameData = getframe(figure(1));
            FrameMat = FrameData.cdata;
            writeVideo(VideoSim, FrameMat)
        end
    end
    %% Saving the steps according to the save period parameter
    if mod(i, cfg.savePeriod) == 0 || i == cfg.N
        %% Saving the latest particle positions in .csv files
        dlmwrite(strcat(cfg.saveFoldername, '/pos_x.csv'),...
            particlePositions(1:sampleInd-1,:,1),...
            '-append');
        dlmwrite(strcat(cfg.saveFoldername, '/pos_y.csv'),...
            particlePositions(1:sampleInd-1,:,2),...
            '-append');
        %writematrix(particlePositions(1:sampleInd-1,:,1),strcat(cfg.saveFoldername, '/pos_x.csv'));%,'append');
        %writematrix(particlePositions(1:sampleInd-1,:,2),strcat(cfg.saveFoldername, '/pos_y.csv'));%,'append');    
           
        particlePositions(1:sampleInd-1,:,:) = 0;

        %% Saving the additional tracked data in a .mat file.
        if exist(strcat(cfg.saveFoldername, '/data.mat'), 'file')
            delete(strcat(cfg.saveFoldername, '/data.mat'));
        end
        save(strcat(cfg.saveFoldername, '/data.mat'), 'addedData');
        
        sampleInd = 1;
        strcat(cfg.saveFoldername, ' - ', num2str(100*(i/cfg.N)), '%');
    end
    %% choose particle and trap to reset and calculate trap positions
    if ~isempty(cfg.resettimes)
        if i == ResetTime
           cfg.useTraps = 1;
           
           %sprintf("Start Laser")
           NumOfParticles2Reset = cfg.NumPar2Reset; 
           ArrivedCheck = zeros(NumOfParticles2Reset, 1);
           [currStepData.resetParticle currStepData.trapidx] = chooseparticleandtrap(cfg.initTrapPositions , currStepData.particlePositions(1:fastparticales,1:2), NumOfParticles2Reset);
           %[currStepData.resetParticle currStepData.trapidx currStepData.traptheta] = chooseparticleandtrap(cfg.initTrapPositions , currStepData.particlePositions(1:fastparticales,1:2), NumOfParticles2Reset);
           ResetpartIdx = currStepData.resetParticle;
           Trapidx = currStepData.trapidx;
           xf = currStepData.particlePositions(ResetpartIdx,1); yf = currStepData.particlePositions(ResetpartIdx,2);
           xt = currStepData.trapPositions(Trapidx,1); yt = currStepData.trapPositions(Trapidx,2); 
           xCurr = xf; yCurr = yf;
           yDisp = yt - yf; xDisp = xt - xf;
           theta = atan2(yDisp,xDisp);
           vx = v .* cos(theta); vy = v .* sin(theta);
           sigX = sign(xDisp); sigY = sign(yDisp);
           currSigX = sigX; currSigY = sigY;
        end
        
        if cfg.useTraps
    
            if ~mod(i - ResetTime, round(dtN))
               for  RP = 1:NumOfParticles2Reset % RP stands for reset particl
                    if (currSigX(RP) == sigX(RP) && currSigY(RP) == sigY(RP))
                        xCurr(RP) = xCurr(RP) + vx(RP); yCurr(RP) = yCurr(RP) + vy(RP);
                        currSigX(RP) = sign(xt(RP) - xCurr(RP)); currSigY(RP) = sign(yt(RP) - yCurr(RP)); 
                        if (currSigX(RP) ~= sigX(RP) || currSigY(RP) ~= sigY(RP)) %if passed the target
                            xCurr(RP) = xt(RP);
                            yCurr(RP) = yt(RP);
                            ArrivedCheck(RP) = 1;
                            if ResetNum <= numOfResets && isempty(setdiff(xCurr, xt))
                            ResetTime = i + ResetTimesStep(ResetNum);
                            end
                        end
                        currStepData.CurtrapPositions = [xCurr yCurr];

                    end
               end
            end
        end
    end
    %% Forces computation
    [fx, fy] = ...
        forcesFunc(cfg, currStepData, addedData);
    %% Check whether to close the laser (relevant for resetting experiment)
  if cfg.useTraps && ~isempty(cfg.resettimes)
      if cfg.releasetogether || NumOfParticles2Reset == 1
          if isempty(setdiff(xCurr, xt)) && isempty(setdiff(yCurr, yt))%xCurr == xt && yCurr == yt
           cfg.useTraps = 0;
           ResetNum = ResetNum + 1;
           currStepData.ResetNum = ResetNum;
           %sprintf("Close laser")
          end
      else
          arrivedidx = find(ArrivedCheck == 1);
          if arrivedidx
             NumOfParticles2Reset = NumOfParticles2Reset - length(arrivedidx);
             currSigX(arrivedidx) = []; currSigY(arrivedidx) = [];
             sigX(arrivedidx) = []; sigY(arrivedidx) = [];
             vx(arrivedidx) = []; vy(arrivedidx) = [];
             xCurr(arrivedidx) = []; xt(arrivedidx) = [];
             yCurr(arrivedidx) = []; yt(arrivedidx) = [];
             ArrivedCheck(arrivedidx) = [];
             currStepData.resetParticle(arrivedidx) = [];
             currStepData.CurtrapPositions = [xCurr yCurr];
             currStepData.trapidx(arrivedidx) = [];
          end
         
      end
  end

    %% Hydrodynamic interactions computation
    if cfg.useHydro
        %DMat =
        %rtd_correct_line_traps_rotne_prager(currStepData.particlePositions,R, D);
       if active
        DMat = rotnePrager_correctedRonActivePassive(currStepData.particlePositions, D, Dfast, cfg);
       else
        DMat = rotnePrager(currStepData.particlePositions, R, D);
       end
        rootMat = chol(DMat,'lower');
        % Note that chol may be used only since it is multiplied by a set
        % of gaussian random variables. See https://stats.stackexchange.com/questions/83850/confused-about-cholesky-and-eigen-decomposition
        DVec = DMat*ones(cfg.numOfParticles*d,1);
        AVec = rootMat*randn(cfg.numOfParticles*d,1);
        Dx = DVec(1:2:end);
        Dy = DVec(2:2:end);
        Ax = AVec(1:2:end);
        Ay = AVec(2:2:end);
    else
        if isempty(cfg.Nfast)
           Ax = randn(cfg.numOfParticles,1).*sqrt(2*D);
           Ay = randn(cfg.numOfParticles,1).*sqrt(2*D);
        else
            Ax_slow = randn(cfg.numOfParticles-fastparticales,1).*sqrt(2*D);
            Ay_slow = randn(cfg.numOfParticles-fastparticales,1).*sqrt(2*D);
            Ax_fast = randn(fastparticales,1).*sqrt(2*Dfast);
            Ay_fast = randn(fastparticales,1).*sqrt(2*Dfast);
            % for FastPar = 1:fastparticales
            %      if cfg.N == 1    % choose on preferable direction
            %          FavorDirx = rand(1);
            %          xDir = round(FavorDirx);
            %          FavorDiry = rand(1);
            %          yDir = round(FavorDiry);
            %      elseif ~mod(cfg.N, cfg.numOfstepsInFavorDir)
            %          FavorDirx = rand(1);
            %          xDir = round(FavorDirx);
            %          FavorDiry = rand(1);
            %          yDir = round(FavorDiry);
            %      end
            %      if xDir
            %          Ax_fast(FastPar,1) = abs(Ax_fast(FastPar, 1));
            %      elseif ~xDir
            %          Ax_fast(FastPar,1) = (-1)*abs(Ax_fast(FastPar, 1));
            %      end
            %      if yDir
            %          Ay_fast(FastPar,1) = abs(Ay_fast(FastPar, 1));
            %      elseif ~yDir
            %          Ay_fast(FastPar,1) = (-1)*abs(Ay_fast(FastPar, 1));
            %     end
            % end
            Ax = [Ax_fast; Ax_slow];
            Ay = [Ay_fast; Ay_slow];
        end
    end
    %% Running the step
    currStepData.particlePositions(:,1) = currStepData.particlePositions(:,1) +...
        Ax.*sqrt(cfg.Dt)+...
        (Dx./(kB.*cfg.T)).*fx.*cfg.Dt;  
    currStepData.particlePositions(:,2) = currStepData.particlePositions(:,2) +...
        Ay.*sqrt(cfg.Dt) +...
        (Dy./(kB.*cfg.T)).*fy.*cfg.Dt;
    %% add favorable direction
    if  ~isempty(cfg.Activespeed)
        currStepData.particlePositions(1:fastparticales,1) = currStepData.particlePositions(1:fastparticales,1) + cfg.Activespeed*cos(phi)*cfg.Dt; 
        currStepData.particlePositions(1:fastparticales,2) = currStepData.particlePositions(1:fastparticales,2) + cfg.Activespeed*sin(phi)*cfg.Dt;
        Aphi = randn(cfg.Nfast,1).*sqrt(2*DrF);
        if isempty(cfg.Omega)
            Omega = 0;
        else
            Omega = cfg.Omega;
        end
        phi = phi + Omega*cfg.Dt + Aphi.*sqrt(cfg.Dt);
     
    end
    %% Checking whether to apply a feedback to the system, and applying it if necessary
    if feedbackCheckFunc(cfg, currStepData, addedData)
        feedbackFunc(cfg, currStepData, addedData);
    end



end
%% checking whether to stop make movie
    if cfg.mkvideo
        close(VideoSim)
    end

%% Reading the output file
   particlePositionsX = dlmread(strcat(cfg.saveFoldername, '/pos_x.csv'));
   particlePositionsY = dlmread(strcat(cfg.saveFoldername, '/pos_y.csv'));
   %particlePositionsX = readmatrix(strcat(cfg.saveFoldername, '/pos_x.csv'));
   %particlePositionsY = readmatrix(strcat(cfg.saveFoldername, '/pos_y.csv'));
   particlePositions = zeros(size(particlePositionsX, 1), cfg.numOfParticles, 2);
   particlePositions(:,:,2) = particlePositionsY;
   particlePositions(:,:,1) = particlePositionsX;

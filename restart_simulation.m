%% Experiment - constant velocisy - Ron Simulation
%cam = webcam('Logitech BRIO');
% Ron Vatash 13.10.22
NumRestart=20; mu=7;
particleRadius = 10; % particle radius in pixels
fps = 30; % frames per second
g=exprnd(mu,1,NumRestart);
F = dir('*_F');
frames = subdir(F(6).name);
%frames = frames(1:2000);
Power=3.0;
ptime=1;dt_v=0.03;
locationsAtRestart = zeros(NumRestart,2);restartTimes = zeros(NumRestart,2);
xy_intial=zeros(NumRestart,2); % Choosen particle position at the time of reseting
resetting_pos=zeros(NumRestart,2); % Choosen trap position
%particleCenter = zeros(2,2);
%vid.TriggerRepeat = Inf;
%src.FrameRate = 30; start(vid);
%tic
currPic = frames(1).name;
currPic = imread(currPic);
particleCenters = findparticle(currPic,1,21,0.3,21,23,0); % find initial positions (trap positions)
xt1=particleCenters(1,1); yt1=particleCenters(1,2);     %xt,yt particle location at the beggining of the experiment.
intial_pos(1,1)=xt1;intial_pos(1,2)=yt1;
xt2=particleCenters(2,1);yt2=particleCenters(2,2);
intial_pos(2,1)=xt2;intial_pos(2,2)=yt2;
xt3=particleCenters(3,1); yt3=particleCenters(3,2);     %xt,yt particle location at the beggining of the experiment.
intial_pos(3,1)=xt3;intial_pos(3,2)=yt3;
xt4=particleCenters(4,1);yt4=particleCenters(4,2);
intial_pos(4,1)=xt4;intial_pos(4,2)=yt4;
xt5=particleCenters(5,1); yt5=particleCenters(5,2);     %xt,yt particle location at the beggining of the experiment.
intial_pos(5,1)=xt1;intial_pos(5,2)=yt5;
xt6=particleCenters(6,1);yt6=particleCenters(6,2);
intial_pos(6,1)=xt6;intial_pos(6,2)=yt6;
xts = particleCenters(:,1); % make a list of all x traps positions 
yts = particleCenters(:,2); %[yt1; yt2; yt3; yt4; yt5; yt6];
%(the (0,0) of the experiment, the particle will be restarted to here)
%loop for running experiment
% close laser
colormap gray %plot 
imagesc(currPic);
hold on
plot(xts, yts, '.','markersize',10,'color','r');
hold off
tic % save the diffusion time between every reset event (before reseting)
reset_times = [];
for i=1:NumRestart
    sprintf(strcat('close laser and Diffuse for - ',num2str(g(i)),'seconds'))
    i; % the number of the restart
    pause(g(i)); %diffusion for g(i) time
    reset_times = [reset_times; g(i)];
    currPic = frames(round(sum(reset_times(:))*fps)).name;
    %currPic = frames(randi([200 1000],1)).name; 
    currPic = imread(currPic);
    centers = findparticle(currPic,1,21,0.3,21,23,0);
   % if length(centers)< 6 %stop experiment if lost one of the particles
    %    break
    %end
    R = randi(length(centers(:,1))); % randomly select the particle to be resetted
    xf= centers(R,1);% corresponds to x coordinate of stk
    yf= centers(R,2);% corresponds to y coordinate of stk
    %  check which center is nearest to xf and yf
    r_trap1=sqrt((xf-xt1).^2+(yf-yt1).^2);
    r_trap2=sqrt((xf-xt2).^2+(yf-yt2).^2);
    r_trap3=sqrt((xf-xt3).^2+(yf-yt3).^2);
    r_trap4=sqrt((xf-xt4).^2+(yf-yt4).^2);
    r_trap5=sqrt((xf-xt5).^2+(yf-yt5).^2);
    r_trap6=sqrt((xf-xt6).^2+(yf-yt6).^2);
    r_traps = [abs(r_trap1); abs(r_trap2); abs(r_trap3); abs(r_trap4); abs(r_trap5); abs(r_trap6)]; % make vector of all distances from final position
    indrmin = find(r_traps == min(r_traps)); % find the index of closest trap
    xt = xts(indrmin); yt = yts(indrmin); % give me the position of the nearest trap
    resetting_pos(i,1)=xt;resetting_pos(i,2)=yt; % document the position of the trap in the i'th reset event
    locationsAtRestart(i,:) = [xf,yf]; %save them - the locations of the choosen particle at the i'th reset event
    colormap gray %plot 
    imagesc(currPic);
    title(strcat('reset event number',num2str(i,'%03d'))); 
    radius = 3;
    viscircles([xf,yf],radius, 'Color','b'); % Show me particle position
    viscircles([xt,yt],radius, 'Color','r'); % Show me trap position
    v=1.5; % This is velocity for restart at constant velocity if using restart at constant velocity please comment at top v_base and r_base.
    %v=vbase*(sqrt((xf-xt)^2+(yf-yt)^2)/rbase); %the speed for this restart (changes each restart because constant time) 
    % actual restart:
    restartTimes(i,1) = toc; %save time before resetting
    xCurr = xf; yCurr = yf;
    yDisp = yt - yf; xDisp = xt - xf;
    theta = atan2(yDisp,xDisp);
    vx = v .* cos(theta); vy = v .* sin(theta);
    sigX = sign(xDisp); sigY = sign(yDisp);
    currSigX = sigX; currSigY = sigY;
    
    CheckifParticleOnTheLine = particleontheline(centers,particleRadius, R, xf, yf, xt, yt);
    if CheckifParticleOnTheLine == 0 % Checking if there is a particle on the way to the trap - if there is one - NO RESET
        sprintf('open laser') % Add to try
        while (currSigX == sigX && currSigY == sigY)
            xCurr = xCurr + vx; yCurr = yCurr + vy;
            viscircles([xCurr,yCurr],radius, 'Color','b');
            viscircles([xt,yt],radius, 'Color','r');
            currSigX = sign(xt - xCurr); currSigY = sign(yt - yCurr); 
            if (currSigX ~= sigX || currSigY ~= sigY) %if passed the target
                xCurr = xt;
                yCurr = yt;
            end
            pause(dt_v); %small pause between each mask
        end
    else
        pause(ptime); %pause time after resetting
        xy_intial(i,:)=[xCurr,yCurr]; % add here the last position coordinates feed to slm
        restartTimes(i,2) = toc; %save time after resetting first columb will be the time reset event start and second columb will be time when reset event end the delta will be draging time
        sprintf('Steric Disturbance - No RESET')
    end
end
sprintf('close laser')
%end of experiment start of save
Experiment_DATA = [g', locationsAtRestart, resetting_pos]; % Matrix contains all data - 1 column: reset time, 2-3: xpos and ypos at restart of the choosen particle, 4-5: xpos ypos of trap

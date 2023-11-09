function [particleidx trapidx theta] = chooseparticleandtrap(initialTrapPos, initialParticlePos, NumOfParticles2reset)
% this function choose and reset one particale to its initial position
% this function should work with Ron_simulation script
%lx = 7.5e-6; ly = lx;
%NumOfParticles2reset = 1;
particleidx = zeros(NumOfParticles2reset, 1);
trapidx = particleidx;
xts = initialTrapPos(:,1); % make a list of all x traps positions 
yts = initialTrapPos(:,2); %[yt1; yt2; yt3; yt4; yt5; yt6];
choosefromvec = [1:length(xts(:,1))]; % initial numer of traps to choose one from

for i = 1:NumOfParticles2reset
    R = randi(length(choosefromvec)); % randomly select the particle to be reset
    choosenparticle = choosefromvec(R); 
    particleidx(i) = choosenparticle;
    choosefromvec = setdiff(choosefromvec, choosenparticle);
    xf= initialParticlePos(choosenparticle,1);% corresponds to x coordinate of stk at the begginig of the resetting event
    yf= initialParticlePos(choosenparticle,2);% corresponds to y coordinate of stk at the begginig of the resetting event
    Xdiff = xf*ones(length(xts), 1) - xts;
    Ydiff = yf*ones(length(yts), 1) - yts;
    rtraps = sqrt(Xdiff.^2+Ydiff.^2);  % distance from particle to traps
    possibletrapidx = find(rtraps == min(rtraps)); 
    %trapidx(i) = % find the index of closest trap
    checkifemptytrap = setdiff(possibletrapidx, trapidx);
    exludedminmatrix = rtraps;
    while isempty(checkifemptytrap)
          exludedminmatrix = setdiff(exludedminmatrix, min(exludedminmatrix));
          possibletrapidx = find(rtraps == min(exludedminmatrix));
          checkifemptytrap = setdiff(possibletrapidx, trapidx);
    end
    trapidx(i) = possibletrapidx;
    %theta = atan2(Ydiff(trapidx),Xdiff(trapidx));
end
    %Rs = zeros(length(xts),1);
%Rs = Rs+radius;
%rtraps = zeros(length(xts),1);



%xt = xts(indrmin); yt = yts(indrmin);  
%resetting_pos = [xt, yt]; % final position of choosen trap
%locationsAtRestart = [xf,yf];
%xCurr = xf; yCurr = yf;
%CheckifParticleOnTheLine = particleontheline_2(Finalpos,radius, R, xf, yf, xt, yt);
%ResetData = [xt, yt, 0, CheckifParticleOnTheLine]; % trap pos x trap pos y reset time (comes from the simulation) interuption or not 1 or zero
%particleidx = R;
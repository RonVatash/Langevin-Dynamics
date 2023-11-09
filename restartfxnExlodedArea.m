function [initpos ResetData] = restartfxnExlodedArea(TrapPos, Finalpos, radius)
% this function choose and reset one particale to its initial position
% this function should work with Ron_simulation script
%lx = 7.5e-6; ly = lx;
xts = TrapPos(:,1); % make a list of all x traps positions 
yts = TrapPos(:,2); %[yt1; yt2; yt3; yt4; yt5; yt6];
R = randi(length(xts(:,1))); % randomly select the particle to be resetted
xf= Finalpos(R,1);% corresponds to x coordinate of stk
yf= Finalpos(R,2);% corresponds to y coordinate of stk
Rs = zeros(length(xts),1);
Rs = Rs+radius;
rtraps = zeros(length(xts),1);

for i = 1:length(rtraps)
    rtraps(i) = sqrt((xf-xts(i)).^2+(yf-yts(i)).^2); % create a distance from traps to choosen particle vector
end
indrmin = find(rtraps == min(rtraps)); % find the index of closest trap
xt = xts(indrmin); yt = yts(indrmin);
resetting_pos = [xt, yt];
locationsAtRestart = [xf,yf];
xCurr = xf; yCurr = yf;
CheckifParticleOnTheLine = particleontheline_2(Finalpos,radius, R, xf, yf, xt, yt);
ResetData = [xt, yt, 0, CheckifParticleOnTheLine]; % trap pos x trap pos y reset time (comes from the simulation) interuption or not 1 or zero

if CheckifParticleOnTheLine
    initpos = Finalpos;
    initpos(R,1) = xt;
    initpos(R,2) = yt;
    sprintf("Reset Event Succeed")
else
    initpos = Finalpos; 
    sprintf('Steric Disturbance - No RESET')
end 


end

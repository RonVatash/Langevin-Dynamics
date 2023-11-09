function [fx,fy]=getTrapForcesRon(currStepData, cfg)
%%Calculate the forces all of the traps apply on all of the particles.
K = cfg.SpringConst;
Radi = cfg.Rfast(1);
theta = currStepData.traptheta;
particlePositions = currStepData.particlePositions;
CurtrapPositions =  currStepData.CurtrapPositions;
ResetParticleIdx = currStepData.resetParticle;
%TrapIdx = currStepData.trapidx;
ChoosenParPos = particlePositions(ResetParticleIdx,:);
%trapPositions = trapPositions(currStepData.trapidx, :);
% CurtrapPosX = ChoosenParPos(:,1) + Radi*cos(theta);  
% CurtrapPosY = ChoosenParPos(:,2) + Radi*sin(theta);  
% CurtrapPos = [CurtrapPosX CurtrapPosY];


numOfParticles = size(particlePositions,1);
numOfTraps = size(CurtrapPositions,1);

%If A,s (depth, size of the trap) are scalars convert them into a matrix.
if numel(K) == 1
    K = K.*ones(numOfTraps,2);
end
%if numel(s) == 1
 %   s = s.*ones(numOfTraps,2);
%end

fx = zeros(numOfParticles,1);
fy = zeros(numOfParticles,1);


total_force = getTrapForcesForSingleParticleRon(ChoosenParPos, CurtrapPositions, K);


fx(ResetParticleIdx) = fx(ResetParticleIdx) + total_force(:,1);
fy(ResetParticleIdx) = fy(ResetParticleIdx) + total_force(:,2);

% for currParticle = 1:numOfParticles
%     checkedParticlePosition = particlePositions(currParticle,:);
% 
%     total_force = getTrapForcesForSingleParticle(checkedParticlePosition, trapPositions, A, s);
% 
%     fx(currParticle) = total_force(1);
%     fy(currParticle) = total_force(2);
% end
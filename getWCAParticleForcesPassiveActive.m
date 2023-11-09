function [fx,fy]=getWCAParticleForcesPassiveActive(particlesMat, cfg)
% getWCAParticleForces computes the forces on the particles from the other
% particles using WCA repulsion (effectively similar to hard sphere)
% particlesMat: The positions of the particles in the current step
% R: The radii of the different particles (Only works for a single R
% currently)
epsilon = cfg.WCAEpsilon;
Rpassive = cfg.R(1);
numOfParticles = size(particlesMat,1);
fx = zeros(numOfParticles,1);
fy = zeros(numOfParticles,1);
if isempty(cfg.Nfast)
    numOfFastParticles = 0;
else
    numOfFastParticles = cfg.Nfast;
end

if isempty(cfg.Rfast)
    Ractive = Rpassive;
else
    Ractive = cfg.Rfast(1);
end

%R = (Rpassive+Ractive)/2; % added by ron 20230928
 % added by ron 20230928
%Rc = 2.5*R; %  cut-off distance - beyond that distance, the force and the potential vanishes. - needs to be decided

for currParticle = 1:numOfParticles
    checkedParticlePosition = particlesMat(currParticle,:);
    if currParticle <= numOfFastParticles
        RcurrParticle = Ractive;
    else
        RcurrParticle = Rpassive;
    end
    %% For each particle, checking its interaction with all others
    for pairedParticle = (currParticle+1):numOfParticles
        if pairedParticle <= numOfFastParticles
            RpairedParticle = Ractive;
        else
            RpairedParticle = Rpassive;
        end

        pairedParticlePosition = particlesMat(pairedParticle,:);
        yDiff = checkedParticlePosition(2)-pairedParticlePosition(2);
        xDiff = checkedParticlePosition(1)-pairedParticlePosition(1);
        
        sigma = RpairedParticle + RcurrParticle;
        Rc = (2^(1/6))*sigma;
        distance = sqrt(xDiff^2+yDiff^2);
        %distance = sqrt(xDiff^2+yDiff^2);  % changed by Ron - deleted the R, I don't think it should be there, the positions are the position of the ceter of the sphere - Ask Yael.
        if distance <= Rc
            %% Computing the forces between the particles
            totalForce = getLJForceRon(distance, epsilon, sigma);
            angle = atan2(yDiff,...
                          xDiff);

            %% Adding the force to both particles with opposite signs
            fx(currParticle) = fx(currParticle) + totalForce.*cos(angle);
            fx(pairedParticle) = fx(pairedParticle) - totalForce.*cos(angle);
            fy(currParticle) = fy(currParticle) + totalForce.*sin(angle);
            fy(pairedParticle) = fy(pairedParticle) - totalForce.*sin(angle);
        end
    end
end
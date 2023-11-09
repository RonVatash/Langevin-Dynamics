function [particlesX,particlesY] = randomizePositionsTest2populations(wallPositionsX, wallPositionsY, numOfParticles, R, trappos, RadiusFast)
% Randomizes the particle positions while ensuring no particles overlap
numOfFastParticles = length(trappos(:,1));
particlesX = zeros(numOfParticles + numOfFastParticles,1);
particlesX(1:numOfFastParticles) = trappos(:,1);
particlesY = zeros(numOfParticles + numOfFastParticles,1);
particlesY(1:numOfFastParticles) = trappos(:,2);
lx = abs(wallPositionsX(2) - wallPositionsX(1));
ly = abs(wallPositionsY(2) - wallPositionsY(1));
for currParticle = 1+numOfFastParticles:numOfParticles+numOfFastParticles
    collision = true;
    %Continuing until there is no collision
    while collision
        collision = false;
        % Randomizing a position for the current particle
        particlesX(currParticle) = (lx-2.*R).*rand(1,1) + wallPositionsX(1) + R;
        particlesY(currParticle) = (ly-2.*R).*rand(1,1) + wallPositionsY(1) + R;
        for collidedParticle = 1:currParticle
            if collidedParticle <= numOfFastParticles
                R2 = RadiusFast;
            else
                R2 = R;
            end
            currCollision = ...
                checkCollision_Gilad(particlesX(currParticle), particlesY(currParticle), R,...
                particlesX(collidedParticle),particlesY(collidedParticle), R2);
            if currCollision
                collision = true;
                break
            end
        end
    end
end
particlesX = particlesX(1+numOfFastParticles:end);
particlesY = particlesY(1+numOfFastParticles:end);

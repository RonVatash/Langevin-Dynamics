function total_force = getTrapForcesForSingleParticleRon(particlePosition, trapPosition, K)
%Calculate force from all optical traps on a particle

f = -K.* (particlePosition-trapPosition);

%f = K .* (particlePosition-trapPositions) .* exp((-0.5*(particlePosition-trapPositions)./s).^2);
%total_force = sum(f,1);
total_force = f; %sum(f,1);
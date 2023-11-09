function [fx, fy] = computeForcesActivePassive(cfg, currStepData, addedData)
% this function is a mofdified function of Gilads function - modified by
% ron 20230928
% the purpose was to calculate the force acting between two particles
% according to WCA potential with a cut-off distance of sigma*2^(1/6) - at
% this distance the effective potential is about zero. 
fx = 0;
fy = 0;
if cfg.useWalls
    if strcmp(cfg.wallRepulsionType, 'WCA')
        [fxWall, fyWall] = getWCAWallForces(currStepData.particlePositions,...
                                            cfg.R(1),...
                                            cfg.WCAEpsilon,...
                                            currStepData.wallPositionsX,...
                                            currStepData.wallPositionsY);
    elseif strcmp(cfg.wallRepulsionType, 'Harmonic')
        [fxWall, fyWall] = getHarmonicWallForces(currStepData.particlePositions,...
                                            cfg.R(1),...
                                            cfg.wallHarmonicK,...
                                            currStepData.wallPositionsX,...
                                            currStepData.wallPositionsY);
    elseif strcmp(cfg.wallRepulsionType, 'Gaussian')
        [fxWall, fyWall] = getGaussianWallForces(currStepData.particlePositions,...
                                            cfg.R(1),...
                                            cfg.wallGaussianA, cfg.wallGaussianS,...
                                            currStepData.wallPositionsX,...
                                            currStepData.wallPositionsY);
    else
        error(strcat('Unknown wall repulsion type ', cfg.wallRepulsionType));
    end
    fx = fx + fxWall;
    fy = fy + fyWall;
end

if cfg.useParticleRepulsion
    [fxParticles, fyParticles] = ...
        getWCAParticleForcesPassiveActive(currStepData.particlePositions, cfg);
    fx = fx + fxParticles;
    fy = fy + fyParticles;
end

if cfg.useTraps
    [fxTraps, fyTraps] = getTrapForcesRon(currStepData, cfg);
    fx = fx + fxTraps;
    fy = fy + fyTraps;
end
classdef simStepData < handle
    properties
        stepNum % The numer of the current step in the sim
        
        particlePositions % The positions of the particles in the current step
        trapPositions % The positions of the traps in the current step
        A % The depth of the trap
        s % The size of the trap
        
        wallPositionsX % The x positions of the walls, if applicable
        wallPositionsY % The y positions of the walls, if applicable
        resetParticle; % Choosen particle index
        trapidx; % Nearest trap index
        SpringConst % stiffness
        traptheta % angle between trap position to particle position
        CurtrapPositions % current trap position - changes with time
        ResetNum
    end
end
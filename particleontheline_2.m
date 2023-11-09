function CheckifParticleOnTheLine = particleontheline_2(centers, radius, choosenPart, xf, yf, xt, yt) 
% This function checks if there is a particle on the way to the trap before
% resetting. 
% 1 - No particle on the way 0 - There is a particle - No reset
% Ron Vatash 17.10.22
v = 1.5;
R = choosenPart;
dr = radius/2;
yCurrs = []; xCurrs = [];
xCurr = xf; yCurr = yf;
yDisp = yt - yf; xDisp = xt - xf;
theta = atan2(yDisp,xDisp);
vx = v .* cos(theta); vy = v .* sin(theta);
rvector = sqrt(xDisp^2+yDisp^2);
rvalues = [0:dr:rvector];
xCurrs = rvalues*cos(theta);
yCurrs = rvalues*sin(theta);
%Xs = xCurrs(10:20)
%Ys = yCurrs(10:20)
ParticleIndex = [1:length(centers)]; % the index of particles
ParticleIndex_NoR = setdiff(ParticleIndex, R); % the indeces of the particles that are no
indicators = [];
    for i = 1:length(ParticleIndex_NoR)
        centers_particle_i = centers(ParticleIndex_NoR(i),:);
        wx = find(xCurrs-2*radius < centers_particle_i(1,1) & centers_particle_i(1,1) < xCurrs+2*radius);
        wy = find(yCurrs-2*radius < centers_particle_i(1,2) & centers_particle_i(1,2) < yCurrs+2*radius);
        if isempty(wx)
            Check = 0;
        elseif isempty(wy)
            Check = 0;
        else
            diff = setdiff(wx,wy);
            if length(diff) == length(wx)
                Check = 0;
            else
      
                for j = 1:length(wx)
                    distance = sqrt((xCurrs(wx(j))-centers_particle_i(1,1))^2+(yCurrs(wx(j))-centers_particle_i(1,2))^2);
                    if distance < 2*radius
                        CheckifParticleOnTheLine = 1;
                        
                    return
                    else 
                        Check = 0;
                    end
                end

            end
           
        end
        indicators = [indicators; Check]; 
    end
    if sum(indicators) == 0
       CheckifParticleOnTheLine = 1; % one if there are NO Particles on The Way - Reset Event Succeed
    else
        CheckifParticleOnTheLine = 0; 
    end
    
end
trackfileAll = dir("tr_all*.mat");
fps = 30;
initposX = [1.12854917132807e-06; 2.23527532474401e-06; 1.13321873770091e-06; -1.17441058653142e-06; -2.19776133851434e-06; -1.10358946176809e-06]; % Xtraps locations from previous experiments
initposY = [-1.99263553440715e-06; -1.23133715505655e-07; 1.96279242978705e-06; -2.03044832971602e-06; -1.19446040672767e-08; 2.05004525027041e-06]; % Ytraps locations
trackfileAll = struct2array(load(trackfileAll(1).name)); 
ParticleNames = unique(trackfileAll(:,6));
NumOfSteps = length(unique(trackfileAll(:,5)));
radius = 0.75; % micrometer
lx = 30*radius;
ResetEventsIdx = find(trackfileAll(:,8) == 1);
ResetEventSteps = trackfileAll(ResetEventsIdx, 5);
ResetEventParticle = trackfileAll(ResetEventsIdx, 6);
TrapLocsAtReset = [trackfileAll(ResetEventsIdx, 1) trackfileAll(ResetEventsIdx, 2)];
ResetData = [TrapLocsAtReset, ResetEventSteps, ResetEventParticle];
for i = 550000:NumOfSteps
    checkResetBefore = find(ResetEventSteps == i+1);
    checkResetAfter = find(ResetEventSteps == i);
    stepframe = trackfileAll(trackfileAll(:,5)==i,:); % bring the ith position for all particles
    if i == 1
        plot(initposX.*10^6,initposY.*10^6,'x','Color','r')
        hold on
        viscircles([stepframe(:,1),stepframe(:,2)], ones(length(stepframe(:,1)),1)*radius,'Color','blue');
        xlim([-lx lx]);
        ylim([-lx lx]);
        hold off
        pause(2)
    else
        if isempty(checkResetBefore) && isempty(checkResetAfter)
            cla; clf;
            viscircles([stepframe(:,1),stepframe(:,2)], ones(length(stepframe(:,1)),1)*radius,'Color','blue');
            xlim([-lx lx]);
            ylim([-lx lx]);
            title(strcat("step number - ",num2str(i)))
        elseif isempty(checkResetAfter)
            cla
            ResetParticleIdx = ResetData(checkResetBefore, 4);
            TrapLocX = ResetData(checkResetBefore, 1);
            TrapLocY = ResetData(checkResetBefore, 2);
            NotResetparticle = setdiff(ParticleNames,ResetParticleIdx);
            plot(TrapLocX,TrapLocY,'x','Color','r')
            hold on
            viscircles([stepframe(NotResetparticle,1),stepframe(NotResetparticle,2)], ones(length(NotResetparticle),1)*radius,'Color','blue');
            viscircles([stepframe(ResetParticleIdx,1),stepframe(ResetParticleIdx,2)], radius,'Color','red');
            
            xlim([-lx lx]);
            ylim([-lx lx]);
            title(strcat("step number - ",num2str(i)))
            subtitle("Reset Event Reset Event",'Color','red')
            hold off
            pause(2)

        else
            cla
            ResetParticleIdx = ResetData(checkResetAfter, 4);
            NotResetparticle = setdiff(ParticleNames,ResetParticleIdx);
            plot(TrapLocX,TrapLocY,'x','Color','r')
            hold on
            viscircles([stepframe(NotResetparticle,1),stepframe(NotResetparticle,2)], ones(length(NotResetparticle),1)*radius,'Color','blue');
            viscircles([stepframe(ResetParticleIdx,1),stepframe(ResetParticleIdx,2)], radius,'Color','red');
            xlim([-lx lx]);
            ylim([-lx lx]);
            title(strcat("step number - ",num2str(i)))
            subtitle("Reset Event Reset Event",'color','red')
            hold off
            pause(2)
        end
    end

        pause(1/fps);
end
%% With Video record
trackfileAll = dir("tr_all*.mat");
fps = 30;
initposX = [1.12854917132807e-06; 2.23527532474401e-06; 1.13321873770091e-06; -1.17441058653142e-06; -2.19776133851434e-06; -1.10358946176809e-06]; % Xtraps locations from previous experiments
initposY = [-1.99263553440715e-06; -1.23133715505655e-07; 1.96279242978705e-06; -2.03044832971602e-06; -1.19446040672767e-08; 2.05004525027041e-06]; % Ytraps locations
trackfileAll = struct2array(load(trackfileAll(1).name)); 
ParticleNames = unique(trackfileAll(:,6));
NumOfSteps = length(unique(trackfileAll(:,5)));
radius = 0.75; % micrometer
lx = 30*radius;
ResetEventsIdx = find(trackfileAll(:,8) == 1);
ResetEventSteps = trackfileAll(ResetEventsIdx, 5);
ResetEventParticle = trackfileAll(ResetEventsIdx, 6);
TrapLocsAtReset = [trackfileAll(ResetEventsIdx, 1) trackfileAll(ResetEventsIdx, 2)];
ResetData = [TrapLocsAtReset, ResetEventSteps, ResetEventParticle];
v = VideoWriter('new.avi');
open(v)
for i = 1:NumOfSteps
    checkResetBefore = find(ResetEventSteps == i+1);
    checkResetAfter = find(ResetEventSteps == i);
    stepframe = trackfileAll(trackfileAll(:,5)==i,:); % bring the ith position for all particles
    if i == 1
        plot(initposX.*10^6,initposY.*10^6,'x','Color','r')
        hold on
        viscircles([stepframe(:,1),stepframe(:,2)], ones(length(stepframe(:,1)),1)*radius,'Color','blue');
        xlim([-lx lx]);
        ylim([-lx lx]);
        hold off
        saveas(figure(1), "A.jpeg")
        A = imread("A.jpeg");
        for k = 1:60 % equivalent to 2sec pause
            writeVideo(v,A)
        end
        %pause(2)
    else
        if isempty(checkResetBefore) && isempty(checkResetAfter)
            cla; clf;
            viscircles([stepframe(:,1),stepframe(:,2)], ones(length(stepframe(:,1)),1)*radius,'Color','blue');
            xlim([-lx lx]);
            ylim([-lx lx]);
            title(strcat("step number - ",num2str(i)))
            saveas(figure(1), "A.jpeg")
            A = imread("A.jpeg");
            writeVideo(v,A)
        elseif isempty(checkResetAfter)
            cla
            ResetParticleIdx = ResetData(checkResetBefore, 4);
            TrapLocX = ResetData(checkResetBefore, 1);
            TrapLocY = ResetData(checkResetBefore, 2);
            NotResetparticle = setdiff(ParticleNames,ResetParticleIdx);
            plot(TrapLocX,TrapLocY,'x','Color','r')
            hold on
            viscircles([stepframe(NotResetparticle,1),stepframe(NotResetparticle,2)], ones(length(NotResetparticle),1)*radius,'Color','blue');
            viscircles([stepframe(ResetParticleIdx,1),stepframe(ResetParticleIdx,2)], radius,'Color','red');
            
            xlim([-lx lx]);
            ylim([-lx lx]);
            title(strcat("step number - ",num2str(i)))
            subtitle("Reset Event Reset Event",'Color','red')
            hold off
            saveas(figure(1), "A.jpeg")
            A = imread("A.jpeg");
            for k = 1:60 % equivalent to 2sec pause
                 writeVideo(v,A)
            end
          %  pause(2)

        else
            cla
            ResetParticleIdx = ResetData(checkResetAfter, 4);
            NotResetparticle = setdiff(ParticleNames,ResetParticleIdx);
            plot(TrapLocX,TrapLocY,'x','Color','r')
            hold on
            viscircles([stepframe(NotResetparticle,1),stepframe(NotResetparticle,2)], ones(length(NotResetparticle),1)*radius,'Color','blue');
            viscircles([stepframe(ResetParticleIdx,1),stepframe(ResetParticleIdx,2)], radius,'Color','red');
            xlim([-lx lx]);
            ylim([-lx lx]);
            title(strcat("step number - ",num2str(i)))
            subtitle("Reset Event Reset Event",'color','red')
            hold off
            saveas(figure(1), "A.jpeg")
            A = imread("A.jpeg");
            for k = 1:60 % equivalent to 2sec pause
                writeVideo(v,A)
            end
            
            %pause(2)
        end
    end

        pause(1/fps);
end
        close(v)
    
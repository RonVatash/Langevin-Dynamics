function printCurrStepActivePassive(cfg, currStepData, addedData)
    cla     
    if cfg.useWalls
        plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(1)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(2)],'-');
        hold on
        plot([currStepData.wallPositionsX(2),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(2)],'-');
        plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(1),currStepData.wallPositionsY(1)],'-');
        plot([currStepData.wallPositionsX(1),currStepData.wallPositionsX(2)],[currStepData.wallPositionsY(2),currStepData.wallPositionsY(2)],'-');
    end
    if isempty(cfg.Nfast)
       xlim(cfg.xlimits)
       ylim(cfg.ylimits)
        viscircles([currStepData.particlePositions(:,1),...
                    currStepData.particlePositions(:,2)],...
            ones(cfg.numOfParticles,1).*cfg.R(1),'Color','b');
        hold off
        pause(0.001);
    else
       cla
       xlim(cfg.xlimits)
       ylim(cfg.ylimits)
       hold on
       title("Active vs Passive diffusion")
       if isempty(cfg.Rfast)
           RF = cfg.R(1);
       else
           RF = cfg.Rfast(1);
       end
       if cfg.useTraps && currStepData.resetParticle(1)
           cla 
           xlim(cfg.xlimits)
           ylim(cfg.ylimits)
           hold on
           title(strcat("Reset event number__", num2str(currStepData.ResetNum)))
           ResetParticlIdx = currStepData.resetParticle;
           plot(currStepData.CurtrapPositions(: ,1),...
                    currStepData.CurtrapPositions(: ,2),...
            '.', 'Color','g', MarkerSize=15);
            Trapidx = currStepData.trapidx;
           plot(currStepData.trapPositions(Trapidx ,1),...
                    currStepData.trapPositions(Trapidx ,2),...
            '.', 'Color','g', MarkerSize=15);
                viscircles([currStepData.particlePositions(ResetParticlIdx ,1),...
                    currStepData.particlePositions(ResetParticlIdx ,2)],...
            ones(length(ResetParticlIdx) ,1).*RF,'Color','m');
                
                FastParticleIdx = [1:1:cfg.Nfast];
                OtherParticle = setdiff(FastParticleIdx, ResetParticlIdx);
                        viscircles([currStepData.particlePositions(OtherParticle ,1),...
                    currStepData.particlePositions(OtherParticle ,2)],...
            ones(length(OtherParticle),1).*RF,'Color','r');
        hold on
             viscircles([currStepData.particlePositions(cfg.Nfast+1:end ,1),...
                    currStepData.particlePositions(cfg.Nfast + 1:end ,2)],...
            ones(cfg.numOfParticles-cfg.Nfast,1).*cfg.R(1),'Color','b');
        hold off
        pause(0.001);

        % plot(currStepData.trapPositions(currStepData.trapidx ,1),...
        %             currStepData.trapPositions(currStepData.trapidx ,2),...
        %     '.', 'Color','r', MarkerSize=10);
       else
        viscircles([currStepData.particlePositions(1:cfg.Nfast ,1),...
                    currStepData.particlePositions(1:cfg.Nfast ,2)],...
            ones(cfg.Nfast,1).*RF,'Color','r');
        hold on
             viscircles([currStepData.particlePositions(cfg.Nfast+1:end ,1),...
                    currStepData.particlePositions(cfg.Nfast + 1:end ,2)],...
            ones(cfg.numOfParticles-cfg.Nfast,1).*cfg.R(1),'Color','b');
        hold off
        pause(0.001);
       end

    end
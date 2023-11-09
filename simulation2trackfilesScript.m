clear; clc;
TopLevel = pwd;
Folders = dir(TopLevel);
Folders = Folders(3:end);
tr_all_mks_dist = [];
for i = 2230:length(Folders)
    Xpos6paricles = load(strcat(Folders(i).name,"\pos_x.csv"));
    Ypos6paricles = load(strcat(Folders(i).name,"\pos_y.csv"));
    NumOfParticles = length(Xpos6paricles(1,:));
    FrameNum = 1:length(Xpos6paricles);
    FrameNum = FrameNum';
    trackfile = [];
    Meter2Micro = 10^6;
    for j = 1:NumOfParticles
        TrParticleJ = zeros(length(FrameNum),7);
        TrParticleJ(:,1) = Xpos6paricles(:,j)*Meter2Micro;
        TrParticleJ(:,2) = Ypos6paricles(:,j)*Meter2Micro;
        TrParticleJ(:,5) = FrameNum;
        TrParticleJ(:,6) = j;
        XposSq = TrParticleJ(:,1).^2; %The position in x squared
        YposSq = TrParticleJ(:,2).^2; %The position in y squared
        TrParticleJ(:,7) = sqrt(XposSq + YposSq); %The Distance from the origin
        trackfile = [trackfile; TrParticleJ];
    end
    if ~ mod(i,100)
        i
    end
    save(strcat('tr',num2str(i,'%06d'),'.mat'),"trackfile");
    tr_all_mks_dist = [tr_all_mks_dist; trackfile];
end

save(strcat('tr_all_dist_',num2str(length(Folders),'%06d'),'_movies.mat'),"tr_all_mks_dist")
%% With reset
clear; clc;
LastFrame = [];
Success = [];
NumOfFramesTot = 0;
FrameNum = [];
TopLevel = pwd;
Folders = dir(TopLevel);
Folders = Folders(3:end);
tr_all_mks_dist = [];
for i = 1:length(Folders)
    Xpos6paricles = load(strcat(Folders(i).name,"\pos_x.csv"));
    Ypos6paricles = load(strcat(Folders(i).name,"\pos_y.csv"));
    if i ~= 1
       ResetData = struct2array(load(strcat(Folders(i-1).name,"\RestData_",num2str(i,'%06d'),".mat")));
       ResetParticleidx = find(round(Xpos6paricles(1,:),8) == round(ResetData(1),8));
    end
  
    NumOfParticles = length(Xpos6paricles(1,:));
    FrameNum = 1+NumOfFramesTot:NumOfFramesTot+length(Xpos6paricles(:,1));
    FrameNum = FrameNum';
    LastFrame = [LastFrame; length(Xpos6paricles(:,1))];
    trackfile = [];
    Meter2Micro = 10^6;
    for j = 1:NumOfParticles
        TrParticleJ = zeros(length(FrameNum),8);
        TrParticleJ(:,1) = Xpos6paricles(:,j)*Meter2Micro;
        TrParticleJ(:,2) = Ypos6paricles(:,j)*Meter2Micro;
        TrParticleJ(:,5) = FrameNum;
        TrParticleJ(:,6) = j;
        XposSq = TrParticleJ(:,1).^2; %The position in x squared
        YposSq = TrParticleJ(:,2).^2; %The position in y squared
        TrParticleJ(:,7) = sqrt(XposSq + YposSq); %The Distance from the origin
        if i ~= 1
            if j == ResetParticleidx
                TrParticleJ(1,8) = ResetData(4);
                Success = [Success, ResetData(4)];
            end
        end
        trackfile = [trackfile; TrParticleJ];
    end
    if ~ mod(i,100)
        i
    end
    save(strcat('tr',num2str(i,'%06d'),'.mat'),"trackfile");
    NumOfFramesTot = sum(LastFrame);
    tr_all_mks_dist = [tr_all_mks_dist; trackfile];
end
tr_all_mks_dist = sortrows(tr_all_mks_dist,6,"ascend");
NumOfResetEvents = length(Folders)-1; % the number of attemps
NumOfSuccesses = sum(Success);
%ResetInformation = [NumOfResetEvents, NumOfSuccesses, alpha];
alpha = NumOfSuccesses/NumOfResetEvents; % the percentage of success
save(strcat('tr_all_dist_',num2str(length(Folders),'%06d'),'_movies.mat'),"tr_all_mks_dist")
save("ResetInformation.mat","alpha","NumOfResetEvents","NumOfSuccesses");
%% With reset from many directories
clear; clc;
LastFrame = [];
Success = [];
NumOfFramesTot = 0;
FrameNum = [];
TopLevel = pwd;
MainFolders = dir(TopLevel);
MainFolders = MainFolders(3:end);
for j = 1:length(MainFolders)
    Folders = dir(MainFolders(j).name);
    Folders = Folders(3:end);
    tr_all_mks_dist = [];
    for i = 1:length(Folders)
        Xpos6paricles = load(strcat(Folders(i).folder,'\',Folders(i).name,"\pos_x.csv"));
        Ypos6paricles = load(strcat(Folders(i).folder,'\',Folders(i).name,"\pos_y.csv"));
        if i ~= 1
           ResetData = struct2array(load(strcat(Folders(i).folder,'\',Folders(i-1).name,"\RestData_",num2str(i,'%06d'),".mat")));
           ResetParticleidx = find(round(Xpos6paricles(1,:),8) == round(ResetData(1),8));
        end
      
        NumOfParticles = length(Xpos6paricles(1,:));
        FrameNum = 1+NumOfFramesTot:NumOfFramesTot+length(Xpos6paricles);
        FrameNum = FrameNum';
        LastFrame = [LastFrame; length(FrameNum)];
        trackfile = [];
        Meter2Micro = 10^6;
        for j = 1:NumOfParticles
            TrParticleJ = zeros(length(FrameNum),8);
            TrParticleJ(:,1) = Xpos6paricles(:,j)*Meter2Micro;
            TrParticleJ(:,2) = Ypos6paricles(:,j)*Meter2Micro;
            TrParticleJ(:,5) = FrameNum;
            TrParticleJ(:,6) = j;
            XposSq = TrParticleJ(:,1).^2; %The position in x squared
            YposSq = TrParticleJ(:,2).^2; %The position in y squared
            TrParticleJ(:,7) = sqrt(XposSq + YposSq); %The Distance from the origin
            if i ~= 1
                if j == ResetParticleidx
                    TrParticleJ(1,8) = ResetData(4);
                    Success = [Success, ResetData(4)];
                end
            end
            trackfile = [trackfile; TrParticleJ];
        end
        if ~ mod(i,100)
            i
        end
        save(strcat(Folders(i).folder,'\tr',num2str(i,'%06d'),'.mat'),"trackfile");
        NumOfFramesTot = sum(LastFrame);
        tr_all_mks_dist = [tr_all_mks_dist; trackfile];
    end
    tr_all_mks_dist = sortrows(tr_all_mks_dist,6,"ascend");
    NumOfResetEvents = length(Folders)-1; % the number of attemps
    NumOfSuccesses = sum(Success);
    %ResetInformation = [NumOfResetEvents, NumOfSuccesses, alpha];
    alpha = NumOfSuccesses/NumOfResetEvents; % the percentage of success
    save(strcat(Folders(i).folder,'\tr_all_dist_',num2str(length(Folders),'%06d'),'_movies.mat'),"tr_all_mks_dist")
    save(strcat(Folders(i).folder,"\ResetInformation.mat"),"alpha","NumOfResetEvents","NumOfSuccesses");
end
%% tracks to trak
tracks = dir('tr*.mat');
tr_all_mks_dist = [];
for i = 2824:2900%length(tracks)
    tr = struct2array(load(tracks(i).name));
    tr_all_mks_dist = [tr_all_mks_dist; tr];
end
tr_all_mks_dist_part =  tr_all_mks_dist;    
save(strcat('tr_all_mks_dist_part_',num2str(21),'.mat'),"tr_all_mks_dist_part");

%save(strcat('tr_all_mks_dist_',num2str(length(tracks),'%06d'),'_movies.mat'),"tr_all_mks_dist");
%% track 2 tracks (if track is too big) 
tillhere = 2824;        
numOftracks = 20;
lngtOftr = (length(tr_all_mks_dist(:,1)))/numOftracks;
start = 1;
for i = 1:numOftracks
    stop = floor(i*lngtOftr);
    tr_all_mks_dist_part = tr_all_mks_dist(start:stop,:);
    start = stop+1;
    save(strcat('tr_all_mks_dist_part_',num2str(i,'%03d'),'.mat'),"tr_all_mks_dist_part");
end
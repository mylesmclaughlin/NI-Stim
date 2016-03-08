function A = NIavgTempDynamics(fileName,dataPath)

if nargin < 2
    dataPath = 'C:\Users\Public\Data\NIData\';
end
binFile = [dataPath fileName];

[p,f,e] = fileparts(binFile);
N = load([p '\' f '.nis'],'-mat');
load([p '\' f '.NISmacro'],'-mat');
N = N.D;
A.repDur = N.stim.burstdur;
A.burstrepperiod = N.stim.burstrepperiod;

% load bin file
extraTime = 100;
sampWin = [-extraTime A.burstrepperiod]; %A.repDur+extraTime]; % ms
D = binread(binFile);
sInd = strfind(D.header,'Fs=');
eInd = strfind(D.header,'Filter=');
fs = str2num(D.header(sInd+3:eInd-3));
sampVec = [round(sampWin(1)/1e3*fs):round(sampWin(2)/1e3*fs)];
trigChInd = 3;
trigLevel = 0.5;

part = 1;
allData = [];
nTrig = 0;
keepData = [];
startTime = 0;
while part
    disp(['Processing PART ' num2str(part) '...'])
    D = binread(binFile,part);
    if ~isempty(keepData)
        D.data = [keepData; D.data];
        keepData = [];
    end
    if D.data(1,trigChInd)>trigLevel
        D.data(1,trigChInd) = 0;
    end
    trigInd = find(D.data(1:end-1,trigChInd)<trigLevel & D.data(2:end,trigChInd)>trigLevel);
    trigTimes = trigInd/fs + startTime;
    for n = 1:length(trigInd)
        if trigInd(n)+sampVec(1)>0 & trigInd(n)+sampVec(end)<length(D.data)
            nTrig = nTrig+1;
            allData(nTrig,:,:) = D.data(trigInd(n)+sampVec,:);
            timeStamps(nTrig) = trigTimes(n);
            %disp(['...found trigger number ' num2str(nTrig)])
        elseif trigInd(n)+sampVec(1)<=0
            nTrig = nTrig+1;
            posInd = find(trigInd(n)+sampVec>0);
            allData(nTrig,:,:) = zeros(size(allData(nTrig-1,:,:)));
            allData(nTrig,posInd,:) = D.data(trigInd(n)+sampVec(posInd),:);
            timeStamps(nTrig) = trigTimes(n);
            disp('WARNING: Trigger error in data processing. Some zero values')
            %return
        elseif trigInd(n)+sampVec(end)>=length(D.data)
            % keep bit of data for next run
            keepData = D.data(trigInd(n)+(round(sampWin(1)/1e3*fs))*2:end,:);
            break
        end
    end
    if D.part == 0
        part = D.part;
    else
        part = part+1;
    end
    
    oldData = D.data;
    if ~isempty(keepData)
        endTime = (length(D.data)-length(keepData))/fs;
    else
        endTime = length(D.data)/fs;
    end
    startTime = endTime;
end


A.fs = fs;
A.nTrig = nTrig;
A.allData = allData;
A.timeStamps = timeStamps;
A.tvec = sampVec/A.fs;

% sort into prestim, stim and poststim
nSeqOn = P.preset.subThresDurOn/(P.stim.burstrepperiod/1000);
nSeqOff = P.preset.subThresDurOff/(P.stim.burstrepperiod/1000);
nSeq = nSeqOff+nSeqOn;
nReps = P.stim.numberreps;
A.preStimInd = [];
A.stimInd = [];
A.postStimInd = [];
if P.preset.subThresOffFirst == 0
    A.preStimInd = [];
    A.stimInd = [];
    A.postStimInd = [];
elseif P.preset.subThresOffFirst == 1
    A.preStimInd = [1:nSeqOff];
    A.stimInd = [nSeqOff+1:nSeq];
    A.postStimInd = [nSeq+1:nSeq+nSeqOff];
    for n = 2:nReps
        A.stimInd = [A.stimInd; A.stimInd(end)+nSeqOff+1:A.stimInd(end)+nSeq];
    end
    for n = 2:nReps-1
        A.postStimInd = [A.postStimInd; A.postStimInd(end)+nSeqOn+1:A.postStimInd(end)+nSeq];
    end
    A.preStimNReps = 1;
    A.stimNReps = nReps;
    A.postStimNReps = nReps-1;
end
A.nSeqOn = nSeqOn;
A.nSeqOff = nSeqOff;

A.stimDataDimensionNames = {'Time','Rep','Data','Channel'};
for n = 1:A.stimNReps
    A.stimData(:,n,:,:) = A.allData(A.stimInd(n,:),:,:);
end
for n = 1:A.postStimNReps
    A.postStimData(:,n,:,:) = A.allData(A.postStimInd(n,:),:,:);
end
A.preStimData(:,1,:,:) = A.allData(A.preStimInd,:,:);

% convert to displacement
accelInd = [4 5 6];
bandpass = [3 500];
filterorder = 2;
[filterb filtera] = butter(filterorder,bandpass/(A.fs/2),'bandpass');
zInd = find(A.tvec<0);
stimStartInd = zInd(end)+1;
trigArtInd = find(A.tvec<15e-3);
stimStartInd = trigArtInd(end)+1;

A.stimDisp = zeros(size(A.stimData));

for n = 1:A.stimNReps
    for m = 1:A.nSeqOn
        data = squeeze(A.stimData(m,n,stimStartInd:end,:));
        calib = mean(data(zInd));
        data = data - calib;
        data = filtfilt(filterb,filtera,data);
        data = cumtrapz(data);
        data = filtfilt(filterb,filtera,data);
        A.stimDisp(m,n,stimStartInd:end,:) = cumtrapz(data);
    end
end
for n = 1:A.postStimNReps
    for m = 1:A.nSeqOff
        data = squeeze(A.postStimData(m,n,:,:));
        calib = mean(data(zInd));
        data = data - calib;
        data = filtfilt(filterb,filtera,data);
        data = cumtrapz(data);
        data = filtfilt(filterb,filtera,data);
        A.postStimDisp(m,n,:,:) = cumtrapz(data);
    end
end
for m = 1:A.nSeqOff
    data = squeeze(A.preStimData(m,1,:,:));
    calib = mean(data(zInd));
    data = data - calib;
    data = filtfilt(filterb,filtera,data);
    data = cumtrapz(data);
    data = filtfilt(filterb,filtera,data);
    A.preStimDisp(m,1,:,:) = cumtrapz(data);
end



% take average
A.stimDataAvgDimensionNames = {'Time','Data','Channel'};

A.stimDataAvg = squeeze(mean(A.stimData,2));
A.postStimDataAvg = squeeze(mean(A.postStimData,2));
A.preStimDataAvg = squeeze(mean(A.preStimData,2));


A.stimDispAvg = squeeze(mean(A.stimDisp,2));
A.postStimDispAvg = squeeze(mean(A.postStimDisp,2));
if ndims(A.postStimDispAvg)<3
    a(1,:,:) = A.postStimDispAvg;
    A.postStimDispAvg = a;
end
A.preStimDispAvg = squeeze(mean(A.preStimDisp,2));
if ndims(A.preStimDispAvg)<3
    a(1,:,:) = A.preStimDispAvg;
    A.preStimDispAvg = a;
end


A.stimDispSTD = squeeze(std(A.stimDisp,0,2));
A.postStimDispSTD = squeeze(std(A.postStimDisp,0,2));
if ndims(A.postStimDispSTD)<3
    a(1,:,:) = A.postStimDispSTD;
    A.postStimDispSTD = a;
end
%A.preStimDispSTD = squeeze(std(A.preStimDisp,0,2));

% find peaks
A.stimRespAmpDimensionNames = {'Time','AccelAxis'};
accelChInd = [4 5 6];
for m = 1:A.nSeqOn
    for i = 1:length(accelChInd)
        [val,ind] = max(abs(A.stimDispAvg(m,:,accelChInd(i))));
        A.stimRespAmp(m,i) = val;
        A.stimRespSTD(m,i) = A.stimDispSTD(m,ind,accelChInd(i));
    end
end
for m = 1:A.nSeqOff
    for i = 1:length(accelChInd)
        [val,ind] = max(abs(A.postStimDispAvg(m,:,accelChInd(i))));
        A.postStimRespAmp(m,i) = val;
        A.postStimRespSTD(m,i) = A.postStimDispSTD(m,ind,accelChInd(i));
    end
end
for m = 1:A.nSeqOff
    for i = 1:length(accelChInd)
        [val,ind] = max(abs(A.preStimDispAvg(m,:,accelChInd(i))));
        A.preStimRespAmp(m,i) = val;
        %A.preStimRespSTD(m,i) = A.preStimDispSTD(m,ind,accelChInd(i));
    end
end

% -- check plot --
% figure; hold on
% for n = 1:A.nTrig
%     plot(A.tvec+A.timeStamps(n),A.allData(n,:,3))
%     plot(A.tvec+A.timeStamps(n),A.allData(n,:,1),'r')
% end

accelCh = 4;
for i = 1:3
    if i == 1
        accelCh = 4;
        figName = 'X-axis';
    elseif  i == 2
        accelCh = 5;
        figName = 'Y-axis';
    elseif i == 3
        accelCh = 6;
        figName = 'Z-axis';
    end
    
    figure('name',figName)
    a = A.preStimDispAvg(:,:,4:6);
    b = A.stimDispAvg(:,:,4:6);
    c = A.postStimDispAvg(:,:,4:6);
    
    yl(1) = min([a(:); b(:); c(:)]);
    yl(2) = max([a(:); b(:); c(:)]);
    
    subplot(1,3,1)
    for n = 1:A.nSeqOff
        plot(A.tvec+A.timeStamps(n),A.preStimDispAvg(n,:,accelCh))
        hold on
    end
    ylim(yl)
    title('Pre stim')
    
    subplot(1,3,2)
    for n = 1:A.nSeqOn
        plot(A.tvec+A.timeStamps(n),A.stimDispAvg(n,:,accelCh),'r')
        hold on
    end
    ylim(yl)
    title('Stim')
    
    subplot(1,3,3)
    for n = 1:A.nSeqOff
        plot(A.tvec+A.timeStamps(n),A.postStimDispAvg(n,:,accelCh),'g')
        hold on
    end
    ylim(yl)
    title('Post stim')
    
    
    figure('name',figName)
    yl(1) = min([A.preStimRespAmp(:); A.stimRespAmp(:)-A.stimRespSTD(:); A.postStimRespAmp(:)-A.postStimRespSTD(:)]);
    yl(2) = max([A.preStimRespAmp(:); A.stimRespAmp(:)+A.stimRespSTD(:); A.postStimRespAmp(:)+A.postStimRespSTD(:)]);
    
    subplot(1,3,1)
    plot(A.timeStamps(1:A.nSeqOff),A.preStimRespAmp(:,accelCh-3),'b.-')
    ylim(yl)
    title('Pre stim')
    
    subplot(1,3,2)
    %plot(A.timeStamps(1:A.nSeqOn),A.stimRespAmp(:,accelCh-3),'b.-')
    errorbar(A.timeStamps(1:A.nSeqOn),A.stimRespAmp(:,accelCh-3),A.stimRespSTD(:,accelCh-3),'r')
    ylim(yl)
    title('Stim')
    
    subplot(1,3,3)
    %plot(A.timeStamps(1:A.nSeqOff),A.postStimRespAmp(:,accelCh-3),'b.-')
    errorbar(A.timeStamps(1:A.nSeqOff),A.postStimRespAmp(:,accelCh-3),A.postStimRespSTD(:,accelCh-3),'g')
    ylim(yl)
    title('Post stim')
end

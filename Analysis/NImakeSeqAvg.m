function A = NImakeSeqAvg(binFile,dataPath,detrendChIndex)

if nargin<2
    dataPath = 'C:\Users\Public\Data\NIData\';
end
binFile = [dataPath binFile];

if nargin<3
    detrendChIndex = [];
end
%load .nis file
[p,f,e] = fileparts(binFile);
try
    N = load([p '\' f '.niso'],'-mat');
catch
    N = load([p '\' f '.nis'],'-mat');
end

N = N.D;
A.repDur = N.stim.burstdur;
A.burstrepperiod = N.stim.burstrepperiod;

%load bin file
extraTime = 100;
sampWin = [-extraTime A.burstrepperiod]; %A.repDur+extraTime]; % ms
D = binread(binFile);
sInd = strfind(D.header,'Fs=');
eInd = strfind(D.header,'Filter=');
fs = str2num(D.header(sInd+3:eInd-3));
sampVec = [round(sampWin(1)/1e3*fs):round(sampWin(2)/1e3*fs)];
trigChInd = 3;
trigLevel = 0.5;

%detrend data
for n = 1:length(detrendChIndex)
    D.data(:,detrendChIndex(n)) = detrend(D.data(:,detrendChIndex(n)));
end



part = 1;
allData = [];
nTrig = 0;
keepData = [];
while part
    disp(['Processing PART ' num2str(part) '...'])
    D = binread(binFile,part);
    if ~isempty(keepData)
        D.data = [keepData; D.data];
        keepData = [];
    end
    trigInd = find(D.data(1:end-1,trigChInd)<trigLevel & D.data(2:end,trigChInd)>trigLevel);
    for n = 1:length(trigInd)
        if trigInd(n)+sampVec(1)>0 & trigInd(n)+sampVec(end)<length(D.data) 
            nTrig = nTrig+1;
            allData(nTrig,:,:) = D.data(trigInd(n)+sampVec,:);
            %disp(['...found trigger number ' num2str(nTrig)])
        elseif trigInd(n)+sampVec(1)<=0
            disp('WARNING: Trigger error in data processing')
            return
        elseif trigInd(n)+sampVec(end)>=length(D.data)
%            keep bit of data for next run
            keepData = D.data(trigInd(n)+(round(sampWin(1)/1e3*fs))*2:end,:);
            break
        end
    end   
    if D.part == 0
        part = D.part;
    else
        part = part+1;
    end
end


A.fs = fs;
A.nTrig = nTrig;
A.allData = allData;
A.tvec = sampVec/A.fs;


if isfield(N.sequence,'seq')
    A.seq = N.sequence.seq; %[1:10];
    A.seqIndex = N.sequence.seqIndex; %repmat(seq,1,4);
    %A.seqIndex =  repmat(A.seq,1,floor(A.nTrig/length(A.seq)));
else
    A.seq = 1:10;
    A.seqIndex = repmat(A.seq,1,4);
end

% convert to displacement
accelInd = [4 5 6];
bandpass = [3 500];
filterorder = 2;
[filterb filtera] = butter(filterorder,bandpass/(A.fs/2),'bandpass');
zInd = find(A.tvec<0);
stimStartInd = zInd(end)+1;
trigArtInd = find(A.tvec<15e-3);
stimStartInd = trigArtInd(end)+1;

A.allDisp = zeros(size(A.allData));
A.nCh = size(A.allData,3);

for n = 1:length(A.seqIndex)
    for m = 1:A.nCh
        data = squeeze(A.allData(n,stimStartInd:end,m));
        calib = mean(data(zInd));
        data = data - calib;
        data = filtfilt(filterb,filtera,data);
        data = cumtrapz(data);
        data = filtfilt(filterb,filtera,data);
        A.allDisp(n,stimStartInd:end,m) = cumtrapz(data);
    end
end



% calculate average and std
for n = 1:length(A.seq)
    eval(['A.allData' num2str(n) ' = A.allData(find(A.seqIndex==n),:,:);']);
    eval(['A.allDisp' num2str(n) ' = A.allDisp(find(A.seqIndex==n),:,:);']);
    if length(find(A.seqIndex==n)) == 1
        eval(['A.avgData(' num2str(n) ',:,:) = squeeze(A.allData(find(A.seqIndex==n),:,:));']);
        eval(['A.stdData(' num2str(n) ',:,:) = squeeze(A.allData(find(A.seqIndex==n),:,:));']);
        eval(['A.avgDisp(' num2str(n) ',:,:) = squeeze(A.allDisp(find(A.seqIndex==n),:,:));']);
        eval(['A.stdDisp(' num2str(n) ',:,:) = squeeze(A.allDisp(find(A.seqIndex==n),:,:));']);
    else
        eval(['A.avgData(' num2str(n) ',:,:) = squeeze(mean(A.allData(find(A.seqIndex==n),:,:)));']);
        eval(['A.stdData(' num2str(n) ',:,:) = squeeze(std(A.allData(find(A.seqIndex==n),:,:)));']);
        eval(['A.avgDisp(' num2str(n) ',:,:) = squeeze(mean(A.allDisp(find(A.seqIndex==n),:,:)));']);
        eval(['A.stdDisp(' num2str(n) ',:,:) = squeeze(std(A.allDisp(find(A.seqIndex==n),:,:)));']);
        
    end
end

A.nSeq = length(A.seq);
A.seqparametername = N.sequence.parametername;
A.seqparametervalues = N.sequence.parametervalues; 

if strcmp(A.seqparametername,'frequency')
    A.frequency = A.seqparametervalues;
    A.amplitude = N.stim.amplitude;
elseif strcmp(A.seqparametername,'amplitude')
    A.frequency = N.stim.frequency;
    A.amplitude = A.seqparametervalues;
elseif strcmp(A.seqparametername,'ampmoddepth')
    A.frequency = N.stim.frequency;
    A.amplitude = N.stim.amplitude;
    A.ampmod = A.seqparametervalues;
    A.amfreq = N.stim.ampmodfreq;
end



function A = NImakeAvg(binFile,detrendChIndex)

if nargin<2
    detrendChIndex = [];
end
% load .nis file
[p,f,e] = fileparts(binFile);
try
    N = load([p '\' f '.niso'],'-mat');
catch
    N = load([p '\' f '.nis'],'-mat');
end

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

% detrend data
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
end


A.fs = fs;
A.nTrig = nTrig;
A.allData1 = allData;
A.avgData = mean(A.allData1);
A.stdData = std(A.allData1);
A.tvec = sampVec/A.fs;
A.nSeq = 1;

A.frequency = N.stim.frequency;
A.amplitude = N.stim.amplitude;
A.amfreq = N.stim.ampmodfreq;
A.seqparametervalues = A.amplitude;
A.seqparametername = 'Amplitude';
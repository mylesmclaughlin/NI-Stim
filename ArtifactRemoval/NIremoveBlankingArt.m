function D = NIremoveBlankingArt(binFile)


dataChInd = [1 2];
trigChInd = 3;
blankWin = [-10 50]; % in samples
downFac = 5;
downData = [];
downStimOn = [];

[p,f,e] = fileparts(binFile);
load([p '\' f '.nis'],'-mat')
if D.stim.continuous
    buffersize = D.stim.buffersize;
else
    buffersize = round(D.stim.burstdur*D.ni.rate);
end
stimFreq = D.stim.frequency;
D = [];
part = 1;
while part~=0
    D = binread(binFile,part);
    if part == 1
        sInd = strfind(D.header,'Fs=');
        eInd = strfind(D.header,'Filter=');
        fs = str2num(D.header(sInd+3:eInd-3));
        downFs = fs/downFac;
    end    
    
    [data,stimOn] = removeart(D.data,dataChInd,trigChInd,blankWin,fs,buffersize,stimFreq);
    locD = [];
    for n = 1:length(dataChInd)
        locD(:,n) = decimate(data(:,n),downFac);
    end
    downData = [downData; locD];
    downStimOn = [downStimOn; stimOn(1:downFac:end)];
    
    figure
    plot(D.data(:,1),'r')
    hold on
    plot(data(:,1),'b')
    plot(stimOn,'m')
    
    part = D.part;
    if part~=0
        part = part+1;
    end
end
D = [];
D.fs = downFs;
D.data = downData;
D.stimOn = downStimOn;
save([p '\' f '.mat'],'D')


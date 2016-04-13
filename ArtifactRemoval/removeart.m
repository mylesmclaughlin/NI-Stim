function [data,stimOn,blankInd,stimPulse] = removeart(data,dataChInd,trigChInd,blankWin,fs,buffersize,stimFreq);

if nargin>4 % make stimPulse from trig channel and stim freq
    [stimPulse,stimOn] = makeStimPulse(data(:,trigChInd),fs,buffersize,stimFreq);
else % use recorded stimulus voltage
    stimChInd = trigChInd;
    stimPulse = data(:,stimChInd); 
end
    
blankInd = getPulseInd(stimPulse,blankWin);
nSampBlank = size(blankInd,2);
data = data(:,dataChInd);

for n = 1:length(blankInd)
    for m = 1:length(dataChInd)
        if blankInd(n,end)<length(data)
            data(blankInd(n,:),m) = linspace(data(blankInd(n,1),m),data(blankInd(n,end),m),nSampBlank);
        end
    end
end

%--------------------------------------------------------------------------
function blankInd = getPulseInd(stimPulse,blankWin);
stimPulse = stimPulse - mean(stimPulse);
trigLevel = 1.5*std(stimPulse);
trigInd = find(stimPulse(1:end-1)<trigLevel & stimPulse(2:end)>=trigLevel);
blankInd = blankWin(1):blankWin(2);
blankInd = repmat(blankInd,length(trigInd),1) + repmat(trigInd,1,length(blankInd));

% figure
% plot(stimPulse)
% hold on
% plot(trigInd,stimPulse(trigInd),'r.')

%--------------------------------------------------------------------------
function [stimPulse,stimOn] = makeStimPulse(data,fs,buffersize,stimFreq);

pulsePeriod = round(fs/stimFreq);
pulseVec = zeros(1,buffersize);
pulseVec(1:pulsePeriod:end) = 1;
sampVec = [0:buffersize-1];
    
trigLevel = 0.5;
trigInd = find(data(1:end-1)<trigLevel & data(2:end)>=trigLevel);

stimPulse = zeros(size(data));
stimOn = zeros(size(data));

for n = 1:length(trigInd)
    stimPulse(trigInd(n)+sampVec) = pulseVec;
    stimOn(trigInd(n)+sampVec) = 1;
end


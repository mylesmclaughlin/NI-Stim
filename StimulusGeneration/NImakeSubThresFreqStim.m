function A = NImakeSubThresFreqStim

% Connect with NIStim
global S NI
A.fs = NI.Rate;

%----  Set local params ---
numberreps = 4;

% super thres probe pulse
supthresAmp = 2.8;
freq = 300;
burstdur = 100;
burstrepperiod = 2000;
phase1pulsewidth = 200;
phase2pulsewidth = 200;

% sub thres params
subthresAmp = 0.5;
subthresFreq = [1 2 4 8 16 32 64 128 256]; %[0.25 0.5 0.75 1 2 3 4 6 8]; %
nseq = length(subthresFreq);
phase = 0.25;
subDurOn = 2; %12;
subDurOff = 1; %8;
randomize = 0;
nprobe = floor(subDurOn/(burstrepperiod/1000));

% Make stim
supthresAmp = supthresAmp/S.current.onevoltequalsXmilliamps;
subthresAmp = subthresAmp/S.current.onevoltequalsXmilliamps;
pulse = makePulse(phase1pulsewidth,phase2pulsewidth,supthresAmp,freq,A.fs);
pulseTrain = makePulseTrain(pulse,burstdur,A.fs);

for n = 1:nseq
    data = makeBaseStim(subthresFreq(n),subthresAmp,subDurOn,subDurOff,A.fs);
    [data,trigger] = addPulseStim(data,pulseTrain,burstrepperiod,subthresFreq(n),phase,nprobe,A.fs);
    A.stim.data(n,:,1) = trigger';
    A.stim.data(n,:,2) = data';
end

% Make struct
A.stim.numberreps = numberreps;
A.stim.randomizesequence = randomize;
A.sequence.on = 1;
A.sequence.nseq = nseq;
A.sequence.parametername = 'subThresFreq';
A.sequence.parametervalues = subthresFreq;
A.trigger.on = 0;

% save data
pathname = 'C:\Users\Public\MATLAB\NI-Stim\Stimuli\';
filename = ['SubAmp=' num2str(subthresAmp) '_SupAmp=' num2str(supthresAmp) '_subFreq=' num2str(subthresFreq) '.mat'];
disp(['Saving ' pathname filename])
save([pathname filename],'A')

%--------------------------------------------------------------------------
function [data,trigger] = addPulseStim(data,pulseTrain,burstrepperiod,subthresFreq,phase,nprobes,fs);

triggerdur = 10; % ms
triggeramp = 1;

burstdelay = (1/subthresFreq)*phase;
delaysamps = round(fs*burstdelay);
sampsperburst = round(fs*(burstrepperiod/1000));
pulsePlaceSamps = 0:sampsperburst:nprobes*sampsperburst;
pulsePlaceSamps = pulsePlaceSamps + delaysamps;
trigger = zeros(length(data),1);
nTrigSamps = round(triggerdur*1e-3*fs);
 
for n = 1:nprobes
    pulseInd = pulsePlaceSamps(n)+1:pulsePlaceSamps(n)+length(pulseTrain);
    data(pulseInd) = pulseTrain;
    
    trigInd = pulsePlaceSamps(n)+1:pulsePlaceSamps(n)+nTrigSamps;
    trigger(trigInd) = triggeramp;
end
    
        
%--------------------------------------------------------------------------     
function [data,tvec] = makeBaseStim(subthresFreq,subthresAmp,subDurOn,subDurOff,fs);

subDurOn = 12;
subDurOff = 8;

sampperburst = round(fs*subDurOn);
tvec = [1/fs:1/fs:sampperburst/fs];
data = subthresAmp*sin(2*pi*subthresFreq*tvec)';
nzerosamp = round(fs*subDurOff);
zvec = zeros(nzerosamp,1);
data = [data; zvec];
 
%--------------------------------------------------------------------------
function pulseTrain = makePulseTrain(pulse,burstdur,fs)

% check that pulse is charge balanced
cb = abs(sum(pulse));
if cb>0.000001
    disp('Pulse is not charge balanced')
    pulse = zeros(size(pulse));
end
sampperburst = round(fs*(burstdur/1000));
nReps = round(sampperburst/length(pulse));
%round(S.stim.sampperperiod/length(pulse));
pulseTrain = repmat(pulse,nReps,1);

%--------------------------------------------------------------------------
function pulse = makePulse(phase1pulsewidth,phase2pulsewidth,supthresAmp,freq,fs)

sampperphase1 = round(fs*(phase1pulsewidth/1e6));
sampperphase2 = round(fs*(phase2pulsewidth/1e6));
part1 =  supthresAmp*ones(sampperphase1,1);
part2 =  -supthresAmp*ones(sampperphase2,1);
pulse = [part1; part2];
pulseperiod = 1/freq;
sampperpulse = round(fs*pulseperiod);
zerosamps = sampperpulse-length(pulse);
zvec = zeros(zerosamps,1);
pulse = [pulse; zvec];
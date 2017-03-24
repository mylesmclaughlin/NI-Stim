function [A,sinwave] = pulsedSine(pulsewidth,pulserate,sinfreq,duration,fs)

if length(duration)==2
    totaldur = duration(2);
    duration = duration(1);
    addzero = 1;
else
    addzero = 0;
end
pulse = makePulse(pulsewidth*1e6,0,1,pulserate,fs);
pulseTrain = makePulseTrain(pulse,duration*1e3,fs);
tvec = [1/fs:1/fs:duration]';
sinwave = sin(2*pi*sinfreq*tvec);
sig = sinwave.*pulseTrain;

if addzero
    totalsamp = totaldur*fs
    zvec = zeros(totalsamp,1);
    zvec(1:length(sig)) = sig;
    sig = zvec;
    tvec = [1/fs:1/fs:totaldur]';
end

triggerdur = 10; % ms
triggeramp = 1;
trigger = zeros(length(sig),1);
nTrigSamps = round(triggerdur*1e-3*fs);
trigger(1:nTrigSamps) = triggeramp;

figure
plot(tvec,sinwave,'r')
hold on
plot(tvec,pulseTrain,'k')
plot(tvec,sig,'b')
hold on
plot(tvec,trigger,'m')

A.fs = fs;
A.stim.data(1,:,1) = trigger;
A.stim.data(1,:,2) = sig;

% save data
pathname = 'C:\Users\Public\MATLAB\NI-Stim\Stimuli\';
filename = ['pulsedSin__sinfreq=' num2str(sinfreq) '_pulsewidth=' num2str(pulsewidth) '_pulserate=' num2str(pulserate) '_fs=' num2str(fs) '.mat'];
disp(['Saving ' pathname filename])
save([pathname filename],'A')

%--------------------------------------------------------------------------
function pulseTrain = makePulseTrain(pulse,burstdur,fs)

% check that pulse is charge balanced
cb = abs(sum(pulse));
if cb>0.000001
    disp('Pulse is not charge balanced')
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
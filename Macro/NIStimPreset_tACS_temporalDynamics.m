function P = NIStimPreset_tACS_temporalDynamics


%% - Setup the marco struture -

% Setup the macro variables
P.preset.basefilename = 'RatTACSTempDynam-'; % set base file name
P.preset.basefilename = [P.preset.basefilename  strrep(strrep(datestr(now),':','-'),' ','-')]; % add date to make filename unique
P.preset.record = 1; % 1 = play stimulus and record data. 0 = play stimulus only;

%% Set fixed variables for easy modification

% sub-threshold (base) settings
subThresAmp = -0.05;
subThresFreq = 4;
subThresDC = 1;
subThresDurOn = 300;
subThresDurOff = 180;
subThresOffFirst = 0;
subThresNumberreps = 1;
disp(['This stimulus will take ' num2str((subThresDurOn + subThresDurOff)/60*subThresNumberreps) ' minutes'])

% supra-threshold (probe pulse train) settings
phaseDelay = 0; % phaseDelay in function of subthreshold stim frequency - 0 = 0; 1 = 2*pi 
supThresBurstdelay = 1000*(1/subThresFreq)*phaseDelay;
supThresAmp = 2.8;
supThresFreq = 300;
supThresBurstdur = 100;
supThresBurstrepperiod = 30000; % ms
supThresWaveformindex = 2; % biphasic pulse
supThresPhase1pulsewidth = 200;
supThresPhase2pulsewidth = 200;

% calculate NIstim parameters
nSeqOn = subThresDurOn/(supThresBurstrepperiod/1000);
nSeqOff = subThresDurOff/(supThresBurstrepperiod/1000);
if subThresOffFirst == 1
    subThresAmpSeq = [zeros(1,nSeqOff) subThresAmp*ones(1,nSeqOn)];
elseif subThresOffFirst == 0
    subThresAmpSeq = [subThresAmp*ones(1,nSeqOn) zeros(1,nSeqOff)];
end

% store settings for stimulus reconstruction 
P.preset.subThresAmpSeq = subThresAmpSeq;
P.preset.subThresFreq = subThresFreq;
P.preset.subThresDC = subThresDC;
P.preset.subThresDurOn = subThresDurOn;
P.preset.subThresDurOff = subThresDurOff;
P.preset.subThresOffFirst = subThresOffFirst;
P.preset.subThresNumberreps = subThresNumberreps;


%--------- make preset stimulus and base stimulus settings ---------
P.stim.stim = 0;
P.stim.amplitude = supThresAmp;
P.stim.frequency = supThresFreq;
P.stim.phase = 0;
P.stim.continuous = 0;
P.stim.ramp = 0;
P.stim.burstdur = supThresBurstdur;
P.stim.burstrepperiod = supThresBurstrepperiod;
P.stim.numberreps = subThresNumberreps;
P.stim.repsplayed = 0;
P.stim.ampmoddepth = 0;
P.stim.ampmodfreq = 2;
P.stim.ampmodphase = 0;
P.stim.waveformindex = supThresWaveformindex;
P.stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
P.stim.phase1pulsewidth = supThresPhase1pulsewidth;
P.stim.phase2pulsewidth = supThresPhase2pulsewidth;
P.stim.phase1amp = 100;
P.stim.phase2amp = -100;
P.stim.phasegap = 0;
P.stim.sameonallchannels = 1;
P.stim.randomizesequence = 0;

P.stim.burstdelay = supThresBurstdelay;

P.basestim.stim = 1;
P.basestim.amplitude = subThresAmpSeq;
P.basestim.frequency = subThresFreq;
P.basestim.phase = 0;
P.basestim.burstdur = supThresBurstrepperiod;
P.basestim.dc = subThresDC;
P.basestim.ampmoddepth = 0;
P.basestim.ampmodfreq = 5;
P.basestim.ampmodphase = 0.5;
P.basestim.waveformindex = 1;
P.basestim.phase1pulsewidth = 50;
P.basestim.phase2pulsewidth = 50;
P.basestim.phase1amp = 100;
P.basestim.phase2amp = -100;
P.basestim.phasegap = 0;

%--- Play preset stimulus and record data
NIStimPreset(P)

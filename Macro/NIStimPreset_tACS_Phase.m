function P = NIStimPreset_tACS_Phase


%% - Setup the marco struture -

% Setup the macro variables
P.preset.basefilename = 'RatTACSPhase-'; % set base file name
P.preset.basefilename = [P.preset.basefilename  strrep(strrep(datestr(now),':','-'),' ','-')]; % add date to make filename unique
P.preset.record = 1; % 1 = play stimulus and record data. 0 = play stimulus only;

%% Set fixed variables for easy modification

% sub-threshold (base) settings
subThresAmp = 0.2;
subThresFreq = 2; %[0.5 1 2 4 8 16 32 64];
subThresDC = 0;
%subThresNcycles = 5;
subThresNumberreps = 6;
subThresDurOn = 2; %subThresNcycles*(1/subThresFreq);
subThresDurOff = 0;

% supra-threshold (probe pulse train) settings
phaseDelay = [0:0.125:1]; % phaseDelay in function of subthreshold stim frequency - 0 = 0; 1 = 2*pi 
supThresBurstDelaySeq = 1000*(1/subThresFreq)*phaseDelay;
supThresAmp = 2.8;
supThresFreq = 300;
supThresBurstdur = 30;
supThresBurstrepperiod = 1000*(subThresDurOn+subThresDurOff); % ms
supThresWaveformindex = 2; % biphasic pulse
supThresPhase1pulsewidth = 50;
supThresPhase2pulsewidth = 50;

% specify how to combine supra and sub stimuli
stimcombinemethod = 'add-zerocenter'; %'add-zerocenter'; %'add';  %'add-zerocenter' , 'stop-insert'
makebasezero = 1;

disp(['This stimulus will take ' num2str((subThresDurOn + subThresDurOff)/60*subThresNumberreps*length(supThresBurstDelaySeq)) ' minutes'])


randomizesequence = 1;
% store settings for stimulus reconstruction 
P.preset.subThresFreq = subThresFreq;
P.preset.subThresDC = subThresDC;
P.preset.subThresDurOn = subThresDurOn;
P.preset.subThresDurOff = subThresDurOff;
P.preset.subThresNumberreps = subThresNumberreps;
P.preset.randomizesequence = randomizesequence;

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
P.stim.randomizesequence = randomizesequence;

P.stim.burstdelay = supThresBurstDelaySeq;

P.basestim.stim = 1;
P.basestim.amplitude = subThresAmp;
P.basestim.frequency = subThresFreq;
P.basestim.phase = 0;
P.basestim.burstdur = subThresDurOn*1000;
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
P.basestim.stimcombinemethod = stimcombinemethod;
P.basestim.makebasezero = makebasezero;

%--- Play preset stimulus and record data
NIStimPreset(P)

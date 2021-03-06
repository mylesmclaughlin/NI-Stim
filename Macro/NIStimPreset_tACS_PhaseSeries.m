function P = NIStimPreset_tACS_PhaseSeries

%% - Setup the marco struture -

% Setup the macro variables
P.preset.basefilename = 'RatTACSPhaseSeries-'; % set base file name
P.preset.basefilename = [P.preset.basefilename  strrep(strrep(datestr(now),':','-'),' ','-')]; % add date to make filename unique
P.preset.record = 1; % 1 = play stimulus and record data. 0 = play stimulus only;

%% Set fixed variables for easy modification

% sub-threshold (base) settings
subThresAmp = 0.4; %[0 0.2 0.4 0.6]; %0.8;
subThresFreq = 2; %[1 2 4 8 16 32 64 128];
subThresDC = 0;
%subThresNcycles = 5;
subThresNumberreps = 5;

% supra-threshold (probe pulse train) settings
phaseDelay =  [0:0.125:1];% [0:0.25:1]; % %[0:0.0025:0.01]; % % phaseDelay in function of subthreshold stim frequency - 0 = 0; 1 = 2*pi 
makebasezero = 0; % set to 1 to collect data with no TACS 
supThresBurstDelaySeq = 1000*(1/subThresFreq)*phaseDelay;
supThresAmp = [3.8 3.8 3.8];
supThresSeriesRepPeriod = 2000; % amplitudes above are presented xxx ms apart

supThresFreq = 300;
supThresBurstdur = 34;
% supThresFreq = 400;
% supThresBurstdur = 25;
% supThresFreq = 500;
% supThresBurstdur = 20;
% supThresFreq = 600;
% supThresBurstdur = 17;
% supThresFreq = 700;
% supThresBurstdur = 14.5;
% supThresFreq = 800;
% supThresBurstdur = 13;
% supThresFreq = 900;
% supThresBurstdur = 11.5;
% supThresFreq = 1000;
% supThresBurstdur = 10;

supThresWaveformindex = 6; % pulse-series
supThresPhase1pulsewidth = 200;
supThresPhase2pulsewidth = 200;

subThresDurOn = length(supThresAmp) * supThresSeriesRepPeriod/1000; %subThresNcycles*(1/subThresFreq);
subThresDurOff = 0;

supThresBurstrepperiod = 1000*(subThresDurOn+subThresDurOff); % ms

% specify how to combine supra and sub stimuli
stimcombinemethod = 'add-zerocenter';  % 'add', 'add-zerocenter' ,'stop-insert';  

disp(['This stimulus will take ' num2str((subThresDurOn + subThresDurOff)/60*subThresNumberreps*length(supThresBurstDelaySeq)*length(subThresAmp)) ' minutes'])

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
P.stim.amplitude = supThresAmp(1);
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
P.stim.waveformlist = {'sine','pulse','triangle','gaussian','custom','pulse-series','triangle-series','gaussian-series'};  % pulse, custom
P.stim.phase1pulsewidth = supThresPhase1pulsewidth;
P.stim.phase2pulsewidth = supThresPhase2pulsewidth;
P.stim.phase1amp = 100;
P.stim.phase2amp = -100;
P.stim.phasegap = 0;
P.stim.sameonallchannels = 1;
P.stim.randomizesequence = randomizesequence;
P.stim.series.amplitude = supThresAmp;
P.stim.series.burstrepperiod = supThresSeriesRepPeriod;

P.stim.burstdelay = supThresBurstDelaySeq;
P.stim.burstphasedelay = 0;

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

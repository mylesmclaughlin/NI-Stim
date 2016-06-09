function M = NIStimMacro_PulseShape


%% - Setup the marco struture -

% Setup the macro variables
M(1).macro.basefilename = 'Rat-PulseShape-'; % set base file name
M(1).macro.basefilename = [M(1).macro.basefilename  strrep(strrep(datestr(now),':','-'),' ','-')]; % add date to make filename unique
M(1).macro.nreps = 5; % - number of time to repeated each stimulus contained in the macro structure
M(1).macro.randomize = 1;  % randomize the order of stimului presentation
M(1).macro.record = 1; % record data (separate file for each rep of macro)

M(1).macro.allinonebinfile = 1;
% if M(1).macro.record == 0
%     M(1).macro.allinonebinfile = 1; % specify if the recorded data will be stored in one continuous .bin file ...
% elseif M(1).macro.record == 1
%     M(1).macro.allinonebinfile = 0; %... or if each element in the array will have its own .bin file
% end

% Set fixed variables for easy modification
ampRange = 1.7:.1:2.6;
frequency = 300;
burstdur = 100;
burstrepperiod = 2000;
phase1pulsewidth = 100;
phase2pulsewidth = 100;
phasegap = 100;
numberreps = 1;

% --Setup stimulation variables in a struct array --
% Stimulus 1 - Biphasic pulse
M(1).stim.seqname = 'BP';
M(1).stim.stim = 0;
M(1).stim.amplitude = ampRange;
M(1).stim.frequency = frequency;
M(1).stim.phase = 0;
M(1).stim.continuous = 0;
M(1).stim.ramp = 0;
M(1).stim.burstdur = burstdur;
M(1).stim.burstrepperiod = burstrepperiod;
M(1).stim.numberreps = numberreps;
M(1).stim.repsplayed = 0;
M(1).stim.ampmoddepth = 0;
M(1).stim.ampmodfreq = 5;
M(1).stim.ampmodphase = 0;
M(1).stim.waveformindex = 2;
M(1).stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
M(1).stim.phase1pulsewidth = phase1pulsewidth;
M(1).stim.phase2pulsewidth = phase2pulsewidth;
M(1).stim.phase1amp = 100;
M(1).stim.phase2amp = -100;
M(1).stim.phasegap = 0;
M(1).stim.sameonallchannels = 1;
M(1).stim.randomizesequence = 0;

% Stimulus 2 - Inter phase gap
M(2).stim.seqname = 'IPG';
M(2).stim.stim = 0;
M(2).stim.amplitude = ampRange;
M(2).stim.frequency = frequency;
M(2).stim.phase = 0;
M(2).stim.continuous = 0;
M(2).stim.ramp = 0;
M(2).stim.burstdur = burstdur;
M(2).stim.burstrepperiod = burstrepperiod;
M(2).stim.numberreps = numberreps;
M(2).stim.repsplayed = 0;
M(2).stim.ampmoddepth = 0;
M(2).stim.ampmodfreq = 5;
M(2).stim.ampmodphase = 0;
M(2).stim.waveformindex = 2;
M(2).stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
M(2).stim.phase1pulsewidth = phase1pulsewidth;
M(2).stim.phase2pulsewidth = phase2pulsewidth;
M(2).stim.phase1amp = 100;
M(2).stim.phase2amp = -100;
M(2).stim.phasegap = phasegap;
M(2).stim.sameonallchannels = 1;
M(2).stim.randomizesequence = 0;

% Stimulus 3 - Psuedo Monophasic
M(3).stim.seqname = 'PM';
M(3).stim.stim = 0;
M(3).stim.amplitude = ampRange;
M(3).stim.frequency = frequency;
M(3).stim.phase = 0;
M(3).stim.continuous = 0;
M(3).stim.ramp = 0;
M(3).stim.burstdur = burstdur;
M(3).stim.burstrepperiod = burstrepperiod;
M(3).stim.numberreps = numberreps;
M(3).stim.repsplayed = 0;
M(3).stim.ampmoddepth = 0;
M(3).stim.ampmodfreq = 5;
M(3).stim.ampmodphase = 0;
M(3).stim.waveformindex = 2;
M(3).stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
M(3).stim.phase1pulsewidth = phase1pulsewidth;
M(3).stim.phase2pulsewidth = phase2pulsewidth*10;
M(3).stim.phase1amp = 100;
M(3).stim.phase2amp = -10;
M(3).stim.phasegap = 0;
M(3).stim.sameonallchannels = 1;
M(3).stim.randomizesequence = 0;

% Stimulus 4 - Triangle
M(4).stim.seqname = 'TRI';
M(4).stim.stim = 0;
M(4).stim.amplitude = ampRange;
M(4).stim.frequency = frequency;
M(4).stim.phase = 0;
M(4).stim.continuous = 0;
M(4).stim.ramp = 0;
M(4).stim.burstdur = burstdur;
M(4).stim.burstrepperiod = burstrepperiod;
M(4).stim.numberreps = numberreps;
M(4).stim.repsplayed = 0;
M(4).stim.ampmoddepth = 0;
M(4).stim.ampmodfreq = 5;
M(4).stim.ampmodphase = 0;
M(4).stim.waveformindex = 3;
M(4).stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
M(4).stim.phase1pulsewidth = phase1pulsewidth;
M(4).stim.phase2pulsewidth = phase2pulsewidth;
M(4).stim.phase1amp = 100;
M(4).stim.phase2amp = -100;
M(4).stim.phasegap = phasegap;
M(4).stim.sameonallchannels = 1;
M(4).stim.randomizesequence = 0;

% Stimulus 5 - Gauss
M(5).stim.seqname = 'GAUS';
M(5).stim.stim = 0;
M(5).stim.amplitude = ampRange;
M(5).stim.frequency = frequency;
M(5).stim.phase = 0;
M(5).stim.continuous = 0;
M(5).stim.ramp = 0;
M(5).stim.burstdur = burstdur;
M(5).stim.burstrepperiod = burstrepperiod;
M(5).stim.numberreps = numberreps;
M(5).stim.repsplayed = 0;
M(5).stim.ampmoddepth = 0;
M(5).stim.ampmodfreq = 5;
M(5).stim.ampmodphase = 0;
M(5).stim.waveformindex = 4;
M(5).stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
M(5).stim.phase1pulsewidth = phase1pulsewidth;
M(5).stim.phase2pulsewidth = phase2pulsewidth;
M(5).stim.phase1amp = 100;
M(5).stim.phase2amp = -100;
M(5).stim.phasegap = phasegap;
M(5).stim.sameonallchannels = 1;
M(5).stim.randomizesequence = 0;

%% - Run the macro -
M = NIStimMacro(M);

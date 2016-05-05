function NIStim(command,input1,input2)

if nargin == 0
    command = 'init';
end

if isobject(command)
    if input2 == 1
        command = 'update';
    elseif input2 == 2
        command = 'queueStim';
    elseif input2 ==3
        command = 'queueListen'
    end
end

switch command
    case 'init'
        NIinit
    case 'startStim'
        NIstartStim
    case 'updateStim'
        NIupdateStim
    case 'stopStim'
        NIstopStim
    case 'repCount'
        NIrepCount
    case 'queueStim'
        NIqueueStim(input1)
    case 'startRec'
        NIstartRec
    case 'stopRec'
        NIstopRec
    case 'startPlay'
        NIstartPlay
    case 'stopPlay'
        NIstopPlay
    case 'loadData'
        NIloadData
    case 'quickRec'
        NIquickRec
    case 'close'
        NIclose
    case 'continuous'
        NIcontinuous
    case 'waveform'
        NIwaveform
    case 'cutsomwaveform'
        NIcustomwaveform
    case 'pulsewidth'
        NIpulsewidth
    case 'frequency'
        NIfrequency
    case 'playbackmode'
        NIplaybackmode
    case 'accelcalib'
        NIaccelcalib
    case 'processmode'
        NIprocessmode
    case 'showrecplot'
        NIshowrecplot
    case 'stimRecProc'
        NIstimRecProc
    case 'stimRec'
        NIstimRec
    case 'stimProc'
        NIstimProc
    case 'stimFileName'
        NIstimFileName
    case 'recFileName'
        NIrecFileName
    case 'changeStimulator'
        NIchangeStimulator
    case 'listenOn'
        if nargin<2
            input1 = [];
        end
        NIlistenOn(input1)
    case 'listenOff'
        NIlistenOff
    case 'queueListen'
        NIqueueListen(input1)
    case 'update'
        NIupdate(input1)
end

%-------------------------------------------------------------------------
function NIinit
global S NI LH RH AO

% Get settings
S = NIStimSettings;

% Create NI object
[NI,LH,RH,S] = NIcreateObj(S);

% Make plot
NIStimFig;

% Connect to current source
if S.current.present == 1
    D = CurrentControl('connect',S.current.name);
    D = CurrentControl('disableOutput',S.current.name);
end

% Start NI object
if S.ni.connected == 1
    if ~NI.IsRunning
        NI.startBackground();
    end
end

% Autozero current source
if S.current.present == 1
    CurrentControl('autoZero',S.current.name);
end

% Make AO for listening to audio output
% AO = NIlistenCreatObj;

% Try parallel processing
% delete(gcp)
% P = parpool;

%-------------------------------------------------------------------------
function [NI,LH,RH,S] = NIcreateObj(S)

if S.ni.connected == 1
    NI = daq.createSession('ni');
    
    % add input channels
    chIn = addAnalogInputChannel(NI,S.ni.devname,S.ni.chin, 'Voltage');
    %set(NI.Channels,'InputType',  S.ni.inputtype);
    for n = 1:S.ni.nchin
        set(NI.Channels(n),'InputType',S.ni.inputtype{n});
        set(NI.Channels(n),'Range',S.ni.voltrange(n,:));
    end
    
    % add output channels
    if S.stim.stim == 1
        chOut = addAnalogOutputChannel(NI,S.ni.devname,S.ni.chout, 'Voltage');
        RH = addlistener(NI,'DataRequired',@(src,event) NIStim(src,event,2));
        NI.NotifyWhenScansQueuedBelow = S.ni.buffersize;
        S.stim.data = zeros(S.ni.buffersize,S.ni.nchout);
        queueOutputData(NI,S.stim.data);
        queueOutputData(NI,S.stim.data);
    else
        RH = [];
    end
    
    NI.Rate = S.ni.rate;
    NI.IsContinuous = 1;
    %LH = addlistener(NI,'DataAvailable',@NIStim);
    LH = addlistener(NI,'DataAvailable',@(src,event) NIStim(src,event,1));
    NI.NotifyWhenDataAvailableExceeds = S.ni.buffersize;
else
    NI.Rate = S.ni.rate;
    LH = [];
    RH = [];
end

%-------------------------------------------------------------------------
function fullfileName = NImakeFileName(fullfileName)

[path,fname,ext] = fileparts(fullfileName);
refFile = [path '\' fname];
extList = {'.bin','.nis','.niso','.proc','.qs'};
if exist([refFile extList{1}],'file') | exist([refFile extList{2}],'file') | exist([refFile extList{3}],'file') | exist([refFile extList{4}],'file')
    uind = strfind(fname,'_');
    if length(uind)>1
        fname = fname(1:uind(end)-1);
    end
    m = 0;
    testfile = [path '\' fname];
    while exist([testfile extList{1}],'file') | exist([testfile extList{2}],'file') | exist([testfile extList{3}],'file') | exist([refFile extList{4}],'file')
        m=m+1;
        testfile = [path '\' fname '_' num2str(m)];
    end
    fullfileName = [testfile ext];
end

%-------------------------------------------------------------------------
function NIstartRec
global NI S
S.rec.rec = 1;
S.rec.filename = get(S.rec.filenamebut,'String');
S.rec.fullfileName = [S.data.dir S.rec.filename '.bin'];

% check file name and append 1 if it already exists
S.rec.fullfileName = NImakeFileName(S.rec.fullfileName);

set(S.fig.fhdl,'name',['NIStim - Recording ' S.rec.fullfileName])
set(S.rec.recbut,'String','Stop Rec','Callback','NIStim(''stopRec'')','backgroundcolor',[1 0 0]);

NIsaveStimParam

NIstopWatch('start'); % restart timer

%-------------------------------------------------------------------------
function NIsaveStimParam(stimOnly)
global NI S

if nargin== 0
    stimOnly = 0;
end
% save file with all settings
%NIparseStimGUI; % make sure latest stim param are saved
D.timestamp = datestr(now);
D.ni = S.ni;
D.amp = S.amp;
D.current = S.current;
D.data = S.data;
D.stim = S.stim;
D.basestim = S.stim;
D.accel = S.accel;
sequence = S.sequence;
sequence.data = [];
D.sequence = sequence;

[path,fname,ext] = fileparts(S.rec.fullfileName);
if stimOnly == 1
    save([path '\' fname '.niso'],'D','-mat');
else
    save([path '\' fname '.nis'],'D','-mat');
end

%-------------------------------------------------------------------------
function NIstopRec
global NI S
S.rec.rec = 0;
set(S.rec.recbut,'String','Record','Callback','NIStim(''startRec'')','backgroundcolor',[0 1 0]);
set(S.fig.fhdl,'name','NIStim')
NIstopWatch('stop')

%-------------------------------------------------------------------------
function NIstartPlay
global NI S T

S.rec.filename = get(S.rec.filenamebut,'String');
S.rec.fullfileName = [S.data.dir S.rec.filename '.bin'];
if ~exist(S.rec.fullfileName,'file')
    disp(['File not found: ' S.rec.fullfileName])
    return
end
set(S.rec.recbut,'string','Stop play','callback','NIStim(''stopPlay'')','backgroundcolor',[1 0 0])
S.rec.playind = 1;
S.rec.rec = 0;
S.stim.stim = 0;

set(S.fig.fhdl,'name',['NIStim - Playing ' S.rec.fullfileName])

T = timer;
T.ExecutionMode = 'fixedRate';
T.period = S.ni.updateperiod;
T.TimerFcn = 'NIStim(''loadData'')';
start(T);

%-------------------------------------------------------------------------
function NIloadData
global S

D = binread(S.rec.fullfileName,[S.rec.playind S.rec.playind+S.ni.buffersize-1]);
if isstruct(D)
    S.rec.playind = S.rec.playind+S.ni.buffersize;
    event.Data = D.data;
    NIupdate(event)
else
    NIstopPlay
end

%-------------------------------------------------------------------------
function NIstopPlay
global NI S T

set(S.rec.recbut,'string','Play','callback','NIStim(''startPlay'')','backgroundcolor',[0 1 0])
set(S.fig.fhdl,'name','NIStim')
stop(T)
delete(T)

%-------------------------------------------------------------------------
function NIstartStim
global NI S RH

% Make stimulus
S.stim.stim = 1;
NIparseStimGUI
NImakeStim
S.stim.firstqueue = 1;

if S.sequence.on == 1 
    if ~strcmpi(S.stim.waveformlist(S.stim.waveformindex),'custom') % custom files should already contain sequenced data
        NImakeSeqStim
    end
end
try
    NIsaveStimParam % make sure correct data is saved
end

% Add channels
NI.stop
chOut = addAnalogOutputChannel(NI,S.ni.devname,S.ni.chout, 'Voltage');
RH = addlistener(NI,'DataRequired',@(src,event) NIStim(src,event,2));
NI.NotifyWhenScansQueuedBelow = S.stim.buffersize; % update buffer size

% reset processing
if S.proc.proc == 1
    NIresetProcessedData
end

% Queue zero data
holdbuffer = S.stim.data;
S.stim.data = zeros(size(S.stim.data));
%if S.stim.buffersize/S.ni.rate<10;
    queueOutputData(NI,S.stim.data);
    queueOutputData(NI,S.stim.data);
%else
%    queueOutputData(NI,S.stim.data);
%end

% Start NI card
NI.startBackground();

% Enable current source
if S.current.present == 1
    D = CurrentControl('autoZero',S.current.name);
    pause(0.5)
    D = CurrentControl('enableOutput',S.current.name);
    pause(0.5)
end

% Queue stim data
S.stim.starttime = datestr(now);
S.stim.data = holdbuffer;

% Update GUI
set(S.stim.startbut,'String','Stop Stim','Callback','NIStim(''stopStim'')','backgroundcolor',[1 0 0])
set(S.rec.startstimbut,'String','Stop Stim','Callback','NIStim(''stopStim'')','backgroundcolor',[1 0 0])
set(S.stim.updatebut,'enable','on')
set(S.stim.stimrecbut,'enable','off')
set(S.stim.stimprocbut,'enable','off')
set(S.stim.stimrecprocbut,'enable','off')

% save current stim parameters
NIwriteCurrentStimParam

if S.rec.rec == 0
    % params normally saved with recording
    % but if no recording save a record anyway
    
    S.rec.filename = get(S.rec.filenamebut,'String');
    S.rec.fullfileName = [S.data.dir S.rec.filename '.bin'];
    % check file name and append 1 if it already exists
    S.rec.fullfileName = NImakeFileName(S.rec.fullfileName);
    
    stimOnly = 1;
    NIsaveStimParam(stimOnly);
end
%-------------------------------------------------------------------------
function NIupdateStim
global NI S

% Make stimulus
NIparseStimGUI
NImakeStim

% Queue stim
queueOutputData(NI,S.stim.transitiondata);
queueOutputData(NI,S.stim.data);

% save current stim parameters
NIwriteCurrentStimParam

%-------------------------------------------------------------------------
function NIstopStim
global NI S LH RH

set(S.stim.startbut,'String','Stopping...','enable','off')
set(S.rec.startstimbut,'String','Stopping...','enable','off')
S.stim.stim = 0;
S.stim.current.amplitude = 0; % set to zero for ramp up on next stim
S.stim.endtime = datestr(now);

% Save any processed data processing
% if S.proc.proc == 1
%     S.proc.proc = 0;
%     set(S.proc.procbut,'value',S.proc.proc);
%     NIsaveProcessedData
% end

% Make sure to finish on zero
% delete(RH);
% queueOutputData(NI,zeros(size(S.stim.data)));
% lastQueued = NI.ScansQueued;
% while NI.ScansQueued > 1
%     disp(['Scans queued = ' num2str(NI.ScansQueued) '; Last Queued = ' num2str(lastQueued)]);
%     pause(0.1)
%     lastQueued = NI.ScansQueued;
% end
% disp(['Scans queued = ' num2str(NI.ScansQueued) '; Last Queued = ' num2str(lastQueued)]);
% disp('Out of while loop')
holdbuffer = S.stim.data;

S.stim.data = zeros(size(S.stim.data));
bufferDur = S.stim.buffersize/NI.Rate;
queueOutputData(NI,zeros(size(S.stim.data))); % just added 02/02/2016 - TEST
pause(bufferDur*2)

% Disable current source
if S.current.present == 1
    D = CurrentControl('disableOutput',S.current.name);
end

% Stop recording
% if S.rec.rec == 1;
%     NIstopRec
% end

%disp('1')
% Delete channels
NI.stop
%disp('2')
NI.IsContinuous = 0;
%disp('3')
delete(RH)
%delete(LH)
%disp('4')

for n = length(NI.Channels):-1:1
    if strmatch('ao',NI.Channels(n).ID)
        removeChannel(NI,n)
    end
end
%disp('5')
NI.IsContinuous = 1;
%LH = addlistener(NI,'DataAvailable',@(src,event) NIStim(src,event,1));
%NI.NotifyWhenDataAvailableExceeds = S.ni.buffersize;
NI.startBackground();
%disp('6')
S.stim.data = holdbuffer;

% Update GUI
set(S.stim.startbut,'String','Stimulate','Callback','NIStim(''startStim'')','backgroundcolor',[0 1 0],'enable','on')
set(S.rec.startstimbut,'String','Stimulate','Callback','NIStim(''startStim'')','backgroundcolor',[0 1 0],'enable','on')
set(S.stim.updatebut,'enable','off')
set(S.stim.stimrecbut,'enable','on')
set(S.stim.stimprocbut,'enable','on')
set(S.stim.stimrecprocbut,'enable','on')

%-------------------------------------------------------------------------
function NIqueueStim(event)
global S NI

% update stim to next sequence
if S.sequence.on == 1
    S.stim.firstqueue = 0;
    S.sequence.thisseq = S.sequence.thisseq+1;
    if S.stim.numberreps==Inf % continuously loop through seqIndex
        S.sequence.loopthisseq = mod(S.sequence.thisseq,S.sequence.nseq);
        if S.sequence.loopthisseq == 0
            S.sequence.loopthisseq = S.sequence.nseq;
        end
    else
        S.sequence.loopthisseq = S.sequence.thisseq;
    end
    if S.sequence.loopthisseq>length(S.sequence.seqIndex)
        S.stim.data = zeros(size(S.sequence.data.stim1));
        S.sequence.thisseq = S.sequence.thisseq-1;
    else
        eval(['S.stim.data = S.sequence.data.stim' num2str(S.sequence.seqIndex(S.sequence.loopthisseq)) ';']);
    end
end

% queue stimulation data
if S.stim.firstqueue == 1
    queueOutputData(NI,S.stim.transitiondata);
    S.stim.firstqueue = 0;
else
    queueOutputData(NI,S.stim.data);
end

% keep track of reps
NIrepCount

%-------------------------------------------------------------------------
function NIrepCount
global S NI RH

% count reps here - more accurate
if S.stim.continuous == 0
    if S.sequence.on == 0
        S.stim.repsplayed = S.stim.repsplayed + 1;
        set(S.stim.repsplayedbut,'String',['rep ' num2str(S.stim.repsplayed)]);
    elseif S.sequence.on == 1
        S.stim.repsplayed = floor(S.sequence.thisseq/S.sequence.nseq);
        if S.sequence.loopthisseq <= length(S.sequence.seqIndex)
            set(S.stim.repsplayedbut,'String',['rep ' num2str(S.stim.repsplayed) ', seq ' num2str(S.sequence.seqIndex(S.sequence.loopthisseq))]);
        else
            set(S.stim.repsplayedbut,'String',['rep ' num2str(S.stim.repsplayed) ', seq ' num2str(S.sequence.seqIndex(end))]);
        end
    end
end

% %stop stim is all reps played
% if S.sequence.on == 0
%     if S.stim.repsplayed >= S.stim.numberreps
%         pause(0.1)
%         %NIstopStim
%     end
% elseif S.sequence.on == 1
%     if S.stim.numberreps ~= Inf
%         if S.sequence.thisseq >= length(S.sequence.seqIndex)
%             pause(0.1)
%             %NIstopStim
%         end
%     end
% end

%-------------------------------------------------------------------------
function NIparseStimGUI
global S

if S.current.present == 1
    S.stim.amplitude = str2num(get(S.stim.ampbut,'string'))/S.current.onevoltequalsXmilliamps;
elseif S.current.present == 0
    S.stim.amplitude = str2num(get(S.stim.ampbut,'string'));
end

S.stim.frequency = str2num(get(S.stim.freqbut,'string'));
S.stim.phase = str2num(get(S.stim.phasebut,'string'));
S.stim.continuous = get(S.stim.contbut,'value');
S.stim.ramp = get(S.stim.rampbut,'value');
S.stim.burstdur = str2num(get(S.stim.burstdurbut,'string'));
S.stim.burstrepperiod = str2num(get(S.stim.burstrepperiodbut,'string'));
S.stim.numberreps = str2num(get(S.stim.numberrepsbut,'string'));
S.stim.ampmoddepth = str2num(get(S.stim.ampmoddepthbut,'string'));
S.stim.ampmodfreq = str2num(get(S.stim.ampmodfreqbut,'string'));
S.stim.sameonallchannels = get(S.stim.samechanbut,'value');
S.stim.waveformindex = get(S.stim.waveformbut,'value');
S.stim.phase1pulsewidth = str2num(get(S.stim.phase1pwbut,'string'));
S.stim.phase2pulsewidth = str2num(get(S.stim.phase2pwbut,'string'));
S.stim.phase1amp = str2num(get(S.stim.phase1ampbut,'string'));
S.stim.phase2amp = str2num(get(S.stim.phase2ampbut,'string'));
S.stim.phasegap = str2num(get(S.stim.phasegapbut,'string'));

if ~strcmpi(S.stim.waveformlist(S.stim.waveformindex),'custom')
    SequenceFields = {'amplitude','frequency','phase','ampmoddepth','ampmodfreq','burstdelay'};
    hit = 0;
    for n = 1:length(SequenceFields)
        eval(['param = S.stim.' SequenceFields{n} ';']);
        fieldL = length(param);
        if fieldL>1
            S.sequence.on = 1;
            S.sequence.thisseq = 0;
            S.sequence.loopthisseq = 1;
            S.sequence.nseq = fieldL;
            S.sequence.parametername = SequenceFields{n};
            S.sequence.basestimseq = 0;
            S.sequence.parametervalues = param;
            eval(['S.stim.' SequenceFields{n} ' = S.sequence.parametervalues(1);']);
            S.sequence.seq = [1:S.sequence.nseq];
            if S.stim.numberreps == Inf
                S.sequence.seqIndex =  S.sequence.seq;
            else
                if S.stim.randomizesequence == 0
                    S.sequence.seqIndex =  repmat(S.sequence.seq,1,S.stim.numberreps);
                elseif S.stim.randomizesequence == 1
                    S.sequence.seqIndex = [];
                    for r = 1:S.stim.numberreps
                        rInd = randperm(S.sequence.nseq);
                        S.sequence.seqIndex =  [S.sequence.seqIndex S.sequence.seq(rInd)];
                    end
                end
            end
            hit = 1;
        end
    end
    
    if hit == 0 % look for sequence in basestimulus
        BaseSequenceFields = {'amplitude','frequency','phase','ampmoddepth','ampmodfreq'};
        for n = 1:length(BaseSequenceFields)
            eval(['param = S.basestim.' BaseSequenceFields{n} ';']);
            fieldL = length(param);
            if fieldL>1
                S.sequence.on = 1;
                S.sequence.thisseq = 0;
                S.sequence.loopthisseq = 1;
                S.sequence.nseq = fieldL;
                S.sequence.parametername = BaseSequenceFields{n};
                S.sequence.basestimseq = 1;
                S.sequence.parametervalues = param;
                eval(['S.basestim.' SequenceFields{n} ' = S.sequence.parametervalues(1);']);
                S.sequence.seq = [1:S.sequence.nseq];
                if S.stim.numberreps == Inf
                    S.sequence.seqIndex =  S.sequence.seq;
                else
                    S.sequence.seqIndex =  repmat(S.sequence.seq,1,S.stim.numberreps);
                end
                hit = 1;
            end
        end
    end
    
    if hit == 0;
        S.sequence.on = 0;
        S.sequence.thisseq = 0;
        S.sequence.loopthisseq = S.sequence.thisseq;
        S.sequence.nseq = 0;
        S.sequence.parametername = '';
        S.sequence.basestimseq = [];
        S.sequence.parametervalues = [];
        S.sequence.data = [];
        S.sequence.seq = [1:S.sequence.nseq];
        S.sequence.seqIndex =  S.sequence.seq;
        S.stim.randomizesequence = 0;
    end
end
%-------------------------------------------------------------------------
function NIwriteCurrentStimParam
global S

S.stim.current.amplitude = S.stim.amplitude;
S.stim.current.frequency = S.stim.frequency;
S.stim.current.phase = S.stim.phase;
S.stim.current.continuous = S.stim.continuous;
S.stim.current.burstdur = S.stim.burstdur;
S.stim.current.burstrepperiod = S.stim.burstrepperiod;
S.stim.current.ampmoddepth = S.stim.ampmoddepth;
S.stim.current.ampmodfreq = S.stim.ampmodfreq;
S.stim.current.waveformindex = S.stim.waveformindex;
S.stim.current.phase1pulsewidth = S.stim.phase1pulsewidth;
S.stim.current.phase2pulsewidth = S.stim.phase2pulsewidth;
S.stim.current.phase1amp = S.stim.phase1amp;
S.stim.current.phase2amp = S.stim.phase2amp;
S.stim.current.phasegap = S.stim.phasegap;

%-------------------------------------------------------------------------
function NImakeStim
global NI S

% determine stimulus period
if S.stim.continuous == 1
    if S.stim.ampmoddepth == 0
        S.stim.period = 1/S.stim.frequency;
    else
        if S.stim.ampmodfreq ~= 0
            S.stim.period = 1/S.stim.ampmodfreq;
        else
            S.stim.period = 1/S.stim.frequency;
        end
    end
else
    S.stim.period = (S.stim.burstrepperiod/1000);
end
S.stim.sampperperiod = round(NI.Rate*S.stim.period);


% adjust ni buffer size to  multiple of stimulus period
if S.stim.sampperperiod <= S.ni.buffersize
    S.stim.stimupfac = ceil(S.ni.buffersize/S.stim.sampperperiod);
    S.stim.sampperperiod =  S.stim.stimupfac*S.stim.sampperperiod;
end

% Make stim
if S.stim.continuous == 1
    tvec = [1/NI.Rate:1/NI.Rate:S.stim.sampperperiod/S.ni.rate];
    sampperburst = S.stim.sampperperiod;
    zvec = [];
    zvec_delay = [];
elseif S.stim.continuous == 0
    sampperburst = round(NI.Rate*(S.stim.burstdur/1000));
    zerosamps = S.stim.sampperperiod - sampperburst;
    if S.stim.burstdelay~= 0
        delaysamps = round(NI.Rate*(S.stim.burstdelay/1000));
        zerosamps = zerosamps - delaysamps;
        zvec_delay = zeros(delaysamps,1);
    else
        zvec_delay = [];
    end
    if S.stim.burstphasedelay ~=0
        phaseperiod = S.stim.burstphasedelay/S.basestim.frequency;
        phasedelaysamps = round(NI.Rate*phaseperiod);
        if S.stim.burstdelay == 0 % just use phase delay
            zerosamps = zerosamps - phasedelaysamps;
            zvec_delay = zeros(phasedelaysamps,1);
        else % combine phase and time delays
            cycleperiod = 1/S.basestim.frequency;
            cycledelaysamps = round(NI.Rate*cycleperiod);
            numcompletecyclesindelay = floor(delaysamps/cycledelaysamps);
            delaysamps = numcompletecyclesindelay*cycledelaysamps;
            delaysamps = delaysamps + phasedelaysamps;
            
            zerosamps = S.stim.sampperperiod - sampperburst;
            zerosamps = zerosamps - delaysamps;
            zvec_delay = zeros(delaysamps,1);
        end
    end
    
    zvec = zeros(zerosamps,1);
    tvec = [1/NI.Rate:1/NI.Rate:sampperburst/S.ni.rate];
end

if strcmpi(S.stim.waveformlist(S.stim.waveformindex),'sine')
    if S.stim.ampmoddepth == 0
        stimData = S.stim.amplitude*sin(2*pi*S.stim.frequency*tvec + S.stim.phase*2*pi)';
        
        %         amplitudeTrans = linspace(S.stim.current.amplitude,S.stim.amplitude,length(tvec));
        %         frequencyTrans = linspace(S.stim.current.frequency,S.stim.frequency,length(tvec));
        %         phaseTrans = linspace(S.stim.current.phase,S.stim.phase,length(tvec));
        %         stimDataTrans = (amplitudeTrans.*sin(2*pi*frequencyTrans.*tvec + phaseTrans*2*pi))';
        stimDataTrans = stimData;
    elseif S.stim.ampmoddepth ~= 0
        %AMphase = 0.5;
        
        stimData = [];
        stimData(:,1) = S.stim.amplitude*(1-S.stim.ampmoddepth/100/2)*sin(2*pi*S.stim.frequency                     *tvec + S.stim.phase*2*pi)';
        stimData(:,2) = S.stim.amplitude*(  S.stim.ampmoddepth/100/2)*sin(2*pi*(S.stim.frequency+S.stim.ampmodfreq) *tvec + S.stim.phase*2*pi + S.stim.ampmodphase*2*pi)';
        
        %         amplitudeTrans = linspace(S.stim.current.amplitude,S.stim.amplitude,length(tvec));
        %         frequencyTrans = linspace(S.stim.current.frequency,S.stim.frequency,length(tvec));
        %         phaseTrans = linspace(S.stim.current.phase,S.stim.phase,length(tvec));
        %
        %         stimDataTrans = [];
        %         stimDataTrans(:,1) = (amplitudeTrans.*(1-S.stim.ampmoddepth/100/2).*sin(2*pi*frequencyTrans.*tvec + phaseTrans*2*pi))';
        %         stimDataTrans(:,2) = (amplitudeTrans.*(  S.stim.ampmoddepth/100/2).*sin(2*pi*(frequencyTrans+S.stim.ampmodfreq).*tvec + phaseTrans*2*pi + S.stim.ampmodphase*2*pi))';
        stimDataTrans = stimData;
        if ~isempty(zvec)
            zvec(:,2) = zvec(:,1);
        end
        if ~isempty(zvec_delay)
            zvec_delay(:,2) = zvec_delay(:,1);
        end
        if S.stim.sameonallchannels == 1 %send both freqs on 1 channel
            stimData = sum(stimData')';
            stimDataTrans = sum(stimDataTrans')';
            if ~isempty(zvec)
                zvec = sum(zvec')';
            end
            if ~isempty(zvec_delay)
                zvec_delay = sum(zvec_delay')';
            end
        end
    end
    
    
elseif sum(strcmpi(S.stim.waveformlist(S.stim.waveformindex),{'pulse','triangle','gaussian'}))
    pulse = NImakepulse;
    
    % check that pulse is charge balanced
    cb = abs(sum(pulse));
    if cb>0.000001
        disp('Pulse is not charge balanced')
        pulse = zeros(size(pulse));
    end
    nReps = round(sampperburst/length(pulse));
    %round(S.stim.sampperperiod/length(pulse));
    stimData = repmat(pulse,nReps,1);
    %stimDataTrans = repmat(pulse,nReps,1);
    if S.ampblank.on == 1;
        S.ampblank.pulsetrain = repmat(S.ampblank.singlepulse,nReps,1);
    end
    
    if S.stim.ampmoddepth ~= 0
        stimData = (1+(S.stim.ampmoddepth/100)*cos(2*pi*S.stim.ampmodfreq*tvec + S.stim.ampmodphase*2*pi))'.*stimData;
        stimData = S.stim.amplitude*stimData/max(stimData);
    end
    stimDataTrans = stimData;
elseif sum(strcmpi(S.stim.waveformlist(S.stim.waveformindex),{'pulse-series','triangle-series','gaussian-series'}))
    
    pulse = NImakepulse;
    
    % check that pulse is charge balanced
    cb = abs(sum(pulse));
    if cb>0.000001
        disp('Pulse is not charge balanced')
        pulse = zeros(size(pulse));
    end
    nReps = round(sampperburst/length(pulse));
    %round(S.stim.sampperperiod/length(pulse));
    pulseTrain = repmat(pulse,nReps,1);
    pulseTrain  = pulseTrain /max(pulseTrain );
    
    %stimDataTrans = repmat(pulse,nReps,1);
    if S.ampblank.on == 1;
        S.ampblank.pulsetrain = repmat(S.ampblank.singlepulse,nReps,1);
    end
    
    sampperseries = round(NI.Rate*(S.stim.series.burstrepperiod/1000)); 
    seriesvec = zeros(sampperseries,1);
    stimData = [];
    for n = 1:length(S.stim.series.amplitude)
        locStimData = seriesvec;
        locStimData(1:length(pulseTrain)) = (S.stim.series.amplitude(n)/S.current.onevoltequalsXmilliamps) * pulseTrain; 
        stimPulseSeriesInd(1,n) = length(stimData)+1;
        stimPulseSeriesInd(2,n) = stimPulseSeriesInd(1,n) + length(pulseTrain);
        stimData = [stimData; locStimData];
    end
    
    stimDataTrans = stimData;
    
elseif strcmpi(S.stim.waveformlist(S.stim.waveformindex),'custom')
    stimData = squeeze(S.stim.customdata(1,:,:));
    if S.ni.nchout > 1
        stimData(:,2) = S.stim.amplitude*(stimData(:,2)/max(stimData(:,2)));
    else
        stimData = stimData';
        stimData(:,1) = S.stim.amplitude*(stimData(:,1)/max(stimData(:,1)));
    end
    stimDataTrans = stimData;
end

% apply ramp
if S.stim.ramp == 1
    amplitudeTrans = linspace(0,1,length(tvec))';
    stimDataTrans = amplitudeTrans.*stimDataTrans;
end

if ~strcmpi(S.stim.waveformlist(S.stim.waveformindex),'custom')
    if S.basestim.stim == 1;
        % add base stimulus
        if ~exist('stimPulseSeriesInd','var'); stimPulseSeriesInd = []; end
        [stimData,stimDataTrans] = NIaddBaseStim(stimData,stimDataTrans,zvec_delay,stimPulseSeriesInd);
    else
        % or add zero buffers
        stimData = [zvec_delay; stimData; zvec];
        stimDataTrans = [zvec_delay; stimDataTrans; zvec];
    end
end

if S.ampblank.on == 1;
    S.ampblank.pulsetrain = [zvec_delay; S.ampblank.pulsetrain; zvec]; 
end

% check charge balance
cb = sum(stimData(:));
if cb > 0.5
    disp(['WARNING: Stimulus is not charge balanced - ' num2str(cb) ' samples or ' num2str((cb*1/NI.Rate)*1000) 'mA offset' ])
end
cb = sum(stimDataTrans(:));
if cb > 0.5
    disp(['WARNING: Stimulus is not charge balanced - ' num2str(cb) ' samples or ' num2str((cb*1/NI.Rate)*1000) 'mA offset' ])
end

S.stim.data = [];
S.stim.transitiondata = [];

if strcmpi(S.stim.waveformlist(S.stim.waveformindex),'custom')
    S.stim.data = stimData;
    S.stim.transitiondata = stimDataTrans;
else
    nD = size(stimData,2);
    if nD == 1
        for n = 1:S.ni.nchout
            S.stim.data = [S.stim.data stimData];
            S.stim.transitiondata = [S.stim.transitiondata stimDataTrans];
        end
        if S.trigger.on == 1
            triggerData = zeros(length(stimData),1);
            nTrigSamps = round(S.trigger.dur*1e-3*NI.Rate);
            triggerData(1:nTrigSamps) = S.trigger.amp;
            S.stim.data(:,1) = triggerData;
            S.stim.transitiondata(:,1) = triggerData;
        end
        if S.ampblank.on == 1;
            for n = 1:length(S.ampblank.pulsechind)
                S.stim.data(:,S.ampblank.pulsechind(n)) = S.ampblank.pulsetrain;
                S.stim.transitiondata(:,S.ampblank.pulsechind(n)) = S.ampblank.pulsetrain;
            end
        end
    else
        if S.trigger.on == 1
            triggerData = zeros(length(stimData),1);
            nTrigSamps = round(S.trigger.dur*1e-3*NI.Rate);
            triggerData(1:nTrigSamps) = S.trigger.amp;
            S.stim.data(:,1) = triggerData;
            S.stim.transitiondata(:,1) = triggerData;
        end
        S.stim.data(:,2:nD+1) = stimData;
        S.stim.transitiondata(:,2:nD+1) = stimDataTrans;
    end
end
S.stim.buffersize = length(S.stim.data);
S.stim.repsplayed = 0;
set(S.stim.repsplayedbut,'String',['rep ' num2str(S.stim.repsplayed)]);

%-------------------------------------------------------------------------
function [stimData,stimDataTrans] = NIaddBaseStim(stimData,stimDataTrans,zvec_delay,stimPulseSeriesInd);
global S NI

sampperrep = round(NI.Rate*(S.stim.burstrepperiod/1000));
sampperburst = round(NI.Rate*(S.basestim.burstdur/1000));
nzerosamp = sampperrep-sampperburst;

basestimamplitude = S.basestim.amplitude/S.current.onevoltequalsXmilliamps;

tvec = [1/NI.Rate:1/NI.Rate:sampperburst/S.ni.rate];
if S.basestim.dc == 0
    if strcmpi(S.stim.waveformlist(S.basestim.waveformindex),'sine')
        if S.basestim.ampmoddepth == 0
            baseData = basestimamplitude*sin(2*pi*S.basestim.frequency*tvec + S.basestim.phase*2*pi)';
            baseDataTrans = baseData;
        elseif S.basestim.ampmoddepth ~= 0
            baseData = [];
            baseData(:,1) = basestimamplitude*(1-S.basestim.ampmoddepth/100/2)*sin(2*pi*S.basestim.frequency                         *tvec + S.basestim.phase*2*pi)';
            baseData(:,2) = basestimamplitude*(  S.basestim.ampmoddepth/100/2)*sin(2*pi*(S.basestim.frequency+S.basestim.ampmodfreq) *tvec + S.basestim.phase*2*pi + S.basestim.ampmodphase*2*pi)';
            baseDataTrans = baseData;
            if S.stim.sameonallchannels == 1 %send both freqs on 1 channel
                baseData = sum(baseData')';
                baseDataTrans = sum(baseDataTrans')';
            end
        end
    end
elseif S.basestim.dc == 1
    % test DC stim
    baseData = basestimamplitude*ones(size(tvec'));
    baseDataTrans = basestimamplitude*ones(size(tvec'));
end

zvec = zeros(nzerosamp,1);
baseData = [baseData; zvec];
baseDataTrans = [baseDataTrans; zvec];

S.basestim.data = baseData;
S.basestim.transitiondata = baseDataTrans;

if strcmpi(S.basestim.stimcombinemethod,'stop-insert')
    if isempty(stimPulseSeriesInd)
        stimInd = length(zvec_delay)+1:length(zvec_delay)+length(stimData);
        baseData(stimInd) = stimData;
        baseDataTrans(stimInd) = stimDataTrans;
    else
        nDelay = length(zvec_delay);
        for n = 1:length(stimPulseSeriesInd)
            baseData(stimPulseSeriesInd(1,n)+nDelay:stimPulseSeriesInd(2,n)+nDelay) = stimData(stimPulseSeriesInd(1,n):stimPulseSeriesInd(2,n));
            baseDataTrans(stimPulseSeriesInd(1,n)+nDelay:stimPulseSeriesInd(2,n)+nDelay) = stimDataTrans(stimPulseSeriesInd(1,n):stimPulseSeriesInd(2,n));
        end
    end
elseif strcmpi(S.basestim.stimcombinemethod,'add')
    if isempty(stimPulseSeriesInd)
        stimInd = length(zvec_delay)+1:length(zvec_delay)+length(stimData);
        baseData(stimInd) = baseData(stimInd)+stimData;
        baseDataTrans(stimInd) = baseDataTrans(stimInd)+stimDataTrans;
    else
        nDelay = length(zvec_delay);
        for n = 1:length(stimPulseSeriesInd)
            baseData(stimPulseSeriesInd(1,n)+nDelay:stimPulseSeriesInd(2,n)+nDelay) = baseData(stimPulseSeriesInd(1,n)+nDelay:stimPulseSeriesInd(2,n)+nDelay) + stimData(stimPulseSeriesInd(1,n):stimPulseSeriesInd(2,n));
            baseDataTrans(stimPulseSeriesInd(1,n)+nDelay:stimPulseSeriesInd(2,n)+nDelay) = baseDataTrans(stimPulseSeriesInd(1,n)+nDelay:stimPulseSeriesInd(2,n) + nDelay)+stimDataTrans(stimPulseSeriesInd(1,n):stimPulseSeriesInd(2,n));
        end
    end
    
elseif strcmpi(S.basestim.stimcombinemethod,'add-zerocenter')
    if isempty(stimPulseSeriesInd)
        stimInd = length(zvec_delay)+1:length(zvec_delay)+length(stimData);
        nonzeroind = find(stimData~=0);
        baseData(stimInd(nonzeroind)) = stimData(nonzeroind);
        nonzeroind = find(stimDataTrans~=0);
        baseDataTrans(stimInd(nonzeroind)) = stimDataTrans(nonzeroind);
    else
        nDelay = length(zvec_delay);
        for n = 1:length(stimPulseSeriesInd)
            locStimData = stimData(stimPulseSeriesInd(1,n):stimPulseSeriesInd(2,n));
            locStimInd = stimPulseSeriesInd(1,n)+nDelay:stimPulseSeriesInd(2,n)+nDelay;
            nonzeroind = find(locStimData~=0);
            baseData(locStimInd(nonzeroind)) = locStimData(nonzeroind);
            
            locStimDataTrans = stimDataTrans(stimPulseSeriesInd(1,n):stimPulseSeriesInd(2,n));
            baseDataTrans(locStimInd(nonzeroind)) = locStimDataTrans(nonzeroind);
        end
    end
end

if S.basestim.makebasezero
    % make the base stim zero
    newbaseData = zeros(size(baseData));
    newbaseDataTrans = zeros(size(baseDataTrans));    
    if isempty(stimPulseSeriesInd)
        nonzeroind = find(stimData~=0);
        newbaseData(stimInd(nonzeroind)) = baseData(stimInd(nonzeroind));
        nonzeroind = find(stimDataTrans~=0);
        newbaseDataTrans(stimInd(nonzeroind)) = baseDataTrans(stimInd(nonzeroind));
    else
        for n = 1:length(stimPulseSeriesInd)
            locStimData = stimData(stimPulseSeriesInd(1,n):stimPulseSeriesInd(2,n));
            locStimInd = stimPulseSeriesInd(1,n)+nDelay:stimPulseSeriesInd(2,n)+nDelay;
            nonzeroind = find(locStimData~=0);
            newbaseData(locStimInd(nonzeroind)) = baseData(locStimInd(nonzeroind));
            locStimDataTrans = stimDataTrans(stimPulseSeriesInd(1,n):stimPulseSeriesInd(2,n));
            newbaseDataTrans(locStimInd(nonzeroind)) = baseDataTrans(locStimInd(nonzeroind));
        end
    end
    baseData = newbaseData;
    baseDataTrans = newbaseDataTrans;            
end


stimData = baseData;
stimDataTrans = baseDataTrans;

%-------------------------------------------------------------------------
function NImakeSeqStim
global S NI

S.sequence.data.stim1 = S.stim.data;
for n = 2:S.sequence.nseq
    if S.sequence.basestimseq == 0;
        eval(['S.stim.' S.sequence.parametername ' = S.sequence.parametervalues(n);']);
    elseif S.sequence.basestimseq == 1;
        eval(['S.basestim.' S.sequence.parametername ' = S.sequence.parametervalues(n);']);
    end
    NImakeStim;
    eval(['S.sequence.data.stim' num2str(n) '  = S.stim.data;']);
end
S.stim.data = S.sequence.data.stim1;

%-------------------------------------------------------------------------
function [pulse,charge1,charge2,interpeakgap] = NImakepulse;
global S NI

    sampperphase1 = round(NI.Rate*(S.stim.phase1pulsewidth/1e6));
    samppergap = round(NI.Rate*(S.stim.phasegap/1e6));
    sampperphase2 = round(NI.Rate*(S.stim.phase2pulsewidth/1e6));

    if sampperphase1 < 2 | sampperphase2 < 2
        disp('WARNING: Pulse width is too low for this sample rate')
    end

    if sum(strcmpi(S.stim.waveformlist(S.stim.waveformindex),{'pulse','pulse-series'})) 
        part1 =  S.stim.amplitude*(S.stim.phase1amp/100)*ones(sampperphase1,1);
        part2 =  S.stim.amplitude*(S.stim.phase2amp/100)*ones(sampperphase2,1);
    elseif sum(strcmpi(S.stim.waveformlist(S.stim.waveformindex),{'triangle','triangle-series'}))
        widthhalfmax1 = sampperphase1;
        widthhalfmax2 = sampperphase2;

        go = 1;
        while go
            sampperphase1 = sampperphase1+1;
            test = fwhm([1:sampperphase1], triang(sampperphase1));
            if test >= widthhalfmax1;
                go = 0;
            end
        end
        samppergap = samppergap - (sampperphase1-widthhalfmax1)/2;
        go = 1;
        while go
            sampperphase2 = sampperphase2+1;
            test = fwhm([1:sampperphase2], triang(sampperphase2));
            if test >= widthhalfmax1;
                go = 0;
            end
        end
        samppergap = samppergap - (sampperphase2-widthhalfmax2)/2;

        part1 = S.stim.amplitude*(S.stim.phase1amp/100)*triang(sampperphase1);
        part2 = S.stim.amplitude*(S.stim.phase2amp/100)*triang(sampperphase2);
    elseif sum(strcmpi(S.stim.waveformlist(S.stim.waveformindex),{'gaussian','gaussian-series'}))
        widthhalfmax1 = sampperphase1;
        widthhalfmax2 = sampperphase2;

        go = 1;
        while go
            sampperphase1 = sampperphase1+1;
            test = fwhm([1:sampperphase1], gausswinzero(sampperphase1));
            if test >= widthhalfmax1;
                go = 0;
            end
        end
        samppergap = samppergap - (sampperphase1-widthhalfmax1)/2;
        go = 1;
        while go
            sampperphase2 = sampperphase2+1;
            test = fwhm([1:sampperphase2], gausswinzero(sampperphase2));
            if test >= widthhalfmax1;
                go = 0;
            end
        end
        samppergap = samppergap - (sampperphase2-widthhalfmax2)/2;

        part1 = S.stim.amplitude*(S.stim.phase1amp/100)*gausswinzero(sampperphase1);
        part2 = S.stim.amplitude*(S.stim.phase2amp/100)*gausswinzero(sampperphase2);
    end

    charge1 = 1e3 * S.current.onevoltequalsXmilliamps * sum(part1)*(1/S.ni.rate); % charge in microC
    charge2 = 1e3 * S.current.onevoltequalsXmilliamps * sum(part2)*(1/S.ni.rate);
    interpeakgap = (sampperphase1/2 + samppergap + sampperphase2/2)/S.ni.rate * 1e6;

    pulse = [part1; zeros(samppergap,1); part2];
    pulseperiod = 1/S.stim.frequency;
    sampperpulse = round(NI.Rate*pulseperiod);

zerosamps = sampperpulse-length(pulse);
if zerosamps < 0
    disp('WARNING: Pulse width is too long for this stimulation frequency')
end

if S.ampblank.on == 0;
    zvec = zeros(zerosamps,1);
    pulse = [pulse; zvec];
elseif S.ampblank.on == 1;
    S.ampblank.prepulsensamp = ceil(S.ampblank.prepulsetime*S.ni.rate);
    S.ampblank.postpulsensamp = ceil(S.ampblank.postpulsetime*S.ni.rate);
   
    blankpulse = S.ampblank.pulseamp*ones(length(pulse)+S.ampblank.prepulsensamp +S.ampblank.postpulsensamp,1);
    
    zvec1 = zeros(S.ampblank.prepulsensamp,1);
    zvec2 = zeros(zerosamps - S.ampblank.prepulsensamp,1);
    pulse = [zvec1; pulse; zvec2];
    
    zerosamps = sampperpulse-length(blankpulse);
    zvec = zeros(zerosamps,1);
    S.ampblank.singlepulse = [blankpulse; zvec];
end

%-------------------------------------------------------------------------
function NIclose
global NI LH RH AO

if ~isempty(LH)
    delete(LH)
end
if ~isempty(RH)
    delete(RH)
end
if ~isstruct(NI)
    NI.stop
    delete(NI)
end
if ~isempty(AO)
    AO.stop
    delete(AO)
end

%-------------------------------------------------------------------------
function NIupdate(event)
global NI S P AO

% start timer
T = tic;

% make event.Data writable
eventData = event.Data;

% make third channel difference
% for n = 1:3
%     eventData(:,n) = eventData(:,n) - mean(eventData(:,n));
% end
% eventData(:,3) = eventData(:,1) - eventData(:,2);

% save raw recording data
if S.rec.rec == 1
    if exist(S.rec.fullfileName,'file')
        S.data.fid = binappend(S.rec.fullfileName,eventData,[],S.data.fid,S.data.datatype);
    else
        S.data.fid = binwrite(S.rec.fullfileName,S.data.header,eventData,S.data.datatype,0);
    end
    NIstopWatch('update');
end

if S.stim.stim & S.impedance.calc
    normVoltage = S.impedance.voltageconversion * (eventData(:,S.impedance.voltagechind) - mean(eventData(:,S.impedance.voltagechind)));
    normCurrent = S.impedance.currentconversion * (eventData(:,S.impedance.currentchind) - mean(eventData(:,S.impedance.currentchind)));
    S.impedance.val = rms(normVoltage)./rms(normCurrent);
    %disp(S.impedance.val)
end

% make raw plot data
S.rec.rawplotdata = [S.rec.rawplotdata(S.rec.rawplotbuffersize+1:end,:); eventData];

% if S.accel.accel == 1
%     S.rec.procplotdata = NIaccelprocess(S.rec.rawplotdata);
% else
S.rec.procplotdata = S.rec.rawplotdata;
%end

% plot data
if S.rec.showplot == 1
    for n = 1:S.ni.nchin
        S.rec.plotdata(:,n) = S.rec.procplotdata(1:S.rec.plotdownsample:end,n);
        if S.impedance.calc
            if n == S.impedance.voltagechind
                set(S.rec.p(n),'ydata',S.rec.plotdata(:,n)*S.impedance.voltageconversion)
            elseif n == S.impedance.currentchind
                set(S.rec.p(n),'ydata',S.rec.plotdata(:,n)*S.impedance.currentconversion * 1e3)
            else
                set(S.rec.p(n),'ydata',S.rec.plotdata(:,n)/S.amp.gain)
            end
        else
            set(S.rec.p(n),'ydata',S.rec.plotdata(:,n)/S.amp.gain)
        end
    end
end

% listen data
if S.rec.listen
    listenData(:,1) = eventData(:,S.rec.listenchind)/10;
    %listenData(:,2) = eventData(:,S.rec.listenchind)/10;
    S.rec.listendata = [S.rec.listendata; listenData];
end
% process data
if S.proc.proc
    % execute on different core??
    NIprocess(eventData)
    %F = parfeval(P,@NIprocess,0,eventData);
end

% quick sample
if S.rec.quickrec
    %     figure;
    %     %plot(S.rec.timevec,S.rec.plotdata/S.amp.gain)
    %     plot([1/S.ni.rate:1/S.ni.rate:S.ni.buffersize/S.ni.rate]*1e3,eventData/S.amp.gain)
    %     xlabel('Time (ms)')
    %     ylabel('Voltage (V)')
    S.rec.quickrec = 0;
    %    NIquickSave(eventData);
    assignin('base','data',eventData)
    assignin('base','fs',S.ni.rate)
    disp('Data avaialbe in workspace')
end

% display timer
T = toc(T);
ppt = round(100*(T/S.ni.updateperiod));
if ppt>90
    disp(['WARNING: Using ' num2str(ppt) '% of processing time. System is too slow'])
elseif ppt>80
    disp(['Warning: Using ' num2str(ppt) '% of processing time. System may be too slow'])
end

%--------------------------------------------------------------------------
function NIprocess(eventData)
global S

eval(['NIprocess_' S.proc.type '(''data'',eventData);']);

%--------------------------------------------------------------------------
function NIquickRec
global S
S.rec.quickrec = 1;
%--------------------------------------------------------------------------
function NIquickSave(data)
global S

% update name and session number
S.rec.filename = get(S.rec.filenamebut,'String');
S.rec.quickName = [S.data.dir S.rec.filename '.qs'];

% check file name and append _1 if it already exists
S.rec.quickName = NImakeFileName(S.rec.quickName);

save(S.rec.quickName,'S','data')

%--------------------------------------------------------------------------
function NIsaveProcessedData
global S

% update name and session number
S.rec.filename = get(S.rec.filenamebut,'String');
S.proc.filename = [S.data.dir S.rec.filename '.proc'];

% check file name and append _1 if it already exists
S.proc.filename = NImakeFileName(S.proc.filename);

P.stim = S.stim;
P.proc = S.proc;

save(S.proc.filename,'P')

%--------------------------------------------------------------------------
function NIresetProcessedData
global S

eval(['S = NIprocess_' S.proc.type '(''settings'',S);']);
eval(['NIprocess_' S.proc.type '(''plot'');']);

%--------------------------------------------------------------------------
function NIStimFig
global S

% make figure
S.fig.fhdl = findobj('tag','NIStimFigure');
if isempty(S.fig.fhdl)
    S.fig.fhdl = figure;
else
    clf(S.fig.fhdl)
end
set(S.fig.fhdl,'DeleteFcn','NIStim(''close'')')
set(S.fig.fhdl,'menubar','none','toolbar','figure','numbertitle','off')
set(S.fig.fhdl,'name','NIStim')
set(S.fig.fhdl,'tag','NIStimFigure')
xb = 0.12;
yb = 0.05;
S.fig.tgroup = uitabgroup('Parent', S.fig.fhdl);

% add stimulating tab
S.fig.tab1 = uitab('Parent', S.fig.tgroup, 'Title', 'Stimulate');
S.stim.updatebut  = uicontrol('Parent',S.fig.tab1,'Style','push','String','Update','Callback','NIStim(''updateStim'')','units','normalized','position',[0.2 0.05 xb yb]);
if S.stim.stim == 0
    set(S.stim.updatebut,'enable','off')
    S.stim.startbut  = uicontrol('Parent',S.fig.tab1,'Style','push','String','Stimulate','Callback','NIStim(''startStim'')','backgroundcolor',[0 1 0],'units','normalized','position',[0.03 0.05 xb yb]);
elseif S.stim.stim == 1
    set(S.stim.updatebut,'enable','on')
    S.stim.startbut  = uicontrol('Parent',S.fig.tab1,'Style','push','String','Stop Stim','Callback','NIStim(''stopStim'')','backgroundcolor',[1 0 0],'units','normalized','position',[0.03 0.05 xb yb]);
end

if S.current.present == 1
    S.stim.ampbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Amplitude (mA)','units','normalized','position',[0.0 0.9 2*xb yb]);
elseif S.current.present == 0
    S.stim.ampbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Amplitude (V)','units','normalized','position',[0.0 0.9 2*xb yb]);
end
S.stim.ampbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.amplitude),'units','normalized','position',[0.25 0.9 xb yb],'backgroundcolor',[1 1 1]);
S.stim.stimulatorbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Stimulator','units','normalized','position',[0.38 0.89 xb yb]);
S.stim.stimulatorbut = uicontrol('Parent',S.fig.tab1,'Style','popup','String',S.current.namelist,'value',S.current.namevalue,'units','normalized','position',[0.5 0.9 xb yb],'backgroundcolor',[1 1 1],'callback','NIStim(''changeStimulator'')');
S.stim.freqbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Frequency (Hz)','units','normalized','position',[0.0 0.83 2*xb yb]);
S.stim.freqbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.frequency),'units','normalized','position',[0.25 0.83 xb yb],'backgroundcolor',[1 1 1],'callback','NIStim(''frequency'')');
S.stim.phasebuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Phase (cyc)','units','normalized','position',[0.38 0.83 xb yb]);
S.stim.phasebut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.phase),'units','normalized','position',[0.5 0.83 xb yb],'backgroundcolor',[1 1 1]);
S.stim.contbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Continuous','units','normalized','position',[0.0 0.76 2*xb yb]);
S.stim.contbut = uicontrol('Parent',S.fig.tab1,'Style','checkbox','units','normalized','position',[0.25 0.76 xb yb],'value',S.stim.continuous,'callback','NIStim(''continuous'')');
S.stim.rampbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Ramp','units','normalized','position',[0.38 0.76 xb yb]);
S.stim.rampbut = uicontrol('Parent',S.fig.tab1,'Style','checkbox','units','normalized','position',[0.5 0.76 xb yb],'value',S.stim.ramp);
S.stim.burstdurbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Burst dur (ms)','units','normalized','position',[0.0 0.69 2*xb yb]);
S.stim.burstdurbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.burstdur),'units','normalized','position',[0.25 0.69 xb yb],'backgroundcolor',[1 1 1]);
S.stim.burstrepperiodbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Burst rep period (ms)','units','normalized','position',[0.0 0.62 2*xb yb]);
S.stim.burstrepperiodbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.burstrepperiod),'units','normalized','position',[0.25 0.62 xb yb],'backgroundcolor',[1 1 1]);
S.stim.numberrepsbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Reps #','units','normalized','position',[0.38 0.62 xb yb]);
S.stim.numberrepsbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.numberreps),'units','normalized','position',[0.5 0.62 xb yb],'backgroundcolor',[1 1 1]);
S.stim.repsplayedbut = uicontrol('Parent',S.fig.tab1,'Style','text','String',['rep ' num2str(S.stim.repsplayed)],'units','normalized','position',[0.63 0.62 xb yb],'enable','off');
S.stim.ampmoddepthbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','AM depth (%)','units','normalized','position',[0.0 0.55 2*xb yb]);
S.stim.ampmoddepthbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.ampmoddepth),'units','normalized','position',[0.25 0.55 xb yb],'backgroundcolor',[1 1 1]);
S.stim.ampmodfreqbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','AM frequency (Hz)','units','normalized','position',[0.0 0.48 2*xb yb]);
S.stim.ampmodfreqbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.ampmodfreq),'units','normalized','position',[0.25 0.48 xb yb],'backgroundcolor',[1 1 1]);
S.stim.samechanbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Same on all channels','units','normalized','position',[0.375 0.4775 1.8*xb yb]);
S.stim.samechanbut = uicontrol('Parent',S.fig.tab1,'Style','checkbox','value',S.stim.sameonallchannels,'units','normalized','position',[0.595 0.48 xb yb]);
S.stim.waveformbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Wavefrom','units','normalized','position',[0.0 0.41 2*xb yb]);
S.stim.waveformbut = uicontrol('Parent',S.fig.tab1,'Style','popupmenu','String',S.stim.waveformlist,'value',S.stim.waveformindex,'units','normalized','position',[0.25 0.41 xb yb],'backgroundcolor',[1 1 1],'callback','NIStim(''waveform'')');
S.stim.custombut = uicontrol('Parent',S.fig.tab1,'Style','edit','String','','units','normalized','position',[0.4 0.41 xb*1.85 yb],'backgroundcolor',[1 1 1],'callback','NIStim(''cutsomwaveform'')');
S.stim.phase1pwbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Phase 1:   PW (us)','units','normalized','position',[0.0 0.34 2*xb yb]);
S.stim.phase1pwbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.phase1pulsewidth),'units','normalized','position',[0.25 0.34 xb yb],'backgroundcolor',[1 1 1],'callback','NIStim(''pulsewidth'')');
S.stim.phasegapbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Inter-phase gap (us)','units','normalized','position',[0.0 0.27 2*xb yb]);
S.stim.phasegapbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.phasegap),'units','normalized','position',[0.25 0.27 xb yb],'backgroundcolor',[1 1 1],'callback','NIStim(''pulsewidth'')');
S.stim.phase2pwbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Phase 2:   PW (us)','units','normalized','position',[0.0 0.20 2*xb yb]);
S.stim.phase2pwbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.phase2pulsewidth),'units','normalized','position',[0.25 0.20 xb yb],'backgroundcolor',[1 1 1],'callback','NIStim(''pulsewidth'')');
S.stim.phase1ampbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Amp (%)','units','normalized','position',[0.38 0.34 xb yb]);
S.stim.phase1ampbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.phase1amp),'units','normalized','position',[0.5 0.34 xb yb],'backgroundcolor',[1 1 1]);
S.stim.phase2ampbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Amp (%)','units','normalized','position',[0.38 0.20 xb yb]);
S.stim.phase2ampbut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',num2str(S.stim.phase2amp),'units','normalized','position',[0.5 0.20 xb yb],'backgroundcolor',[1 1 1]);

S.stim.phase1Cperphasebuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Charge/phase:','units','normalized','position',[0.62 0.34 1.5*xb yb]);
S.stim.phase1Cperphasebut = uicontrol('Parent',S.fig.tab1,'Style','text','String','','units','normalized','position',[0.8 0.34 xb yb]);
S.stim.phase2Cperphasebuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Charge/phase:','units','normalized','position',[0.62 0.20 1.5*xb yb]);
S.stim.phase2Cperphasebut = uicontrol('Parent',S.fig.tab1,'Style','text','String','','units','normalized','position',[0.8 0.20 xb yb]);
S.stim.interpeakgapbuttxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','Inter peak gap:','units','normalized','position',[0.62 0.27 1.5*xb yb]);
S.stim.interpeakgapbut = uicontrol('Parent',S.fig.tab1,'Style','text','String','','units','normalized','position',[0.8 0.27 xb yb]);

S.stim.stimrecbut  = uicontrol('Parent',S.fig.tab1,'Style','push','String','Stim & Rec','Callback','NIStim(''stimRec'')','units','normalized','position',[0.6 0.05 xb yb]);
S.stim.stimprocbut  = uicontrol('Parent',S.fig.tab1,'Style','push','String','Stim & Proc','Callback','NIStim(''stimProc'')','units','normalized','position',[0.725 0.05 xb yb]);
S.stim.stimrecprocbut  = uicontrol('Parent',S.fig.tab1,'Style','push','String','Stim & R & P','Callback','NIStim(''stimRecProc'')','units','normalized','position',[0.85 0.05 xb yb]);

S.stim.filenamebut = uicontrol('Parent',S.fig.tab1,'Style','edit','String',S.rec.filename,'Callback','NIStim(''stimFileName'')','units','normalized','position',[0.725 0.11 xb*2 yb],'backgroundcolor',[1 1 1]);
S.stim.filenametxt = uicontrol('Parent',S.fig.tab1,'Style','text','String','File name','units','normalized','position',[0.605 0.11 xb yb]);

NIwriteCurrentStimParam
S.stim.current.amplitude = 0; % ramp up on initial switch on

% add recording tab
S.fig.tab2 = uitab('Parent', S.fig.tgroup, 'Title', 'Record');
S.rec.a = axes('parent', S.fig.tab2);
set(S.rec.a,'position',[0.1 0.25 0.85 0.7])
S.rec.p = plot(S.rec.timevec, S.rec.plotdata);
S.rec.legend = legend(S.ni.chilabel);
legend boxoff
ylim(S.rec.ylims)
xlim(S.rec.xlims)
xlabel('Time (s)')
if S.impedance.calc == 1
    ylabel('Voltage (V) | Current (mA)')
else
    ylabel('Voltage (V)')
end

if S.rec.rec == 0
    S.rec.recbut  = uicontrol('Parent',S.fig.tab2,'Style','push','String','Record','Callback','NIStim(''startRec'')','backgroundcolor',[0 1 0],'units','normalized','position',[0.03 0.13 xb yb]);
elseif S.rec.rec == 1
    S.rec.recbut  = uicontrol('Parent',S.fig.tab2,'Style','push','String','Stop Rec','Callback','NIStim(''stopRec'')','backgroundcolor',[1 0 0],'units','normalized','position',[0.03 0.13 xb yb]);
end
S.rec.startstimbut  = uicontrol('Parent',S.fig.tab2,'Style','push','String','Stimulate','Callback','NIStim(''startStim'')','backgroundcolor',[0 1 0],'units','normalized','position',[0.03 0.05 xb yb]);
S.rec.quickbut  = uicontrol('Parent',S.fig.tab2,'Style','push','String','Quick Samp','Callback','NIStim(''quickRec'')','units','normalized','position',[0.17 0.13 xb yb]);
S.rec.filenamebut = uicontrol('Parent',S.fig.tab2,'Style','edit','String',S.rec.filename,'Callback','NIStim(''recFileName'')','units','normalized','position',[0.7 0.05 xb*2 yb],'backgroundcolor',[1 1 1]);
S.rec.filenametxt = uicontrol('Parent',S.fig.tab2,'Style','text','String','File name','units','normalized','position',[0.58 0.05 xb yb]);

S.rec.playbacktxt = uicontrol('Parent',S.fig.tab2,'Style','text','String','Play back mode','units','normalized','position',[0.75 0.11 1.2*xb yb]);
S.rec.playbackbut = uicontrol('Parent',S.fig.tab2,'Style','checkbox','units','normalized','position',[0.91 0.12 xb yb],'value',S.rec.playbackmode,'callback','NIStim(''playbackmode'')');

S.rec.showplottxt = uicontrol('Parent',S.fig.tab2,'Style','text','String','Display on','units','normalized','position',[0.3 0.04 1.2*xb yb]);
S.rec.showplotbut = uicontrol('Parent',S.fig.tab2,'Style','checkbox','units','normalized','position',[0.43 0.05 xb yb],'value',S.rec.showplot,'callback','NIStim(''showrecplot'')');

if S.accel.accel == 1
    S.rec.accelcalibbut = uicontrol('Parent',S.fig.tab2,'Style','push','units','normalized','String','Calibrate','position',[0.17 0.05 xb yb],'callback','NIStim(''accelcalib'')');
end

NIcontinuous
NIwaveform
NIpulsewidth
NIfrequency
NIplaybackmode

% add processing tab
S.fig.tab3 = uitab('Parent', S.fig.tgroup, 'Title', 'Process');
eval(['NIprocess_' S.proc.type '(''plot'');']);
S.proc.proctxt = uicontrol('Parent',S.fig.tab3,'Style','text','String','Processing on','units','normalized','position',[0.75 0.02 1.2*xb yb]);
S.proc.procbut = uicontrol('Parent',S.fig.tab3,'Style','checkbox','units','normalized','position',[0.91 0.03 xb yb],'value',S.proc.proc,'callback','NIStim(''processmode'')');

if strcmp(S.proc.type,'closedLoop')
    % add closed loop button
    S.proc.proctxt = uicontrol('Parent',S.fig.tab3,'Style','text','String','Processing on','units','normalized','position',[0.75 0.02 1.2*xb yb]);
    S.proc.procbut = uicontrol('Parent',S.fig.tab3,'Style','checkbox','units','normalized','position',[0.91 0.03 xb yb],'value',S.proc.proc,'callback','NIStim(''processmode'')');
    S.proc.closedloopbut =uicontrol('Parent',S.fig.tab3,'Style','push','String','Start CL','Callback','NIprocess_closedLoop(''startStim'')','backgroundcolor',[0 1 0],'units','normalized','position',[0.03 0.05 xb yb]);
    S.proc.closetheloop = 0;
end

%--------------------------------------------------------------------------
function NIcontinuous
global S

cval = get(S.stim.contbut,'value');
if cval == 1 % disable burst options
    set(S.stim.burstdurbuttxt,'enable','off');
    set(S.stim.burstdurbut,'enable','off');
    set(S.stim.burstrepperiodbuttxt,'enable','off');
    set(S.stim.burstrepperiodbut,'enable','off');
    set(S.stim.numberrepsbuttxt,'enable','off');
    set(S.stim.numberrepsbut,'enable','off');
    set(S.stim.rampbut,'value',1);
    set(S.stim.rampbut,'enable','on');
    set(S.stim.rampbuttxt,'enable','on');
elseif cval == 0 % enable burst options
    set(S.stim.burstdurbuttxt,'enable','on');
    set(S.stim.burstdurbut,'enable','on');
    set(S.stim.burstrepperiodbuttxt,'enable','on');
    set(S.stim.burstrepperiodbut,'enable','on');
    set(S.stim.numberrepsbuttxt,'enable','on');
    set(S.stim.numberrepsbut,'enable','on');
    set(S.stim.rampbut,'value',0);
    set(S.stim.rampbut,'enable','off');
    set(S.stim.rampbuttxt,'enable','off');
end

%--------------------------------------------------------------------------
function NIwaveform
global S

windex = get(S.stim.waveformbut,'value');
if sum(strcmpi(S.stim.waveformlist(windex),{'pulse','triangle','gaussian'})) % disable pulse options
    set(S.stim.contbut,'enable','on');
    set(S.stim.custombut,'enable','off')
    
    set(S.stim.freqbut,'enable','on');
    set(S.stim.phasebut,'enable','on');
    set(S.stim.ampmoddepthbut,'enable','on');
    set(S.stim.ampmodfreqbut,'enable','on');
    
    set(S.stim.phase1pwbuttxt,'enable','on');
    set(S.stim.phase1pwbut,'enable','on');
    set(S.stim.phasegapbuttxt,'enable','on');
    set(S.stim.phasegapbut,'enable','on');
    set(S.stim.phase2pwbuttxt,'enable','on');
    set(S.stim.phase2pwbut,'enable','on');
    set(S.stim.phase1ampbuttxt,'enable','on');
    set(S.stim.phase1ampbut,'enable','on');
    set(S.stim.phase2ampbuttxt,'enable','on');
    set(S.stim.phase2ampbut,'enable','on');
    
    set(S.stim.phase1Cperphasebuttxt,'enable','on');
    set(S.stim.phase1Cperphasebut,'enable','on');
    set(S.stim.phase2Cperphasebuttxt,'enable','on');
    set(S.stim.phase2Cperphasebut,'enable','on');
    set(S.stim.interpeakgapbuttxt,'enable','on');
    set(S.stim.interpeakgapbut,'enable','on');
elseif strcmpi(S.stim.waveformlist(windex),'sine')
    set(S.stim.contbut,'enable','on');
    
    set(S.stim.custombut,'enable','off')
    
    set(S.stim.freqbut,'enable','on');
    set(S.stim.phasebut,'enable','on');
    set(S.stim.ampmoddepthbut,'enable','on');
    set(S.stim.ampmodfreqbut,'enable','on');
    
    set(S.stim.phase1pwbuttxt,'enable','off');
    set(S.stim.phase1pwbut,'enable','off');
    set(S.stim.phasegapbuttxt,'enable','off');
    set(S.stim.phasegapbut,'enable','off');
    set(S.stim.phase2pwbuttxt,'enable','off');
    set(S.stim.phase2pwbut,'enable','off');
    set(S.stim.phase1ampbuttxt,'enable','off');
    set(S.stim.phase1ampbut,'enable','off');
    set(S.stim.phase2ampbuttxt,'enable','off');
    set(S.stim.phase2ampbut,'enable','off');
    
    
    set(S.stim.phase1Cperphasebuttxt,'enable','off');
    set(S.stim.phase1Cperphasebut,'enable','off');
    set(S.stim.phase2Cperphasebuttxt,'enable','off');
    set(S.stim.phase2Cperphasebut,'enable','off');
    set(S.stim.interpeakgapbuttxt,'enable','off');
    set(S.stim.interpeakgapbut,'enable','off');
    
elseif strcmpi(S.stim.waveformlist(windex),'custom')
    set(S.stim.contbut,'enable','off','value',0);
    NIcontinuous
    
    set(S.stim.custombut,'enable','on')
    
    set(S.stim.freqbut,'enable','off');
    set(S.stim.phasebut,'enable','off');
    set(S.stim.ampmoddepthbut,'enable','off');
    set(S.stim.ampmodfreqbut,'enable','off');
    
    set(S.stim.phase1pwbuttxt,'enable','off');
    set(S.stim.phase1pwbut,'enable','off');
    set(S.stim.phasegapbuttxt,'enable','off');
    set(S.stim.phasegapbut,'enable','off');
    set(S.stim.phase2pwbuttxt,'enable','off');
    set(S.stim.phase2pwbut,'enable','off');
    set(S.stim.phase1ampbuttxt,'enable','off');
    set(S.stim.phase1ampbut,'enable','off');
    set(S.stim.phase2ampbuttxt,'enable','off');
    set(S.stim.phase2ampbut,'enable','off');
    
    set(S.stim.phase1Cperphasebuttxt,'enable','off');
    set(S.stim.phase1Cperphasebut,'enable','off');
    set(S.stim.phase2Cperphasebuttxt,'enable','off');
    set(S.stim.phase2Cperphasebut,'enable','off');
    set(S.stim.interpeakgapbuttxt,'enable','off');
    set(S.stim.interpeakgapbut,'enable','off');
    
end
NIpulsewidth

%--------------------------------------------------------------------------
function NIcustomwaveform
global S NI

[S.stim.customfilename, S.stim.stimdir] = uigetfile(S.stim.stimdir, 'Select a stimulus file');
S.stim.stimdir = [S.stim.stimdir '\'];
set(S.stim.custombut,'String',S.stim.customfilename)

load([S.stim.stimdir S.stim.customfilename])
if A.fs ~= NI.Rate
    S.stim.customdata = resample(A.stim.data,1,A.fs/NI.Rate);
else
    S.stim.customdata = A.stim.data;
end
nCD = length(S.stim.customdata);
set(S.stim.burstdurbut,'String',(nCD/NI.Rate)*1e3);
set(S.stim.burstrepperiodbut,'String',(nCD/NI.Rate)*1e3);
if isfield(A,'sequence')
    set(S.stim.numberrepsbut,'String',num2str(A.stim.numberreps));
    S.stim.randomizesequence = A.stim.randomizesequence;
end

%S.stim.ramp = 0;
%S.stim.continuous = 0;
S.basestim.stim = 0;
NIparseStimGUI

if isfield(A,'sequence') 
    S.sequence.on = A.sequence.on;
    S.sequence.thisseq = 0;
    S.sequence.loopthisseq = 1;
    S.sequence.nseq = A.sequence.nseq;
    S.sequence.parametername = A.sequence.parametername;
    S.sequence.basestimseq = 0;
    S.sequence.parametervalues = A.sequence.parametervalues;
    S.sequence.seq = [1:S.sequence.nseq];
    if S.stim.numberreps == Inf
        S.sequence.seqIndex =  S.sequence.seq;
    else
        if S.stim.randomizesequence == 0
            S.sequence.seqIndex =  repmat(S.sequence.seq,1,S.stim.numberreps);
        elseif S.stim.randomizesequence == 1
            S.sequence.seqIndex = [];
            for r = 1:S.stim.numberreps
                rInd = randperm(S.sequence.nseq);
                S.sequence.seqIndex =  [S.sequence.seqIndex S.sequence.seq(rInd)];
            end
        end
    end
    for n = 1:S.sequence.nseq
        eval(['S.sequence.data.stim' num2str(n) '  = squeeze(S.stim.customdata(n,:,:));']);
    end
end

if S.ni.nchout<2
    S.stim.customdata = S.stim.customdata(:,:,2);
end

if isfield(A,'trigger')
    S.trigger.on = A.trigger.on;
end
%sound(S.stim.customdata,NI.Rate)

%--------------------------------------------------------------------------
function NIpulsewidth
global S NI

NIparseStimGUI

sampperphase1 = round(NI.Rate*(S.stim.phase1pulsewidth/1e6));
samppergap = round(NI.Rate*(S.stim.phasegap/1e6));
sampperphase2 = round(NI.Rate*(S.stim.phase2pulsewidth/1e6));

S.stim.phase1pulsewidth = (sampperphase1/NI.Rate)*1e6;
S.stim.phase2pulsewidth =  (sampperphase2/NI.Rate)*1e6;
S.stim.phasegap =  (samppergap/NI.Rate)*1e6;

set(S.stim.phase1pwbut,'string',num2str(S.stim.phase1pulsewidth));
set(S.stim.phase2pwbut,'string',num2str(S.stim.phase2pulsewidth));
set(S.stim.phasegapbut,'string',num2str(S.stim.phasegap));

% calculate and display charge per phase
NIparseStimGUI
if sum(strcmpi(S.stim.waveformlist(S.stim.waveformindex),{'pulse','triangle','gaussian'}))
    [pulse,charge1,charge2,interpeakgap] = NImakepulse;
    set(S.stim.phase1Cperphasebut,'String',[num2str(charge1) ' uC'])
    set(S.stim.phase2Cperphasebut,'String',[num2str(charge2) ' uC'])
    set(S.stim.interpeakgapbut,'String',[num2str(interpeakgap) ' uS'])
    pf = findobj('tag','pulseFig');
    if isempty(pf)
        figure('tag','pulseFig')
    else
        figure(pf)
    end
    plot([1/S.ni.rate:1/S.ni.rate:length(pulse)/S.ni.rate]*1e6,pulse,'.-')
    hold on
    assignin('base','pulse',pulse)
    assignin('base','time',[1/S.ni.rate:1/S.ni.rate:length(pulse)/S.ni.rate]*1e6)
    
    if S.ampblank.on == 1;
        plot([1/S.ni.rate:1/S.ni.rate:length(pulse)/S.ni.rate]*1e6,S.ampblank.singlepulse,'r')
    end
    xlabel('Time (\muS)')
    xlim([0 2*(S.stim.phase1pulsewidth + S.stim.phase2pulsewidth + S.stim.phasegap)])
    ylim([min(pulse)*1.1 max(pulse)*1.1])
    figure(S.fig.fhdl)
end
%--------------------------------------------------------------------------
function NIfrequency
global S NI

NIparseStimGUI
if S.stim.frequency > NI.Rate/2
    set(S.stim.freqbut, 'string', num2str(round(NI.Rate/2)))
end

% adjust stim freq to give exact value base on number of samples per period
S.stim.frequency = round(NI.Rate)/round(NI.Rate/S.stim.frequency);
set(S.stim.freqbut,'String',num2str(S.stim.frequency));
        

%--------------------------------------------------------------------------
function NIplaybackmode
global S NI

S.rec.playback = get(S.rec.playbackbut,'value');
if S.rec.playback == 1
    set(S.rec.recbut,'string','Play','callback','NIStim(''startPlay'')')
    set(S.stim.startbut,'enable','off')
    set(S.rec.startstimbut,'enable','off')
    if ~isstruct(NI)
        if NI.IsRunning
            NI.stop;
        end
    end
else
    set(S.rec.recbut,'string','Record','callback','NIStim(''startRec'')')
    set(S.stim.startbut,'enable','on')
    set(S.rec.startstimbut,'enable','on')
    if ~isstruct(NI)
        if ~NI.IsRunning
            NI.startBackground();
        end
    end
end

%--------------------------------------------------------------------------
function NIaccelcalib
global S NI

disp('Hold accelerometer still with the x-axis pointing up')
a = input('Press ENTER when ready');
D = NIrecord10sec('x_calib');
x_still = D.data(:,S.accel.chinds);

disp('Hold accelerometer still with the y-axis pointing up')
a = input('Press ENTER when ready');
D = NIrecord10sec('y_calib');
y_still = D.data(:,S.accel.chinds);

disp('Hold accelerometer still with the z-axis pointing up')
a = input('Press ENTER when ready');
D = NIrecord10sec('z_calib');
z_still = D.data(:,S.accel.chinds);

[Calibration] = calibration(x_still, y_still, z_still, S.ni.rate);
S.accel.callibration.g0 = [Calibration.x_0g Calibration.y_0g Calibration.z_0g];
S.accel.callibration.g1 = [Calibration.g_x Calibration.g_y Calibration.g_z];

g0 = S.accel.callibration.g0;
g1 = S.accel.callibration.g1;
save([S.path 'AccelCalib\CalibrationData_' date],'g0','g1')

S.accel.donecalib = 1;

%--------------------------------------------------------------------------
function D = NIrecord10sec(filename)
global S

% add date string with no spaces or colons
filename = [filename strrep(strrep(datestr(now),':','-'),' ','-')];
disp('Recording 10 seconds of data. Please wait...')
set(S.rec.filenamebut,'String',filename)
NIStim('startRec')
pause(10)
NIStim('stopRec')
disp('Done')
D = binread([S.data.dir filename '.bin']);

%--------------------------------------------------------------------------
function data = NIaccelprocess(data)
global S NI

for n = 1:3
    data(:,S.accel.chinds(n)) = 10*data(:,S.accel.chinds(n)) - 13;
end

% if S.accel.donewcalib == 1;
%     S.accel.calib = mean(data(:,S.accel.chinds));
%     S.accel.lasteventdata = [0 0 0];
%     S.accel.donewcalib = 0;
% end
%
% for n = 1:3
%     data(:,n) = data(:,n) - S.accel.calib(n);
% end
%
% data(:,S.accel.chinds) = filter(S.accel.filterb,S.accel.filtera,data(:,S.accel.chinds));
% data(:,S.accel.chinds) = cumtrapz(data(:,1:3));
% data(:,S.accel.chinds) = filter(S.accel.filterb,S.accel.filtera,data(:,S.accel.chinds));
% data(:,S.accel.chinds) = cumtrapz(data(:,S.accel.chinds));
%
% data(:,S.accel.chinds) =  data(:,S.accel.chinds)/1000000;
% for n = 1:3
%     data(:,S.accel.chinds(n)) = S.accel.lasteventdata(n) + data(:,S.accel.chinds(n));
% end
% S.accel.lasteventdata = data(S.rec.rawplotbuffersize+1,S.accel.chinds);
%
% %data(:,S.accel.chinds) = filter(S.accel.filterb,S.accel.filtera,data(:,S.accel.chinds));

%--------------------------------------------------------------------------
function NIprocessmode
global S NI

S.proc.proc = get(S.proc.procbut,'value');

%--------------------------------------------------------------------------
function NIshowrecplot
global S NI

S.rec.showplot = get(S.rec.showplotbut,'value');
if S.rec.showplot == 1
    set(S.rec.a,'visible','on')
    set(S.rec.p,'visible','on')
elseif S.rec.showplot == 0
    set(S.rec.a,'visible','off')
    set(S.rec.p,'visible','off')
end

%--------------------------------------------------------------------------
function NIstimRecProc
global S NI
NIstartRec
S.proc.proc = 1;
set(S.proc.procbut,'value',S.proc.proc);
NIstartStim
%--------------------------------------------------------------------------
function NIstimRec
global S NI
NIstartRec
NIstartStim

%--------------------------------------------------------------------------
function NIstimProc
global S NI
S.proc.proc = 1;
set(S.proc.procbut,'value',S.proc.proc);
NIstartStim

%--------------------------------------------------------------------------
function NIstimFileName
global S
set(S.rec.filenamebut,'String',get(S.stim.filenamebut,'String'))

%--------------------------------------------------------------------------
function NIrecFileName
global S
set(S.stim.filenamebut,'String',get(S.rec.filenamebut,'String'))

%--------------------------------------------------------------------------
function NIlistenOn(chind)
global S AO NI

queueOutputData(AO,rand(AO.Rate,1));
queueOutputData(AO,rand(AO.Rate,1));
S.rec.listen = 1;
if ~isempty(chind)
    S.rec.listenchind = chind;
end
AO.startBackground();
AO.IsRunning

%--------------------------------------------------------------------------
function   NIlistenOff
global S AO

S.rec.listen = 0;
AO.stop
S.rec.listendata = [];

%--------------------------------------------------------------------------
function   NIqueueListen(event)

queueData = S.rec.listendata(1:NI.rate);
size(queueData)
queueOutputData(NI,queueData);
S.rec.listendata = S.rec.listendata(NI.rate+1:end);

%--------------------------------------------------------------------------
function AO = NIlistenCreatObj
global NI

AO = daq.createSession('directsound');
try
    addAudioOutputChannel(AO,'Audio3',1);
catch
    addAudioOutputChannel(AO,'Audio1',1);
end
AOL = addlistener(AO,'DataRequired',@(src,event) NIStim(src,event,3));
AO.NotifyWhenScansQueuedBelow = NI.Rate;
AO.UseStandardSampleRates = false;
AO.Rate = NI.Rate;
AO.IsContinuous = 1;

%--------------------------------------------------------------------------
function NIchangeStimulator
global S

S.current.namevalue = get(S.stim.stimulatorbut,'value');
S.current.name = S.current.namelist(S.current.namevalue);
S.current.onevoltequalsXmilliamps = S.current.onevoltequalsXmilliampslist(S.current.namevalue);
set(S.stim.ampbuttxt,'String',['Amplitude (' S.current.namelistunits{S.current.namevalue} ')']);

%--------------------------------------------------------------------------
function NIstopWatch(command)
global S

switch command
    case 'start'
        tic
    case 'update'
        a=toc; sec = rem(a,60); min = (a-sec)/60;
        if floor(sec)<10
            timeStr = [num2str(min) ':0' num2str(floor(sec))];
        else
            timeStr = [num2str(min) ':' num2str(floor(sec))];
        end
        set(S.rec.recbut,'string',timeStr)
    case 'stop'
        timeStr = ['0:00'];
end



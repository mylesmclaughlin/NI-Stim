function NIclosedloop(command,input1,input2)

if isobject(command)
    if input2 == 2
        command = 'queueStim';
    end
end

switch command
    case 'init'
        NICLinit;
    case 'processdata'
        NICLdata
    case 'queueStim'
        NICLqueueStim(input1)
    case 'startStim'
        NICLstartStim
    case 'stopStim'
        NICLstopStim
end


%--------------------------------------------------------------------------
function NICLinit
global CL S OUT

device = daq.getDevices;
devInd = [];
for n = 1:length(device)
    if strcmp(device(n).Description,'National Instruments USB-6343')
        devInd = n;
        break
    end
end
disp(['Connected to ' device(devInd).Description])
CL.ni.connected = 1;
CL.ni.devname = device(devInd).ID;
CL.ni.description = device(devInd).Description;
CL.ni.chout = [0 1]; %[0 1 2 3]; % 0 out is always trigger
CL.ni.nchout = length(CL.ni.chout);
CL.ni.rate = 1000;%200e3; %200e3; %

OUT = daq.createSession('ni');
chOut = addAnalogOutputChannel(OUT,CL.ni.devname,CL.ni.chout, 'Voltage');
RH = addlistener(OUT,'DataRequired',@(src,event) NIclosedloop(src,event,2));
OUT.IsContinuous = 1;
OUT.NotifyWhenScansQueuedBelow = S.ni.buffersize;

triggerData = zeros(S.ni.buffersize,1);
nTrigSamps = round(S.trigger.dur*1e-3*OUT.Rate);
triggerData(1:nTrigSamps) = S.trigger.amp;
S.stim.data(:,1) = triggerData;
S.stim.data(:,2) = zeros(S.ni.buffersize,1);

queueOutputData(OUT,S.stim.data);
queueOutputData(OUT,S.stim.data);

%--------------------------------------------------------------------------
function NICLdata
global S
% test phase extraction
accel = sqrt(S.rec.procplotdata(:,4).^2 + S.rec.procplotdata(:,5).^2 + S.rec.procplotdata(:,6).^2);
accel = accel-mean(accel);
tremorphase = angle(hilbert(S.rec.procplotdata(:,5) - mean(S.rec.procplotdata(:,5))));
S.rec.procplotdata(:,2) = tremorphase;
%tvec = [1/NI.Rate:1/NI.Rate:length(eventData(:,2))/S.ni.rate];
closedloopstim = S.stim.amplitude*cos(tremorphase')';
S.rec.procplotdata(:,3) = closedloopstim;

S.stim.data(:,2) = closedloopstim(end-S.ni.buffersize+1:end);
%--------------------------------------------------------------------------
function NICLqueueStim(input1)
global OUT S
queueOutputData(OUT,S.stim.data);

%--------------------------------------------------------------------------
function NICLstartStim
global OUT 
OUT.startBackground();
disp('Starting closed-loop stimulation')

%--------------------------------------------------------------------------
function NICLstopStim
global OUT 
stop(OUT)
disp('Stopping closed-loop stimulation')

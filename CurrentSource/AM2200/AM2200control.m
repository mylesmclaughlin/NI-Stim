function D = AM2200control(command)

if nargin==0
    command = 'connect';
end

switch command
    case 'connect'
        D = AM2200connect;
    case 'autoZero'
        D = AM2200autoZero;
    case 'enableOutput'
        D = AM2200enableOutput;
    case 'disableOutput'
        D = AM2200disableOutput;
end

%--------------------------------------------------------------------------
function D = AM2200connect
global X B

device = daq.getDevices;
% 
% if strcmp(device(1).Description,'National Instruments USB-6216 (BNC)')
%     devInd = 1;
% elseif strcmp(device(2).Description,'National Instruments USB-6216 (BNC)')
%     devInd = 2;
% end
devInd = strmatch('National Instruments USB-6216 (BNC)',{device.Description});

X.devname = device(devInd).ID;
X.description = device(devInd).Description;

B = daq.createSession('ni');
X.chD = addDigitalChannel(B,X.devname,'Port0/Line1:2', 'OutputOnly');
disp(['Connected to ' device(devInd).Description ' for AM 2200 current source control'])
D = 1;

%--------------------------------------------------------------------------
function   D = AM2200enableOutput;
global B

outputSingleScan(B,[1,1])
D = 1;
disp('AM 2200 output enabled')
%--------------------------------------------------------------------------
function D = AM2200disableOutput;
global B

outputSingleScan(B,[0,0])
D = 1;
disp('AM 2200 output disabled')


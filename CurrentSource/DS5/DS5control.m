function D = DS5control(command)

if nargin==0
    command = 'connect';
end

switch command
    case 'connect'
        D = DS5connect;
    case 'autoZero'
        D = DS5autoZero;
    case 'enableOutput'
        D = DS5enableOutput;
    case 'disableOutput'
        D = DS5disableOutput;
end

%--------------------------------------------------------------------------
function D = DS5connect
global A IDevice

% iControl Object
A = actxserver('DS5.Control');

% iDevice Object
B = actxserver('DS5.Application');
IDevices = B.Devices;
if IDevices.count == 0
    disp('WARNING: DS5 not connected')
    D = 0;
else
    IDevice = Item(IDevices,'0');
    D = 1;
end

%--------------------------------------------------------------------------
function D = DS5autoZero
global A IDevice

IControl = A.Clone;
IControl.AutoZero = 'True';
set(IDevice,'Control',IControl);
D = 1;

%--------------------------------------------------------------------------
function D = DS5enableOutput
global A IDevice

IControl = A.Clone;
IControl.OutputEnable = 'True';
set(IDevice,'Control',IControl);
D = 1;

%--------------------------------------------------------------------------
function D = DS5disableOutput
global A IDevice

IControl = A.Clone;
IControl.OutputEnable = 'False';
set(IDevice,'Control',IControl);
D = 1;

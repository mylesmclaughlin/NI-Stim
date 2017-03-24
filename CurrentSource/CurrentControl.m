function D = CurrentControl(command,name)


if ischar(name)
    hit = strfind(name,'DS5');
    if ~isempty(hit)
        name = 'Digitimer DS5';
    end
elseif iscell(name)
    name = name{1};
end

switch command
    case 'connect'
        D = CCconnect(name);
    case 'autoZero'
        D = CCautoZero(name);
    case 'enableOutput'
        D = CCenableOutput(name);
    case 'disableOutput'
        D = CCdisableOutput(name);
end

%--------------------------------------------------------------------------
function   D = CCconnect(name);

name
switch name
    case 'Digitimer DS5'
        D = DS5control('connect');
    case 'AM 2200'
        D = AM2200control('connect');
end
%--------------------------------------------------------------------------
function  D = CCautoZero(name);

switch name
    case 'Digitimer DS5'
        D = DS5control('autoZero');
    case 'AM 2200'
        D = 0;
        % manual zero
end

%--------------------------------------------------------------------------
function D = CCenableOutput(name);

switch name
    case 'Digitimer DS5'
        D = DS5control('enableOutput');
    case 'AM 2200'
        D = AM2200control('enableOutput');
end
%--------------------------------------------------------------------------
function D = CCdisableOutput(name);

switch name
    case 'Digitimer DS5'
        D = DS5control('disableOutput');
    case 'AM 2200'
        D = AM2200control('disableOutput');
end

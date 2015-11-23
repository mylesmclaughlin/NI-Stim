function OUT = NIprocess_trigSpike(command,input)

switch command
    case 'settings'
        OUT = NIPsettings(input);
    case 'plot'
        NIPplot
    case 'data'
        NIPdata(input);
    case 'clear'
        NIPclear
    case 'triglevel'
        NIPtriglevel
end

%--------------------------------------------------------------------------
function S = NIPsettings(S)

S.proc.rawdata = [];
S.proc.trigchind = 2;
S.proc.spikechind = 2;
S.proc.spikeconversion = 1;
S.proc.triglevel = 1;
S.proc.nspikesinbuffer = 10;
S.proc.bufferposition = 0;
S.proc.bufferfilled = 0;
S.proc.spikerate = 0;

S.proc.timewindow = [-100 500];
S.proc.trigdeadtime = 200e-3;
S.proc.trigdeadsamp = round(S.proc.trigdeadtime * S.ni.rate);

S.proc.sampvec = [round(S.proc.timewindow(1)/1e3 * S.ni.rate)+1 : round(S.proc.timewindow(2)/1e3 * S.ni.rate)];
S.proc.timevec = 1e3*(S.proc.sampvec/S.ni.rate); % [-S.proc.timewindow(1)/1e3+1/S.ni.rate:1/S.ni.rate:S.proc.timewindow(2)/1e3];
S.proc.bufferdata = zeros(length(S.proc.timevec),S.proc.nspikesinbuffer);
S.proc.bufferrelativetime = zeros(S.proc.nspikesinbuffer,1);

S.proc.avgdata = zeros(length(S.proc.timevec),1);
S.proc.txlims = S.proc.timewindow;

%--------------------------------------------------------------------------
function NIPplot
global S NI

if isfield(S.proc,'a1')
    delete(S.proc.a1)
end
S.proc.a1 = axes('parent', S.fig.tab3);
set(S.proc.a1,'position',[0.1 0.2 0.85 0.75])
S.proc.p1 = plot(S.proc.timevec,S.proc.bufferdata,'k');
hold on
S.proc.p2 = plot(S.proc.timevec,S.proc.avgdata,'r');
xlabel('Time (ms)')
ylabel('Amp (V)')
xlim(S.proc.txlims)
S.proc.lh = legend([num2str(S.proc.spikerate) ' spikes/sec']);
legend boxoff

xb = 0.12;
yb = 0.05;
S.proc.clearbut  = uicontrol('Parent',S.fig.tab3,'Style','push','String','Clear','Callback','NIprocess_trigSpike(''clear'')','units','normalized','position',[0.03 0.05 xb yb]);

S.proc.trigleveltxt  = uicontrol('Parent',S.fig.tab3,'Style','text','String','Trigger Level','units','normalized','position',[0.2 0.05 xb yb]);
S.proc.triglevelbut  = uicontrol('Parent',S.fig.tab3,'Style','edit','String',num2str(S.proc.triglevel),'backgroundcolor',[1 1 1],'Callback','NIprocess_trigSpike(''triglevel'')','units','normalized','position',[0.33 0.05 xb yb]);

%--------------------------------------------------------------------------
function NIPdata(eventData)
global S NI

S.proc.rawdata = [S.proc.rawdata; eventData];
%disp(['Size of raw data = ' num2str(size(S.proc.rawdata))]);

%trigInd = find(S.proc.rawdata(1:end-1,S.ni.nchin)<S.proc.triglevel & S.proc.rawdata(2:end,S.ni.nchin)>S.proc.triglevel);
trigInd = find(S.proc.rawdata(1:end-1,S.proc.trigchind)<S.proc.triglevel & S.proc.rawdata(2:end,S.proc.trigchind)>S.proc.triglevel);
%disp(['Number of triggers = ' num2str(length(trigInd))]);

lastSampUsed = -S.proc.trigdeadsamp;
for n = 1:length(trigInd)
    % apply dead time
    if trigInd(n)>lastSampUsed+S.proc.trigdeadsamp 
        % check that sample fits in raw data sample
        if trigInd(n)+S.proc.sampvec(1)>0 & trigInd(n)+S.proc.sampvec(end)<length(S.proc.rawdata)
            % we have a valid spike! 
            
            % store relative time base on last spike time
            if S.proc.bufferposition == 0
                S.proc.bufferrelativetime(1) = 0;
            elseif S.proc.bufferposition > S.proc.nspikesinbuffer
                S.proc.bufferrelativetime(1) = S.proc.bufferrelativetime(S.proc.bufferposition) + (trigInd(n)-lastSampUsed)/S.ni.rate;
            else
                S.proc.bufferrelativetime(S.proc.bufferposition + 1) = S.proc.bufferrelativetime(S.proc.bufferposition) + (trigInd(n)-lastSampUsed)/S.ni.rate;
            end
            
            %... now adjust buffer position
            S.proc.bufferposition = S.proc.bufferposition +1;
            if S.proc.bufferposition > S.proc.nspikesinbuffer
                S.proc.bufferposition = 1;
                S.proc.bufferfilled = 1;
            end
            
            % store spike in buffer
            S.proc.bufferdata(:,S.proc.bufferposition) = S.proc.rawdata(S.proc.sampvec+trigInd(n),S.proc.spikechind) * S.proc.spikeconversion;
            % update last spike used
            lastSampUsed = trigInd(n);
        end
    end
end
if lastSampUsed<0
    lastSampUsed = 0;
end

if S.proc.bufferfilled
    S.proc.avgdata(:,1) = mean(S.proc.bufferdata,2);
    relativetime = sort(S.proc.bufferrelativetime);
    S.proc.spikerate = mean(diff(relativetime));
elseif S.proc.bufferposition == 0
    S.proc.avgdata(:,1) = S.proc.avgdata(:,1);
    S.proc.spikerate =  S.proc.spikerate;
else
    S.proc.avgdata(:,1) = mean(S.proc.bufferdata(:,1:S.proc.bufferposition),2);
    relativetime = sort(S.proc.bufferrelativetime(1:S.proc.bufferposition));
    S.proc.spikerate = mean(diff(relativetime));
end

%disp(['lastSampUsed = ' num2str(lastSampUsed)]);
if isempty(trigInd) % no triggers in this data - throw away
    S.proc.rawdata = S.proc.rawdata(end+S.proc.sampvec(1)-1:end,:);
else
    S.proc.rawdata = S.proc.rawdata(lastSampUsed+1:end,:);
end

% update plot
for n = 1:S.proc.nspikesinbuffer
    set(S.proc.p1(n),'ydata',S.proc.bufferdata(:,n))
end
set(S.proc.p2,'ydata',S.proc.avgdata)
legend(S.proc.a1,[num2str(1/S.proc.spikerate) ' spikes/sec'])
legend boxoff

%--------------------------------------------------------------------------
function NIPclear
global S

S.proc.bufferposition = 0;
S.proc.bufferfilled = 0;
S.proc.spikerate = 0;
S.proc.bufferdata = zeros(length(S.proc.timevec),S.proc.nspikesinbuffer);
S.proc.bufferrelativetime = zeros(S.proc.nspikesinbuffer,1);
S.proc.avgdata = zeros(length(S.proc.timevec),1);

%--------------------------------------------------------------------------
function NIPtriglevel
global S

S.proc.triglevel = str2num(get(S.proc.triglevelbut,'string'));

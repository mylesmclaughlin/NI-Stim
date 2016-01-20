function [A,D,M] = NImakeMacroAvg(macroFileName,dataPath)

if nargin<2
    dataPath = 'C:\Users\Public\Data\NIData\';
end

load([dataPath macroFileName],'-mat')
seqparametername = 'amplitude';
%seqparametername = 'phasegap';

for n = 1:length(M)
    M(n).macro.seqname =  strrep(M(n).macro.seqname,'-','');
    M(n).macro.seqname =  strrep(M(n).macro.seqname,'.','');
end

trigChInd = 3;
trigLevel = 0.5;
D.macroList = {};
D.seqparametername = seqparametername;
eval(['D.' D.seqparametername ' = M(1).stim.' D.seqparametername ';']);
eval(['D.seqparametervalues = M(1).stim.' D.seqparametername ';']);
D.nSeq = length(D.seqparametervalues);

if strcmpi(D.seqparametername,'amplitude')
    D.frequency = M(1).stim.frequency;
    D.phasegap = M(1).stim.phasegap;
elseif strcmpi(D.seqparametername,'frequency')
    D.amplitude = M(1).stim.amplitude;
    D.phasegap = M(1).stim.phasegap;
elseif strcmpi(D.seqparametername,'phasegap')
    D.amplitude = M(1).stim.amplitude;
    D.frequency = M(1).stim.frequency;
end

extraTime = 100;
D.sampWin = [-extraTime M(1).stim.burstrepperiod]; 

for n = 1:length(M)
    if M(n).macro.done
        [allData,fs,sampVec] = loadBinFileMacroData(M(n),dataPath,seqparametername,D.sampWin,n);
    end
    D.fs = fs;
    if isempty(D.macroList)
        D.macroList{1} = M(n).macro.seqname;
        if strcmpi(D.seqparametername,'amplitude')
            D.frequency(1) = M(n).stim.frequency;
        end
        thisSeqInd = 1;
        newSeq = 1;
    else 
        hitList = strcmp(D.macroList,M(n).macro.seqname);
        if sum(hitList) == 0
            D.macroList{end+1} = M(n).macro.seqname;
            if strcmpi(D.seqparametername,'amplitude')
                D.frequency(end+1) = M(n).stim.frequency;
            end
            thisSeqInd = length(D.macroList);
            newSeq = 1;
        else
            thisSeqInd = find(hitList==1);
            newSeq = 0;
        end
    end
    
    if newSeq == 1
        eval(['D.allData' D.macroList{thisSeqInd} '(1,:,:,:) = allData;']);
        eval(['D.nRep'  D.macroList{thisSeqInd} ' = 1;']);
    else
        eval(['nRep = D.nRep'  D.macroList{thisSeqInd} ' + 1;']); 
        eval(['D.nRep'  D.macroList{thisSeqInd} ' = nRep;']);
        eval(['D.allData' D.macroList{thisSeqInd} '(' num2str(nRep) ',:,:,:) = allData;']);
    end
    
end

if exist('nRep')==0
    nRep=1;
end

D.nReps = nRep;

    
D.nMacro = length(D.macroList);
D.allDataMatrixInfo = {'Reps','Sequence','Time','Channels'};
D.avgDataMatrixInfo = {'Sequence','Time','Channels'};
D.tvec = sampVec/D.fs;

% average and std data
for m = 1:D.nMacro
    eval(['D.avgData' D.macroList{m} ' = squeeze(mean(D.allData' D.macroList{m} '(:,:,:,:),1));']);
    eval(['D.stdData' D.macroList{m} ' = squeeze(std(D.allData' D.macroList{m} '(:,:,:,:),0,1));']);
    eval(['A(' num2str(m)  ').avgData = D.avgData' D.macroList{m} ';']);
    eval(['A(' num2str(m)  ').stdData = D.stdData' D.macroList{m} ';' ]);
    eval(['A(' num2str(m)  ').allData = D.allData' D.macroList{m} ';' ]);
    eval(['A(' num2str(m)  ').macro = D.macroList{' num2str(m) '}  ;' ]);
    eval(['A(' num2str(m)  ').tvec = D.tvec;']);
    eval(['A(' num2str(m)  ').fs = D.fs;']);
    eval(['A(' num2str(m)  ').nSeq = D.nSeq;']);
    eval(['A(' num2str(m)  ').seqparametervalues = D.seqparametervalues;']);
    eval(['A(' num2str(m)  ').seqparametername = D.seqparametername;']);
    eval(['A(' num2str(m)  ').amplitude = D.amplitude;']);
    eval(['A(' num2str(m)  ').frequency = D.frequency(m);']);
    eval(['A(' num2str(m)  ').phasegap = D.phasegap;']);
    for n = 1:D.nSeq
        eval(['A(' num2str(m)  ').allData' num2str(n) ' = squeeze(D.allData' D.macroList{m} '(:,n,:,:));' ]);
    end
end

sortList = {'BP','IPG','PM','TRI','GAUS'};
sum(strcmpi(A(1).macro,sortList))
if sum(strcmpi(A(1).macro,sortList))
    for n = 1:length(sortList)
        hitInd = find(strcmpi(sortList{n},D.macroList) == 1);
        newA(n) = A(hitInd);
    end
    A = newA;
end

if strfind(A(1).macro,'IPG');
    for n = 1:length(A)
        ipg(n) = str2num(A(n).macro(4:end));
    end
    [s,sind] = sort(ipg);
    A = A(sind);
end

%--------------------------------------------------------------------------
function [allData,fs,sampVec] = loadBinFileMacroData(M,dataPath,seqParam,sampWin,N)

trigChInd = 3;
trigLevel = 0.7;
eval(['nSeq = length(M.stim.' seqParam ');']);
nRep = M.macro.nreps;
nMacro = 5;
if isfield(M.macro,'allinonebinfile')
    M.macro.allinonebinfile = 0;
end

if M.macro.allinonebinfile == 0
    D = binread([dataPath M.macro.localfilename '.bin']);
    sInd = strfind(D.header,'Fs=');
    eInd = strfind(D.header,'Filter=');
    fs = str2num(D.header(sInd+3:eInd-3));
    sampVec = [round(sampWin(1)/1e3*fs):round(sampWin(2)/1e3*fs)];
    trigInd = find(D.data(1:end-1,trigChInd)<trigLevel & D.data(2:end,trigChInd)>trigLevel);
    nTrig = length(trigInd);
    if nTrig<nSeq
        disp('Could not find all triggers')
    elseif nTrig>nSeq
        disp('Found too many triggers')
    end
    
    j=0;
    for n = 1:nTrig
        if trigInd(n)+sampVec(1)>0 & trigInd(n)+sampVec(end)<length(D.data)
            j = j+1;
            if isfield(M,'sequence')
                allData(M.sequence.seqIndex(j),:,:) = D.data(trigInd(n)+sampVec,:);
            else
                allData(j,:,:) = D.data(trigInd(n)+sampVec,:);
            end
        elseif trigInd(n)+sampVec(1)<=0
            disp('WARNING: Trigger error in data processing')
            return
        end
    end
else
    D = binread([dataPath M.macro.basefilename '.bin']);
    sInd = strfind(D.header,'Fs=');
    eInd = strfind(D.header,'Filter=');
    fs = str2num(D.header(sInd+3:eInd-3));
    sampVec = [round(sampWin(1)/1e3*fs):round(sampWin(2)/1e3*fs)];
    trigInd = find(D.data(1:end-1,trigChInd)<trigLevel & D.data(2:end,trigChInd)>trigLevel);
    nTrig = length(trigInd);
    expectedTrig = nSeq * nRep * nMacro;
    if nTrig<expectedTrig
        disp('Could not find all triggers')
        disp(['Expected: '  num2str(expectedTrig) ' but found: ' num2str(nTrig)])
    elseif nTrig>expectedTrig
        disp('Found too many triggers')
        disp(['Expected: '  num2str(expectedTrig) ' but found: ' num2str(nTrig)])
    end
   
    trigRange = [1:nSeq] + nSeq*(N-1);
    j=0;
    for n = trigRange
        if trigInd(n)+sampVec(1)>0 & trigInd(n)+sampVec(end)<length(D.data)
            j = j+1;
            if isfield(M,'sequence')
                allData(M.sequence.seqIndex(j),:,:) = D.data(trigInd(n)+sampVec,:);
            else
                allData(j,:,:) = D.data(trigInd(n)+sampVec,:);
            end
        elseif trigInd(n)+sampVec(1)<=0
            disp('WARNING: Trigger error in data processing')
            return
        end
    end
end
function [A,D] = NImakeMacroAvg(macroFileName)

dataPath = 'C:\Users\Public\Data\NIData\';
load(macroFileName,'-mat')
seqparametername = 'amplitude';

trigChInd = 3;
trigLevel = 0.5;
D.macroList = {};
D.seqparametername = seqparametername;
eval(['D.' D.seqparametername ' = M(1).stim.' D.seqparametername ';']);
eval(['D.seqparametervalues = M(1).stim.' D.seqparametername ';']);
D.nSeq = length(D.seqparametervalues);

if strcmpi(D.seqparametername,'amplitude')
    D.frequency = M(1).stim.frequency;
elseif strcmpi(D.seqparametername,'frequency')
    D.amplitude = M(1).stim.amplitude;
end

extraTime = 100;
D.sampWin = [-extraTime M(1).stim.burstrepperiod]; 

for n = 1:length(M)
    if M(n).macro.done
        [allData,fs,sampVec] = loadBinFileMacroData(M(n),dataPath,seqparametername,D.sampWin);
    end
    D.fs = fs;
    if isempty(D.macroList)
        D.macroList{1} = M(n).macro.seqname;
        thisSeqInd = 1;
        newSeq = 1;
    else 
        hitList = strcmp(D.macroList,M(n).macro.seqname);
        if sum(hitList) == 0
            D.macroList{end+1} = M(n).macro.seqname;
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
    eval(['A(' num2str(m)  ').frequency = D.frequency;']);
    
    for n = 1:D.nSeq
        eval(['A(' num2str(m)  ').allData' num2str(n) ' = squeeze(D.allData' D.macroList{m} '(:,n,:,:));' ]);
    end
end


%--------------------------------------------------------------------------
function [allData,fs,sampVec] = loadBinFileMacroData(M,dataPath,seqParam,sampWin)

trigChInd = 3;
trigLevel = 0.5;
eval(['nSeq = length(M.stim.' seqParam ');']);

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
            allData(j,:,:) = D.data(trigInd(n)+sampVec,:);
        elseif trigInd(n)+sampVec(1)<=0
            disp('WARNING: Trigger error in data processing')
            return
        end
    end
        
end
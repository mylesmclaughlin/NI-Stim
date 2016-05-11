function newA = NIseries2macro(A,macroFile)

if nargin<2
    macroFile = A.macroFile;
end

load(macroFile,'-mat')
repperiod = P.stim.series.burstrepperiod;
amplitude = P.stim.series.amplitude;
nrep = length(amplitude);


%load bin file
sampWin = [1/A.fs*1000 repperiod]; 
nshift = round(sampWin(2)/1e3*A.fs);
for n = 1:length(A.seqparametervalues)
    avgData = [];
    stdData = [];
    allData = [];
    sampVec = [round(sampWin(1)/1e3*A.fs):round(sampWin(2)/1e3*A.fs)];
    
    for m = 1:nrep
        locsampvec = sampVec+(m-1)*nshift;
        avgData(m,:,:) = A.avgData(n,locsampvec,:);
        stdData(m,:,:) = A.stdData(n,locsampvec,:);
        eval(['allData(:,m,:,:) = A.allData' num2str(n) '(:,locsampvec,:);']);
    end
    newA(n).avgData = avgData;
    newA(n).stdData = stdData;
    newA(n).allData = allData;
    newA(n).fs = A.fs;
    newA(n).tvec = A.tvec(sampVec);
    newA(n).seqparametervalues = amplitude;
    newA(n).seqparametername = 'amplitude';
    newA(n).amplitude = amplitude;
    newA(n).phase = P.basestim.frequency * (A.seqparametervalues(n)/1e3);
    newA(n).nSeq = nrep;
    newA(n).frequency = P.stim.frequency;
    newA(n).macro = num2str(newA(n).phase);
    newA(n).sampWin = A.sampWin;
    for m = 1:nrep
        eval(['newA(n).allData' num2str(m) ' = squeeze(newA(n).allData(:,' num2str(m) ',:,:));' ]);
    end
end


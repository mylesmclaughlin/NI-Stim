function A = NItwoSeqParamPhasePlot(A,freq)

if nargin<2
    freq = 2;
end

% get series info
load(A.macroFile,'-mat')
try
    repperiod = P.stim.series.burstrepperiod;
    amplitude = P.stim.series.amplitude;
    nrep = length(amplitude);
    sampWin = [1/A.fs*1000 repperiod];
    nshift = round(sampWin(2)/1e3*A.fs);
catch
    nrep = 1;
    sampWin = A.sampWin;
    nshift = 1;
end

seqValues1 = unique(A.seqparametervalues);
seqValues2 = unique(A.seqparametervalues2);
nSeq1 = length(seqValues1);
nSeq2 = length(seqValues2);

A.phase = freq *seqValues1/1000;
[b,a] = butter(2,[30 1000]/(A.fs/2),'bandpass');

for m = 1:nSeq2
    for n = 1:nSeq1
        o = (m-1)*nSeq1 + n;
        sampVec = [round(sampWin(1)/1e3*A.fs):round(sampWin(2)/1e3*A.fs)];
        sampVec = [1:length(A.tvec)];
        for p = 1:nrep
            locsampvec = sampVec+(p-1)*nshift;
            
            [val,ind] = max(abs(A.avgDisp(o,locsampvec,4)));
            A.respAmpX(n,m,p) = val;
            A.respSTDX(n,m,p) = A.stdDisp(o,locsampvec(ind),4);
            
            [val,ind] = max(abs(A.avgDisp(o,locsampvec,5)));
            A.respAmpY(n,m,p) = val;
            A.respSTDY(n,m,p) = A.stdDisp(o,locsampvec(ind),5);
            
            [val,ind] = max(abs(A.avgDisp(o,locsampvec,6)));
            A.respAmpZ(n,m,p) = val;
            A.respSTDZ(n,m,p) = A.stdDisp(o,locsampvec(ind),6);
            
            eval(['MEP = A.allData' num2str(o) '(:,locsampvec,2);']);
            [nAvg,dum] = size(MEP);
            for q = 1:nAvg
                MEP(q,:) = filter(b,a,MEP(q,:));
                MEParea(q) = sum(abs(MEP(q,:)));
            end
            A.MEParea(n,m,p) = mean(MEParea);
            A.MEPstd(n,m,p) = std(MEParea);
        end
    end
end
accelInd = [4 5 6];
copo = {'b','r','g','m','k','y','c','b','r','g','m','k','y','c'};

for p = 1:nrep
    figure('position',[159 295 1304 420])
    subplot(1,4,1)
    hold on
    for m = 1:nSeq2
        errorbar(A.phase,A.respAmpX(:,m,p)+1000*m,A.respSTDX(:,m,p),copo{m})
    end
    ylabel('Peak Limb Displacement')
    xlabel('Phase')
    xlim([min(A.phase) max(A.phase)])
    
    subplot(1,4,2)
    hold on
    for m = 1:nSeq2
        errorbar(A.phase,A.respAmpY(:,m,p)+1000*m,A.respSTDY(:,m,p),copo{m})
    end
    ylabel('Peak Limb Displacement')
    xlabel('Phase')
    xlim([min(A.phase) max(A.phase)])
    
    subplot(1,4,3)
    hold on
    for m = 1:nSeq2
        errorbar(A.phase,A.respAmpZ(:,m,p)+1000*m,A.respSTDZ(:,m,p),copo{m})
    end
    ylabel('Peak Limb Displacement')
    xlabel('Phase')
    xlim([min(A.phase) max(A.phase)])
    
    subplot(1,4,4)
    hold on
    for m = 1:nSeq2
        errorbar(A.phase,A.MEParea(:,m,p)+10*m,A.MEPstd(:,m,p),copo{m})
    end
    ylabel('Peak Limb Displacement')
    xlabel('Phase')
    xlim([min(A.phase) max(A.phase)])
   
end
function A = NIrepPlot(A)

nRep = size(A.allData1,1);
nCol = 3;
nSeq = A.nSeq
curChInd = 1;
volChInd = 2;
accelChInd = [4 5 6]; 
xl = [A.tvec(1)*1e3 A.repDur*3];
[dum,eInd] = min(abs(A.tvec-A.repDur/1e3));
[dum,sInd] = min(abs(A.tvec-0));    

for i = 1:nSeq;
    figure('name',(['Seq ' num2str(i)]))
    nP = 0;
    ycur = [];
    yvol = [];
    ymov = [];
    impedanceAmp = [];
    for n = 1:nRep
        
        % current
        nP = nP+1;
        subplot(nRep,nCol,nP)
        eval(['ydata = squeeze(A.allData' num2str(i) '(n,:,curChInd));'])
        ydata = ydata - mean(ydata);
        plot(A.tvec*1e3,ydata,'r')
        xlim(xl)
        set(gca, 'XTick', []);
        ycur(n,:) = [min(ydata) max(ydata)];
        
        % voltage
        nP = nP+1;
        subplot(nRep,nCol,nP)
        eval(['ydata = squeeze(A.allData' num2str(i) '(n,:,volChInd));'])
        plot(A.tvec*1e3,ydata,'b')
        xlim(xl)
        set(gca, 'XTick', []);
        yvol(n,:) = [min(ydata) max(ydata)];
        voldata = ydata(sInd:eInd);
        voldata = voldata-mean(voldata);
        impedanceAmp(n)  = max(mfft(voldata,A.fs,[],0))/(A.amplitude(i)*1e-3);
        
        % movement
        nP = nP+1;
        subplot(nRep,nCol,nP)
        eval(['ydata = squeeze(A.allData' num2str(i) '(n,:,accelChInd));'])
        ydata = sqrt(ydata(:,1).^2 + ydata(:,2).^2 + ydata(:,3).^2);
        plot(A.tvec*1e3,ydata,'g')
        xlim(xl)
        set(gca, 'XTick', []);
        ymov(n,:) = [min(ydata) max(ydata)];
    end
    ycur = [min(ycur(:,1)) max(ycur(:,2))];
    yvol = [min(yvol(:,1)) max(yvol(:,2))];
    ymov = [min(ymov(:,1)) max(ymov(:,2))];
    
    nP = 0;
    for n = 1:nRep
        % current
        nP = nP+1;
        subplot(nRep,nCol,nP)
        ylim(ycur)
        
        % voltage
        nP = nP+1;
        subplot(nRep,nCol,nP)
        ylim(yvol)
        text(xl(2)*.6,yvol(2)*.6,[num2str(round(impedanceAmp(n))) '\ohm'])
     
        % movement
        nP = nP+1;
        subplot(nRep,nCol,nP)
        ylim(ymov)
        
    end
end
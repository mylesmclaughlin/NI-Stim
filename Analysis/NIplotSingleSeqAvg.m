function A = NIplotSingleSeqAvg(A)

accelInd = [4 5 6];
bandpass = [3 500];
filterorder = 2;
[filterb filtera] = butter(filterorder,bandpass/(A.fs/2),'bandpass');
A.avgDataDisp = zeros(size(A.avgData));
zInd = find(A.tvec<0);
trigArtInd = find(A.tvec<15e-3);
stimStartInd = trigArtInd(end)+1;
axisString = {'X','Y','Z'};
for m = 1:3
    for n = 1:A.nSeq
        data = squeeze(A.avgData(n,:,accelInd(m)));
        calib = mean(data(zInd));
        data = data - calib;
        dataAccelFilt = filtfilt(filterb,filtera,data(stimStartInd:end));
        data = cumtrapz(dataAccelFilt);
        data = filtfilt(filterb,filtera,data);
        A.avgDataDisp(n,stimStartInd:end,accelInd(m)) = cumtrapz(data);
        ymin(n) = min(A.avgDataDisp(n,:,accelInd(m)));
        ymax(n) = max(A.avgDataDisp(n,:,accelInd(m)));
    end
    
    yl(1) = min(ymin);  
    yl(2) = max(ymax);
    figure('Name',[axisString{m} '-axis'])
    for n = 1:A.nSeq
        subplot(3,4,n)
        plot(A.tvec,yl(2)*0.5*normalize(A.avgData(n,:,1)),'color',[0.7 0.7 0.7])
        hold on
        plot(A.tvec,yl(2)*0.5*normalize(A.avgData(n,:,accelInd(m))),'k','linewidth',0.2)
        plot(A.tvec,(A.avgDataDisp(n,:,accelInd(m))),'g','linewidth',2)
        ylim(yl);
        xlim([min(A.tvec) max(A.tvec)])
        if n==1
            legend(['Current=' num2str(A.amplitude(n)) 'mA'] ,'Acceleration','Displacement')
            legend boxoff
        else
            legend([num2str(A.amplitude(n)) 'mA'])
            legend boxoff
        end
    end
end

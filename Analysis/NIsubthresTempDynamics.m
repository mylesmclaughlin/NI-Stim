function A = NIsubthresTempDynamics(A)

nreps = size(A(1).allData,2);
nseq = size(A(1).allData,1);

accelInd = [4 5 6];
bandpass = [3 500];
filterorder = 2;
[filterb filtera] = butter(filterorder,bandpass/(A.fs/2),'bandpass');
A.allDataDisp = zeros(size(A.allData));
zInd = find(A.tvec<0);
trigArtInd = find(A.tvec<15e-3);
stimStartInd = trigArtInd(end)+1;
axisString = {'X','Y','Z'};
copo = {'b','r','g','m','k','y','c','b','r','g','m','k','y','c'};

for m = 1:length(accelInd)
    for n = 1:nreps
        for seqInd = 1:nseq;
            data = squeeze(A.allData(seqInd,n,:,accelInd(m)));
            calib = mean(data(zInd));
            data = data - calib;
            dataAccelFilt = filtfilt(filterb,filtera,data(stimStartInd:end));
            data = cumtrapz(dataAccelFilt);
            data = filtfilt(filterb,filtera,data);
            A.allDataDisp(seqInd,n,stimStartInd:end,accelInd(m)) = cumtrapz(data);
            ymin(n,seqInd) = min(A.allDataDisp(seqInd,n,:,accelInd(m)));
            ymax(n,seqInd) = max(A.allDataDisp(seqInd,n,:,accelInd(m)));
            peakDisp(n,seqInd) = max(abs(A.allDataDisp(seqInd,n,:,accelInd(m))));
        end
    end
    yl(1) = min(ymin(:));
    yl(2) = max(ymax(:));
    figure('Name',[axisString{m} '-axis'])
    for n = 1:nreps
        for seqInd = 1:nseq
            subplot(1,nreps,n)
            plot(A.tvec,squeeze(A.allDataDisp(seqInd,n,:,accelInd(m))),copo{seqInd})
            hold on
            %plot(A.tvec,yl(2)*0.5*normalize(squeeze(A.allData(seqInd,n,:,accelInd(m)))),'k','linewidth',0.2)
            ylim(yl);
            xlim([min(A.tvec) max(A.tvec)])
        end 
    end
    
    if m == 1
        eval(['A.peakDispX = peakDisp;'])
    elseif m == 2
        eval(['A.peakDispY = peakDisp;'])
    elseif m == 3
        eval(['A.peakDispZ = peakDisp;'])
    end
    
end
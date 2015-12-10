function A = NIgetAccelAmp(A,accelChInd,doplot);

if nargin<3
    doplot = 0;
end

axisstr = {'X','Y','Z'};
if doplot == 1
    figure
    xmin = A.tvec(1);
    xmax = A.tvec(end)/4;
end

for i = 1:length(accelChInd)
    for n = 1:A.nSeq
        avgData = squeeze(A.avgData(n,:,accelChInd(i)));
        avgData = avgData-mean(avgData);
        stdData = squeeze(A.stdData(n,:,accelChInd(i)));
        [A.posAmp(i),A.posInd(i)] = max(avgData);
        [A.negAmp(i),A.negInd(i)] = min(avgData);
        A.posSTD(i) = stdData(A.posInd(i));
        A.negSTD(i) = stdData(A.negInd(i));
        
        eval(['A.respAmp' axisstr{i} '(' num2str(n) ') = A.posAmp(i)-A.negAmp(i);']);
        eval(['A.respSTD' axisstr{i} '(' num2str(n) ') = A.posSTD(i) + A.negSTD(i);']);
        
        if doplot == 1
            subplot(1,3,i)
            plot(A.tvec,avgData,'b','linewidth',2)
            hold on
            plot(A.tvec,avgData+stdData,'k')
            plot(A.tvec,avgData-stdData,'k')
            plot(A.tvec(A.posInd(i)),A.posAmp(i),'ro')
            plot(A.tvec(A.negInd(i)),A.negAmp(i),'ro')
            xlim([xmin xmax])
            title([axisstr{i} ' axis'])
        end
    end
end

if doplot == 1
    [ymax,ind] = max(A.posAmp);
    ymax = ymax + A.posSTD(ind);
    [ymin,ind] = min(A.negAmp);
    ymin = ymin - A.negSTD(ind);
    for i = 1:length(accelChInd)
        subplot(1,3,i)
        ylim([ymin ymax])
    end
end
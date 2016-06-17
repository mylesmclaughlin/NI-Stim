

% filename = {'RatTACSPhase-31-May-2016-13-35-36.bin'...
%             'RatTACSPhase-31-May-2016-13-45-05.bin',...
%             'RatTACSPhase-31-May-2016-13-54-05.bin',...
%             'RatTACSPhase-31-May-2016-14-04-22.bin'...
%             'RatTACSPhase-31-May-2016-14-14-49.bin'...
%             };
filename = {'RatTACSPhase-31-May-2016-14-36-10.bin'...
            'RatTACSPhase-31-May-2016-14-45-57.bin',...
            'RatTACSPhase-31-May-2016-14-53-10.bin',...
            'RatTACSPhase-31-May-2016-15-00-24.bin'...
            'RatTACSPhase-31-May-2016-15-09-52.bin'...
            };
        
    
A = NImakeSeqAvg(filename{1});        
for m = 1:9
    eval(['P.allDisp' num2str(m) ' = zeros(size(A.allDisp' num2str(m) '));']);
end

for  n = 1:length(filename)
    A = NImakeSeqAvg(filename{n});
    for m = 1:9
        eval(['P.allDisp' num2str(m) ' = P.allDisp' num2str(m) ' + A.allDisp' num2str(m) ';'])
    end
end
% 
% for m = 1:9
%     eval(['P.allDisp' num2str(m) ' = P.allDisp' num2str(m) '/length(filename);'])
%     eval(['P.avgData(' num2str(n) ',:,:) = squeeze(mean(P.allDisp' num2str(n) '(:,:,:)));']);
%     eval(['P.stdData(' num2str(n) ',:,:) = squeeze(std(P.allDisp' num2str(n) '(:,:,:)));']);
% end


accelInd = [4 5 6];
for n = 1:15
    for m = 1:9
        eval(['[val,ind] = max(abs(P.allDisp' num2str(m) '(n,:,accelInd(1))));']);
        P.ampX(n,m) = val;
        
        eval(['[val,ind] = max(abs(P.allDisp' num2str(m) '(n,:,accelInd(2))));']);
        P.ampY(n,m) = val;
        
        eval(['[val,ind] = max(abs(P.allDisp' num2str(m) '(n,:,accelInd(3))));']);
        P.ampZ(n,m) = val;
    end
end

figure
for n = 1:15
    plot(P.ampX(n,:)+n*5000)
    hold on
end

figure
for n = 1:15
    plot(P.ampY(n,:)+n*5000)
    hold on
end

figure
for n = 1:15
    plot(P.ampZ(n,:)+n*5000)
    hold on
end



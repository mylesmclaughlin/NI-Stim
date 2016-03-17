function NInetworkPlot(A)
nSeq = 12;

figure; hold on
for n = 1:nSeq/2
    plot(A.tvec*1e3,A.avgData(n,:,4),'b')
end

for n = nSeq/2+1:nSeq
    plot(A.tvec*1e3,A.avgData(n,:,4),'r')
end

figure; hold on
for n = 1:nSeq/2
    plot(A.tvec*1e3,A.avgDisp(n,:,4),'b')
end

for n = nSeq/2+1:nSeq
    plot(A.tvec*1e3,A.avgDisp(n,:,4),'r')
end

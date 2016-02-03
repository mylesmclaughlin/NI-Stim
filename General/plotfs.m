function p = plotfs(fs,data)

tvec = [1:length(data)]/fs;
p = plot(tvec*1e3,data);
xlabel('Time (ms)')
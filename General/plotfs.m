function p = plotfs(fs,data,color)

if nargin<3
    color = 'b';
end
tvec = [1:length(data)]/fs;
p = plot(tvec,data,color);
xlabel('Time (ms)')
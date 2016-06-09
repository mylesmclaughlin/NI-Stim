function A = NIStim_MakeChirp(dur,f0,f1,t1)

pathname = 'C:\Users\Public\MATLAB\NI-Stim\Stimuli\';
A.fs = 200e3;
phi = -90;
method = 'linear';
t = [1/A.fs:1/A.fs:dur];
A.stim.data = chirp(t,f0,t1,f1,'linear',phi)

% charge balance
offset = sum(A.stim.data)/length(A.stim.data);
A.stim.data = A.stim.data-offset;
A.stim.data(2,:) = A.stim.data(1,:);
A.stim.data(1,:) = zeros(size(A.stim.data(2,:)));
A.stim.data = A.stim.data';

filename = ['Chirp_type=' method '_dur=' num2str(dur) '_f0=' num2str(f0) '_f1=' num2str(f1) '_t1=' num2str(t1) '.mat'];
disp(['Saving ' pathname filename])
save([pathname filename],'A')
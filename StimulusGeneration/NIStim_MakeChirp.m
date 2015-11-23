function A = NIStim_MakeChirp(dur,f0,f1,t1)

pathname = 'C:\Users\u0043883\Google Drive\Work\MATLAB\NI-Stim\Stimuli\';
A.fs = 100e3;
phi = -90;
method = 'linear';
t = [1/A.fs:1/A.fs:dur];
A.stim = chirp(t,f0,t1,f1,'linear',phi)

% charge balance
offset = sum(A.stim)/length(A.stim);
A.stim = A.stim-offset;

filename = ['Chirp_type=' method '_dur=' num2str(dur) '_f0=' num2str(f0) '_f1=' num2str(f1) '_t1=' num2str(t1) '.mat'];
disp(['Saving ' pathname filename])
save([pathname filename],'A')
function data = accel3disp(data,fs);

% filter
S.ni.rate = fs;
S.accel.bandpass = [1 50];
S.accel.filterorder = 2;
[S.accel.filterb S.accel.filtera] = butter(S.accel.filterorder,S.accel.bandpass/(S.ni.rate/2),'bandpass');


S.accel.calib = mean(data(:,1:3));

for n = 1:3
    data(:,n) = data(:,n) - S.accel.calib(n);
end

data(:,1:3) = filter(S.accel.filterb,S.accel.filtera,data(:,1:3));
data(:,1:3) = cumtrapz(data(:,1:3));
data(:,1:3) = filter(S.accel.filterb,S.accel.filtera,data(:,1:3));
data(:,1:3) = cumtrapz(data(:,1:3));



clc;clear;close all;

%read the raw data .wave file here
[Y,FS] = audioread('Serialdata_final_test4.wav');
%constants
c = 3E8; %(m/s) speed of light

%radar parameters
Tp = 100E-3; %(s) pulse time variable in hardware and code
N = Tp*FS; %# of samples per pulse variable in hardware and code
fstart = 2260E6; %(Hz) LFM start frequency for example based on our VCO
fstop = 2590E6; %(Hz) LFM stop frequency for example based on our VCO
%fstart = 2402E6; %(Hz) LFM start frequency for ISM band
%fstop = 2495E6; %(Hz) LFM stop frequency for ISM band
BW = fstop-fstart; %(Hz) transmti bandwidth
f = linspace(fstart, fstop, N/2); %instantaneous transmit frequency

%range resolution
rr = c/(2*BW);
max_range = rr*N/2;
%[Y,FS] = audioread('if22.wav');
%the input appears to be inverted
trig = -3*(Y(:,1));
plot(trig)
s = -1*((Y(:,2)));

%parse the data here by triggering off rising edge of sync pulse
count = 0;
thresh = 0;
start = (trig > thresh);
for ii = 100:(size(start,1)-N)
    if start(ii) == 1 & mean(start(ii-11:ii-1)) == 0
        %start2(ii) = 1;
        count = count + 1;
        sif(count,:) = s(ii:ii+N-1);
%         disp('00000000000000')
        time(count) = ii*1/FS;
    end
end
% plot(sif)
% return
%check to see if triggering works
figure(1),plot(trig/1);
hold on;
plot(s,'r');
hold off;
grid on;

%subtract the average
ave = mean(sif,1);
for ii = 1:size(sif,1);
    sif(ii,:) = sif(ii,:) - ave;
end

zpad = 8*N/2;

%RTI plot
figure(10);
v = dbv(ifft(sif,zpad,2));
S = v(:,1:size(v,2)/2);
m = max(max(v));
imagesc(linspace(0,max_range,zpad),time,(1*S)-m,[-35, 0]);
colorbar;
ylabel('time (s)');
xlabel('range (m)');
title('RTI without clutter rejection');
xlim([0 50]);



sif2 = sif(2:size(sif,1),:)-sif(1:size(sif,1)-1,:);
v = ifft(sif2,zpad,2);
S=v;
R = linspace(0,max_range,zpad);
S = dbv(S(:,1:size(v,2)/2));
m = max(max(S));
figure(21);
imagesc(R,time,(1*S)-m,[-35, 0]);
colorbar;
ylabel('time (s)');
xlabel('range (m)');
title('RTI with 2-pulse cancelor clutter rejection');
xlim([0 50]);

save('sif2','sif2')
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
        time(count) = ii*1/FS;
    end
end
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
% %% Compressive sensing tbp<winsize/2;winsize is max windowing of data, i.e., I take 700 out of 763
% % change K and m_cs for optimum values
% %%%%%%%%%%%%%%% Changing params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K=5; m_cs=150;winsize=500;tbp=120;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% CS_basis='Gaussian'; %'';%'Fourier';%
% CS_slct='DGSR';
% % addpath(genpath('Wavelet'))
% sif00=sif;
% max_sif=max(max(abs(sif)));
% sif00=(sif00/max_sif);  %scaling for positive values for log
% % for ii=1:size(sif00,2)
% % [A,D] = dwt(sif00(1:700,ii)','sym4');
% % x = idwt(A,D,'sym4');    
% %     sif00(1:700,ii)=x';
% % % [a,d] = dwt_lifting1D(sif00(:,ii)','LeGall_5x3');
% % % return
% % end
% for ii=1:size(sif00,1)
%     % g_x0=(255/log(256)).*log(sif00(:,ii)'+ones(size(sif00(:,ii)'))); % multiplicative noise to additive noise conversion
% %     g_x0=fft();
% g_x0=sif00(ii,:);
% %     [g_x1R]= process_CS(real(g_x0),K,m_cs,CS_basis,CS_slct,winsize,tbp);
% %     [g_x1I]= process_CS(imag(g_x0),K,m_cs,CS_basis,CS_slct,winsize,tbp);
%     [g_x1]= process_CS((g_x0),K,m_cs,CS_basis,CS_slct,winsize,tbp);
% 
% %     g_x1(isnan(g_x1))=0; g_x1(isinf(g_x1))=0;
%     %  g_x2=(256.^(g_x1/255))-1; %back from log domain to time domain
% %     g_x2=real(ifft(g_x1R+g_x1I));
% %     g_x2(isnan(g_x2))=0; g_x2(isinf(g_x2))=0;
%     
% %     snr00(ii)=10*log10(mean(sif00(:,ii))/std(sif00(:,ii)));
%     sif00(ii,1:length(g_x1))=g_x1;
%     %   sif00=(sif00*max_sif);
% %     snr01(ii)=10*log10(mean(sif00(:,ii))/std(sif00(:,ii)));
%     ii
% end
% 
% num_iter = 30;
%   delta_t = 1/7;
%   kappa = 40;
%   option = 2;
%   ad = anisodiff2D(sif,num_iter,delta_t,kappa,option);
% %   figure, subplot 121, imshow(sif,[]), subplot 122, imshow(ad,[])





% snr00(isnan(snr00))=0; snr00(isinf(snr00))=0;
% SNR
% [mean(abs(snr01))+(10*log(max_sif)) mean(abs(snr00))]
% figure(11);
% v = dbv(ifft(ad,zpad,2));
% S = v(:,1:size(v,2)/2);
% m = max(max(v));
% imagesc(linspace(0,max_range,zpad),time,(1*S)-m,[-35, 0]);
% colorbar;
% ylabel('time (s)');
% xlabel('range (m)');
% title('RTI without clutter rejection with CS');
% xlim([0 50]);
%%
%2 pulse cancelor RTI plot

% sif2 = ad(2:size(ad,1),:)-ad(1:size(ad,1)-1,:);
% v = ifft(sif2,zpad,2);
% S=v;
% R = linspace(0,max_range,zpad);
% S = dbv(S(:,1:size(v,2)/2));
% m = max(max(S));
% figure(20);
% imagesc(R,time,(1*S)-m,[-25, 0]);
% colorbar;
% ylabel('time (s)');
% xlabel('range (m)');
% title('RTI with 2-pulse cancelor clutter rejection with CS');
% xlim([0 50]);


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

% [Ychk,FS] = audioread('sa2.wav');
% csvwrite('sif.csv',sif)
% M = csvread('data.csv');
% imshow(M)

%%
t = 0 : 1/FS : Tp-( 1/FS ) ;
hlength=100;

%N=7.
N=15;
hlength=hlength+1-rem(hlength,2);
h = tftb_window(hlength);
h=h/norm(h);

for ii=1:count
tfr = GLCT((sif(ii,:))',N,FS,hlength);
%Reconstruction of Signal
[m, n] = size(sif(ii,:));
[v, I] = max(abs(tfr),[],1);
for j=1:n
rsig(j)=real(tfr(I(j),j));
end

rsig=rsig/(sum(h)/2);
%'rsig' is the reconstructed signal form GLCT, by ridge reconstruction method.
sif3(ii,:)=rsig;
ii
end
figure(15);
v = dbv(ifft(sif3,zpad,2));
S = v(:,1:size(v,2)/2);
m = max(max(v));
imagesc(linspace(0,max_range,zpad),time,(1*S)-m,[-35, 0]);
colorbar;
ylabel('time (s)');
xlabel('range (m)');
title('RTI without clutter rejection with GLCT');
xlim([0 50])
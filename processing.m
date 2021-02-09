%% Data
%%import data from PyCorder
load_data_installer
HDR1 = ScouseTom_getHDR();
HDR = sopen('thumbeegh12.vhdr');
%% Centre Frequency automated detection
%Detect Fc
data = sread(HDR,Inf,0);
Fs=100000;%sampling frequency
L=length(data(1:3*Fs,1)); %length of singal
t= (1/Fs)*(1:L); %duration of signal
Xk=abs(fft(detrend(data(1:3*Fs,1)))); %FFT and amplitude using abs
Xk = Xk(1:L/2); %first half of the tranform due symmetry
f=Fs*(1:L/2)/L; 

%plot waveform
figure('color','white')
plot(data,'color','blue')
%plot the fourier transform
figure('color','white')
plot(f, Xk)
xlim([0 16000])

%detect max value index
[val,idx]= max(Xk(:));
Fc=round(f(idx));
%% Set x-axis in respect to time
%t = dZ_d/ dZ_fn;
T=[0:size(data,1)-1]/Fs;
plot(T,data(:,1))
T_d=(T(:,2500:end-1001));
xlabel('Time (ms)')

%%50hz removal
[c,e]=butter(1,[49,50]/(HDR.SampleRate/2));
%% Extract EMG
%bandpass filter 30-100Hz

[b,a]=butter(3,[150,300]/(HDR.SampleRate/2));
sdata=filtfilt(b,a,data);
sdata_d=(detrend(sdata(2500:end-1001,:)));
plot(T_d,sdata_d(:,:))
title("EMG");
hold on

findchangepts(sdata_d(:,1),'MaxNumChanges',10,'Statistic','rms');
%% RMS
findpeaks(sdata_d(:,2),50000)
envelope(sdata_d(:,:),50000,'rms')
findchangepts(sdata_d(:,1),'MaxNumChanges',10,'Statistic','rms');
RMS_emg= sqrt(1/Fs*abs(sdata_d).^2);

wl1 = 3000;
[up1,lo1] = envelope(sdata_d(:,:),wl1,'rms');
wl2 = 5000;
[up2,lo2] = envelope(sdata_d(:,:),wl2,'rms');
wl3 = 2^15;
[up3,lo3] = envelope(sdata_d(:,:),wl3,'rms');

plot_param = {'Color', [0.6 0.1 0.2],'Linewidth',2}; 
param_small = {'Color',[0.9 0.4 0.1],'Linewidth',2};
param_large = {'Color',[0 0.4 0],'Linewidth',2};


plot(T_d,sdata_d(:,:))
hold on
%p1 = plot(T_d,up1);
%plot(T_d,lo1);
%p2 = plot(T_d,up2);
%plot(T_d,lo2);

p3 = plot(T_d,up3);
hold on
plot(T_d,lo3)
hold off

legend([p1 p2 p3],'wl = 3000','wl = 5000','wl = 50000')
xlim([0 0.04])
title('RMS Envelope')
plot(T_d,RMS_emg)

L=2500;
SV=round(L/2);
Tms=T_d*1000;
W=floor(Tms/(L-SV));

EMGE_MAV(W) = 0;
EMGE_RMS(W) = 0;
EMGE_IAV(W) = 0;
Start=1;
End=L;

for i = 1:W
    EMGE_MAV(i) = mean(abs(sdata_d(Start:End)));
    EMGE_RMS(i) = rms(sdata_d(Start:End));
    EMGE_IAV(i) = sum(abs(sdata_d(Start:End)));
    Start=Start+SV;
    End=End+SV;
end

EMGE_MAV=EMGE_MAV/max(EMGE_MAV);
EMGE_RMS=EMGE_RMS/max(EMGE_RMS);
EMGE_IAV=EMGE_IAV/max(EMGE_IAV);

figure
T=linspace(0,T_d,W);
plot(T,EMGE_MAV,'-g*')
hold on
plot(T,EMGE_RMS,'-mo')
plot(T,EMGE_IAV,'-k')
title('MAV vs RMS vs IAV')
xlabel('Time(s)'),ylabel('Normalized Magnitude'),grid on
xlim([0 time])
legend('MAV','RMS','IAV')
%% Standard deviation of the baseline to find SNR
t = find(T_d>14 &T_d<17);
STDb=std(sdata_d(t,1),0,1);
%% Standard deviation of the activation to find SNR
t = find(T_d>17.7 &T_d<21.8);
STDa=std(sdata_d(t,1),0,1);
%% Modulation of signal
%Bandpass filter to modulate signal
[b,a]=butter(3,[Fc-500,Fc+500]/(HDR.SampleRate/2));
fdata=filtfilt(b,a,data);
plot(T,fdata(:,:))

%% Hilbert transform
%demodulated signal to get combined dZ signal (slow+fast)
dZ=abs(hilbert(fdata));
plot(T,dZ(:,:));

%% Delete artefacts
%get rid of dc offset
plot(detrend(dZ(2500:end-1001,:)));
dZ_d=(detrend(dZ(2500:end-1001,:)));

%% separate fast and slow signals
%Fast signals
%apply high pass filter
[b,a]=butter(3,[0+100,0+500]/(HDR.SampleRate/2));
dZ_fn=filtfilt(b,a,dZ_d);
plot(T_d,dZ_fn(:,:))
title("Fast dZ");
hold on

tb = find(T_d>10 &T_d<11);
STDFb=std(abs(dZ_fn(tb,2)));
meanb=mean(abs(dZ_fn(tb,2)));

ta = find(T_d>25 &T_d<30);
STDFa=std(abs(dZ_fn(ta,4)));
meana=mean(abs(dZ_fn(ta,2)));

%apply low pass filter
%dZ_fn=(detrend(dZ(1000:end-1,:)));
RMS_dZ_fn=dZ_fn.^2;
[b,a]=butter(1, [0.1,30]/(HDR.SampleRate/2));
RMS_dZ_fn1 = filtfilt(b,a,(RMS_dZ_fn.^2));
plot(T_d,RMS_dZ_fn1(:,:))
title("Slow dZ");
hold on

%% PSD
figure;
pwelch(sdata_d(:,1),[],[],[],Fs)
colormap jet
colorbar;
hold on

pwelch(dZ_fn1(:,:),[],[],[],Fs)
hold on

pwelch([sdata_d(:,1),dZ_fn(:,1)],[],[],[],Fs);

%notch filter
d=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1', 45,'HalfPowerFrequency2',65,'DesignMethod','butter','SampleRate',Fs);

dZ_fn= filtfilt(d,dZ_d);

wlen = 2^15;                        % window length (recomended to be power of 2)
nfft = 4*wlen;                      % number of fft points (recomended to be power of 2)
hop = wlen/4;

%emg spectrogram
w1 = hanning(wlen, 'periodic');
[s,f,t]=spectrogram(sdata(:,2),w1,wlen-hop,[300:1:600],Fs);
imagesc(t,f,abs(s))
xlabel('Time, s')
ylabel('Frequency, Hz')

%fast spectrogram
figure
subplot(4,1,1)
[s,f,t]=spectrogram(dZ_fn(:,1),w1,wlen-hop,[-100:0.1:500],Fs);
imagesc(t,f,abs(s))
xlabel('Time, s')
ylabel('Frequency, Hz')
subplot(4,1,2)
[s,f,t]=spectrogram(dZ_fn(:,2),w1,wlen-hop,[-100:0.1:500],Fs);
imagesc(t,f,abs(s))
xlabel('Time, s')
ylabel('Frequency, Hz')
subplot(4,1,3)
[s,f,t]=spectrogram(dZ_fn(:,3),w1,wlen-hop,[-100:0.1:500],Fs);
imagesc(t,f,abs(s))
xlabel('Time, s')
ylabel('Frequency, Hz')
subplot(4,1,4)
[s,f,t]=spectrogram(dZ_fn(:,4),w1,wlen-hop,[-100:0.1:500],Fs);
imagesc(t,f,abs(s))
xlabel('Time, s')
ylabel('Frequency, Hz')
strip(s)
%% subplot
figure
subplot(2,1,1)
[s,f,t]=spectrogram(dZ_fn(:,2),w1,wlen-hop,[50:0.1:600],Fs);
imagesc(t,f,abs(s))
xlabel('Time, s')
ylabel('Frequency, Hz')
subplot(2,1,2)
[s,f,t]=spectrogram(dZ_fn(:,3),w1,wlen-hop,[50:0.1:600],Fs);
imagesc(t,f,abs(s))
xlabel('Time, s')
ylabel('Frequency, Hz')

x = sdata_d(:,:);   % load data
x = x(:, 1);        % get the first channel

% determine the signal parameters
xlen = length(x);                   % signal length
t = (0:xlen-1)/Fs; 
% analysis parameters
wlen = 1024;                        % window length (recomended to be power of 2)
nfft = 4*wlen;                      % number of fft points (recomended to be power of 2)
hop = wlen/4;                       % hop size

TimeRes = wlen/Fs;                  % time resulution of the analysis (i.e., window duration), s
FreqRes = 2*Fs/wlen;                % frequency resolution of the analysis (using Hanning window), Hz


% time-frequency grid parameters
TimeResGrid = hop/Fs;               % time resolution of the grid, s
FreqResGrid = Fs/nfft; 
% perform STFT
w1 = hanning(wlen, 'periodic');
[~, fs, tS, PSD] = spectrogram(x, w1, wlen-hop, nfft, Fs);
Samp = 20*log10(sqrt(PSD.*enbw(w1, Fs))*sqrt(2));

% perform spectral analysis
w2 = hanning(xlen, 'periodic');
[PS, fX] = periodogram(x, w2, nfft, Fs, 'power');
Xamp = 20*log10(sqrt(PS)*sqrt(2));

% plot the signal waveform
figure(1)
subplot(3, 3, [3 2])
plot(t,x)
grid on
xlim([0 max(t)])
ylim([-1.1*max(abs(x)) 1.1*max(abs(x))])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Time, s')
ylabel('Amplitude')
title('The signal in the time domain')
hold on
% plot the spectrum
subplot(3, 3, [4 7])
plot(fX, Xamp)
grid on
xlim([0-10 max(fX)+10])
ylim([min(Xamp)-10 max(Xamp)+10])
view(-90, 90)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Frequency, Hz')
ylabel('Magnitude, dB')
title('Amplitude spectrum of the signal')
hold on
% plot the spectrogram
subplot(3, 3, [5 6 8 9])
surf(tS, fs, Samp)
shading interp
axis tight
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Time, s')
ylabel('Frequency, Hz')
title('Amplitude spectrogram of the signal')
view(0, 90)
hcol = colorbar('East');
set(hcol, 'FontName', 'Times New Roman', 'FontSize', 12)
ylabel(hcol, 'Magnitude, dB')
% display some analysis paramaters
disp(['Frequency resolution of the analysis: ' num2str(FreqRes) ' Hz'])
disp(['Time resolution of the analysis: ' num2str(TimeRes) ' s'])
disp(['Resolution of the frequency grid: ' num2str(FreqResGrid) ' Hz'])
disp(['Resolution of the time grid: ' num2str(TimeResGrid) ' s'])

% % mark the dominant frequencies in the spectrogram
 [~, inds] = max(Samp, [], 1);
 fmax = fs(inds);
 hold on
 plot3(tS, fmax, zeros(length(tS)), 'r')

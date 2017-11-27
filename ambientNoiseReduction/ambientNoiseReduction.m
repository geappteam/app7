clear all
close all
clc

%% Constants and useful variables
voiceFreqBand = [300 4000];
[yNoise,FsNoise] = audioread('bruit_ambiant_16kHz.wav');
[yVoice,FsVoice] = audioread('parole_propre_12kHz.wav');

% Order of the low-pass filter
N = 2;
% Ripple in dB
Rp = 3;

%RSBs ask for the APP7
RSB = [0 5 10];

%% Increasing voice sampling frequency to noise sampling frequency
% 1- Inserting zeros
[K, L] = rat(FsNoise/FsVoice);

yVoiceBuffer = yVoice;
yVoice = zeros(length(yVoice)*K,1);

y = 1;
z = 1;
for i = 1:length(yVoice)
    if(y == K+1 || i == 1)
        y = 1;
        yVoice(i) = yVoiceBuffer(z); 
        z = z + 1;
        if(z == length(yVoiceBuffer)+1)
            z = 1;  
        end
    else
        yVoice(i) = 0; 
    end
    y = y + 1;
end

% 2 - Low-pass (fc = pi/K)
[b,a] = cheby1(N/2,Rp,pi/K);
yVoice = filter(b*K,a,yVoice);

% 3 - Decimating (L)
yVoice = decimate(yVoice,L);

%Verifying that the sample is ok
audiowrite('parole_propre_16kHz.wav',yVoice,FsNoise);

%% Mixing the noise and the voice
gainRBS = zeros(length(RSB),1);
yMixed = zeros(length(yVoice),length(RSB));
for i = 1:length(RSB)
    gainRBS(i) = sqrt(((10^(RSB(i)/10))^(-1))*sum(yVoice.^2)/sum(yNoise.^2));
    yMixed(:,i) = gainRBS(i)*yNoise(1:length(yVoice)) + yVoice;
    
    %Verifying mixing is ok
    audiowrite(strcat('mixed_16kHz_',num2str(RSB(i)),'dB_RSB.wav'),yMixed(:,i),FsNoise);
end


%% Conceiving FIR filters with the inverse FTDS 
%First version (highpass and lowpass)
Ham = hamming(length(yMixed));
Han = hann(length(yMixed));

N = -((length(yMixed)-1)/2):((length(yMixed)-1)/2);
%Lowpass
ThetaCLow = voiceFreqBand(2)*2*pi/FsNoise;
HLow = (ThetaCLow/pi)*sinc(ThetaCLow*N/pi);
HLowHam = HLow.*Ham';
HLowHan = HLow.*Han';

% figure()
% freqz(HLowHan)
% hold on
% freqz(conv(HLowHan,HLowHan))
% hold off

%Highpass
ThetaCHigh = voiceFreqBand(1)*2*pi/FsNoise;
HHigh = (sin(pi.*N)-sin(ThetaCHigh.*N))./(pi*N);
HHighHam = HHigh.*Ham';
HHighHan = HHigh.*Han';

% figure()
% freqz(HHighHan)
% hold on
% freqz(conv(HHighHan,HHighHan))
% hold off

%Resulting frequency response of the bandpass filter
V1_FIR = conv(HLowHan,HHighHan);
figure()
freqz(V1_FIR)

%Second version (bandpass)
%Filter to compare with 
fir1_filter = fir1(1,voiceFreqBand/(FsNoise/2),'bandpass');

%Filter made by equations of transformation
Theta0 = ThetaCHigh+(ThetaCLow-ThetaCHigh)/2;
V2_FIR = 2*HLowHan.*cos(Theta0.*N);

figure()
freqz(fir1_filter)
hold on
freqz(V2_FIR)
hold off

%% Conceiving IIR filters
%Tolerances/Specs
GaindB = 0;
MinGaindB = -0.5;
MaxGaindB = 0.5;

MaxGaindB_below150Hz_beyond5000Hz = 40;

%Constants
N = 1;
Rp = MaxGaindB;
Rs = MaxGaindB_below150Hz_beyond5000Hz;
Fs = [150 5000];

[B_BUTTER,A_BUTTER] = butter(N,voiceFreqBand/(FsNoise/2),'bandpass');
[B_CHEBY1,A_CHEBY1] = cheby1(N,Rp,voiceFreqBand/(FsNoise/2),'bandpass');
[B_CHEBY2,A_CHEBY2] = cheby2(N,Rs,Fs/(FsNoise/2),'bandpass');
[B_ELLIP,A_ELLIP] = ellip(N,Rp,Rs,voiceFreqBand/(FsNoise/2),'bandpass');

figure()
freqz([B_BUTTER,A_BUTTER])
hold on
freqz([B_CHEBY1,A_CHEBY1])
freqz([B_CHEBY2,A_CHEBY2])
freqz([B_ELLIP,A_ELLIP])
legend('BUTTER','CHEBY1','CHEBY2','ELLIP')

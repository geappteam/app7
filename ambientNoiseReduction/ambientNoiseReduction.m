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
%yVoice = decimate(yVoice,L);
[b_dec,a_dec] = cheby1(N/2,0.05,0.8/L);
yVoice = filter(b_dec,a_dec,yVoice);
yVoiceBUFFER = yVoice(1:L:end);
yVoice = yVoiceBUFFER;

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
Order = 200; %Order 40 for first spec (with enough attenuation in stopband) , Order 100 , Order 200 , Order 300
nbPoints = Order+1;

Ham = hamming(nbPoints);
Han = hann(nbPoints);

N = -((nbPoints-1)/2):((nbPoints-1)/2);
%Lowpass
ThetaCLow = voiceFreqBand(2)*2*pi/FsNoise;
HLow = (ThetaCLow/pi)*sinc(ThetaCLow*N/pi);
HLowHam = HLow.*Ham';
HLowHan = HLow.*Han';

figure()
freqz(HLowHan)
hold on
plot(voiceFreqBand(2)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
hold off
title('FIR Lowpass filter frequency response')

%Highpass
ThetaCHigh = voiceFreqBand(1)*2*pi/FsNoise;
HHigh = -sin(ThetaCHigh*N)./(pi*N);
HHigh(Order/2+1) = 1/(2*pi)*(((-ThetaCHigh)-(-pi))+(pi - ThetaCHigh));

HHighHam = HHigh.*Ham';
HHighHan = HHigh.*Han';

figure()
freqz(HHighHan)
hold on
plot(voiceFreqBand(1)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
title('FIR Highpass filter frequency response')
hold off


%Resulting frequency response of the bandpass filter
V1_FIR = conv(HLowHan,HHighHan);
figure(100)
customFreqz(V1_FIR,1)
hold on
plot(voiceFreqBand(1)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(voiceFreqBand(2)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
hold off
title('FIR Bandpass filter (made of highpass and lowpass) frequency response ')

%Second version (bandpass)
%Filter to compare with 
fir1_filter = fir1(Order,voiceFreqBand/(FsNoise/2),'bandpass');

%Filter made by equations of transformation
Theta0 = ThetaCHigh+(ThetaCLow-ThetaCHigh)/2;
BW = 1850;
Theta1 = 2*pi*BW/FsNoise;
HL = (Theta1/pi)*sinc(Theta1*N/pi);
V2_FIR = 2*HL.*cos(Theta0.*N);

figure()
freqz(V2_FIR)
hold on
plot(voiceFreqBand(1)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(voiceFreqBand(2)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
hold off
title('FIR Bandpass filter (by equations of transformation) frequency response')

%Comparison FIR filters V1 and V2
figure(200)
customFreqz(fir1_filter,1)
hold on
customFreqz(V2_FIR,1)
lines = findall(gcf,'type','line');
set(lines(1),'color','red')
set(lines(2),'color','blue')
plot(voiceFreqBand(1)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(voiceFreqBand(2)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
title('FIR Bandpass filter (by equations of transformation) frequency response compared to fir1() Matlab function')
legend('fir()','FIR V2')
hold off

%% Conceiving IIR filters
%Tolerances/Specs
GaindB = 0;
MinGaindB = -0.5;
MaxGaindB = 0.5;

MaxGaindB_below150Hz_beyond5000Hz = 40;

%Constants
N = 5;
Rp = MaxGaindB*0.9;
Rs = MaxGaindB_below150Hz_beyond5000Hz;
Fs = [150 5000];
MaxGaindB2Linear = 10^(MaxGaindB/20);

[B_BUTTER,A_BUTTER] = butter(N,voiceFreqBand/(FsNoise/2),'bandpass');
[B_CHEBY1,A_CHEBY1] = cheby1(N,Rp,voiceFreqBand/(FsNoise/2),'bandpass');
[B_CHEBY2,A_CHEBY2] = cheby2(N,Rs,Fs/(FsNoise/2),'bandpass'); %Ordre 7
[B_ELLIP,A_ELLIP] = ellip(N,Rp,Rs,voiceFreqBand/(FsNoise/2),'bandpass'); %Ordre 5

B_BUTTER = B_BUTTER*MaxGaindB2Linear;
B_CHEBY1 = B_CHEBY1*MaxGaindB2Linear;
B_CHEBY2 = B_CHEBY2*MaxGaindB2Linear;
B_ELLIP = B_ELLIP*MaxGaindB2Linear;

figure(300)
customFreqz(B_BUTTER,A_BUTTER)
hold on
customFreqz(B_CHEBY1,A_CHEBY1)
customFreqz(B_CHEBY2,A_CHEBY2)
customFreqz(B_ELLIP,A_ELLIP)
legend('BUTTER','CHEBY1','CHEBY2','ELLIP')
lines = findall(gcf,'type','line');
set(lines(1),'color','red')
set(lines(2),'color','blue')
set(lines(3),'color','magenta')
set(lines(4),'color','green')
title('IIR Bandpass filters frequency response')

plot(Fs(1)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(Fs(2)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(voiceFreqBand(1)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(voiceFreqBand(2)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(0:1, MinGaindB*ones(length(0:1),1),'k')
plot(0:1, MaxGaindB*ones(length(0:1),1),'k')
plot(0:1, -MaxGaindB_below150Hz_beyond5000Hz*ones(length(0:1),1),'k')

hold off
%% Real specs obtained for each IIR filters
display(tf(B_BUTTER,A_BUTTER))
display(tf(B_CHEBY1,A_CHEBY1))
display(tf(B_CHEBY2,A_CHEBY2))
display(tf(B_ELLIP,A_ELLIP))

%ANSWER: Filtre elliptique odre 5 à utiliser!

%% Ambient noise reduction demonstration

for i = 1:length(RSB)  
    %FIR
    yMixedFiltered = conv(V1_FIR,yMixed(:,i));
    audiowrite(strcat('mixedFiltered_FIR_V1_16kHz_',num2str(RSB(i)),'dB_RSB.wav'),yMixedFiltered,FsNoise);
    yMixedFiltered = conv(V2_FIR,yMixed(:,i));
    audiowrite(strcat('mixedFiltered_FIR_V2_16kHz_',num2str(RSB(i)),'dB_RSB.wav'),yMixedFiltered,FsNoise);
    %IIR
    yMixedFiltered = filter(B_ELLIP,A_ELLIP,yMixed);
    audiowrite(strcat('mixedFiltered_IIR_16kHz_',num2str(RSB(i)),'dB_RSB.wav'),yMixedFiltered,FsNoise);
end
%% Comparison of the FIR filters
figure()
freqz(V1_FIR)
hold on
freqz(V2_FIR)
freqz(B_ELLIP,A_ELLIP)
legend('FIR V1','FIR V2','IIR')
lines = findall(gcf,'type','line');
set(lines(1),'color','red')
set(lines(2),'color','blue')
set(lines(3),'color','green')
plot(Fs(1)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(Fs(2)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(voiceFreqBand(1)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(voiceFreqBand(2)/(FsNoise/2)*ones(513,1), -(512/2):(512/2),'k')
plot(0:1, MinGaindB*ones(length(0:1),1),'k')
plot(0:1, MaxGaindB*ones(length(0:1),1),'k')
plot(0:1, -MaxGaindB_below150Hz_beyond5000Hz*ones(length(0:1),1),'k')
hold off

title('Chosen FIR (version 1 and 2) and IIR Bandpass filters frequency response comparison')
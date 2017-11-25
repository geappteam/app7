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



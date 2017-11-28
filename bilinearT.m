clc
clear all
close all

fc = 2000;
Fs = 16000;
wc = fc*2*pi;

% Filter to transfer to z domain
fonctionEnS = tf(1,  [1/wc^2, sqrt(2)/wc, 1])
bode(fonctionEnS)
title('Analog filter')

% Expected result after transfer as defined by requirements
Wn = fc / (Fs /2);
[expect_b, expect_a] = butter(2, Wn)
figure()
freqz(expect_b,expect_a)
title('Expected filter from butter()');

%% Biliear transform


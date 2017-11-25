clc
clear all
close all

% In this script, we identify the effect of a room on a sound signal
% and we apply it onto someone`s voice

% Impulse response length (To be adjusted)
m = 4001;

% White noice injected into the room
[pureWN, Fs] = audioread('res/bruit_emis_32kHz.wav');

% White noise with the effect of the room applied
[taintedWN, ~] = audioread('res/bruit_capte_32kHz.wav');



% Solving for h
h = identify(pureWN, taintedWN, m);

prediction = filter(h, 1, pureWN);
E = evalErrorPower(taintedWN(m:end), prediction(m:end))

% Applying the filter coefficients to a human voice
[hum, ~] = audioread('res/parole_emise_32khz.wav');

out = filter(h, 1, hum);

bot = audioplayer(out, Fs);
play(bot);
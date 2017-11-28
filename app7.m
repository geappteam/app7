clc
close all
clear all
%% Identification de la salle
disp('Identification de la salle (N obtenu par recherche binaire pour E < 0.5)...');
identification

%% R�duction du bruit ambient
disp('R�duction du bruit ambient');
addpath('ambientNoiseReduction');
ambientNoiseReduction

%% Transformee bilineaire
disp('Transformee bilineaire');
bilinearT
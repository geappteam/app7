clc
close all
clear all
%% Identification de la salle
disp('Identification de la salle (N obtenu par recherche binaire pour E < 0.5)...');
identification

%% Réduction du bruit ambient
disp('Réduction du bruit ambient');
addpath('ambientNoiseReduction');
ambientNoiseReduction

%% Transformee bilineaire
disp('Transformee bilineaire');
bilinearT
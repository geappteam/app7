clc
close all
clear all
%% Identification de la salle
disp('Identification de la salle (N obtenu par recherche binaire pour E < 0.5)...');
identification

%% R�duction du bruit ambient
disp('\n\n  R�duction du bruit ambient...');
ambientNoiseReduction

%% Transformee bilineaire
disp('\n\n  Transformee bilineaire');
bilinearT
close all
clear all
clc
audio_colors

%% Params

% read file
% path = './proyecto/';
% file = '100a120step5';
% path = '';
% file = 'train16';
% beats = beat_track([path file '.wav']);

for i=1:20
    wavfile = sprintf('train%0.0d.wav',i);
    beats   = beat_track(wavfile);
end
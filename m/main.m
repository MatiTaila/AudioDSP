close all
clear all
clc
audio_colors

%% Params

% read file
% path = './proyecto/';
% file = '100a120step5';
path = '';
file = 'train13';

beats = beat_track([path file '.wav']);
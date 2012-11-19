% startup

close all
clear all
home
audio_colors

%% Params

% read file
% [x,fs] = wavread('./proyecto/orishas_wav.wav'); x = x(:,1);
[x,fs] = wavread('train13.wav');
L      = length(x);

% time windows for spectral flux
n_win = 1024;
n_hop = n_win/2;
win   = n_win/fs;
hop   = n_hop/fs;
t     = win/2:hop:L/fs;         % frame times [s]

% induction window
ind_win   = 5;
n_ind_win = ind_win*fs;
t_ind_w   = win/2:hop:ind_win;  % frame times in induction window [s]
N         = length(t_ind_w);    % frames quantity in induction window

n_bins = 2048;

opt.show_plots = 3;

%% Feature detection: Spectral Flux

SFx = SF(x,n_win,n_hop,n_bins,'hamming',opt);

% Filtrado
[B,A]   = butter(10,.4,'low');
SFx_filt = filtfilt(B,A,SFx);

if opt.show_plots >= 1
    figure;
    plot(t,SFx_filt)
    title('\fontsize{16}Spectral Flux filtrada forward y backward')
end

%% Pre-Tracking

t_ind_win = t(1:N);
SFx_filt_ind_win = SFx_filt(1:N);
[agents] =  pre_tracking(SFx_filt_ind_win,n_win,n_hop,fs,opt);

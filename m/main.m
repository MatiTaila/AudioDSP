close all
clear all
clc
audio_colors

%% Params

% read file
% path = './proyecto/';
% file = '90BPM';
% opt.sintetica  = 1;
path = '';
file = 'train13';
opt.sintetica  = 0;

[x,fs] = wavread([path file '.wav']);
if size(x,2)>1, x = x(:,1); end
L      = length(x);

% time windows for spectral flux
n_win = 1024;
% n_hop = n_win/2;
n_hop = 256;
win   = n_win/fs;
hop   = n_hop/fs;
t     = win/2:hop:L/fs;         % frame times [s]

% induction window
ind_win   = 5;
n_ind_win = ind_win*fs;
t_ind_w   = win/2:hop:ind_win;  % frame times in induction window [s]
N         = length(t_ind_w);    % frames quantity in induction window

n_bins = 4096;
% n_bins = 2048;

opt.show_plots = 2;
opt.save_plots = 0;
opt.log        = 1;
opt.wav_write  = 1;
opt.txt_write  = 1;
opt.compu_mati = 1;

if opt.compu_mati
    click_path = '../../../../matlab/audio/beatroot/audio/31-sticks.wav';
else
    click_path = '..\audio\31-sticks.wav';
end

%% Feature detection: Spectral Flux

SFx = SF(x,n_win,n_hop,n_bins,'hamming',opt);
W_SFx=2*pi*fs/n_hop;
W_c=0.28;
% Filtering
% [B,A]   = butter(2,W_c/4,'low');
[B,A]   = butter(2,W_c,'low');
% [B,A]   = butter(10,.4,'low');

SFx_filt = filtfilt(B,A,SFx);
% SFx_filt = moving_avg(SFx,5);

%% Pre-Tracking

t_ind_win = t(1:N);
SFx_filt_ind_win = SFx_filt(1:N);
agents =  pre_tracking(SFx_filt_ind_win,n_win,n_hop,fs,opt);
for i=1:length(agents)
    agents(i).age = 0;
    agents(i).loss = 0;
end

%% Tracking

MaxTabSF   = peak_filt(SFx_filt);
NSFx       = round(n_win/2):n_hop:L;

if opt.show_plots >= 1
    figure(10)
    hold on
    plot(abs(SFx),'color',red2)
    plot(SFx_filt,'color',blue1)
    h = plot(MaxTabSF(:,1),MaxTabSF(:,2),'o','color',red2);
    %set(get(h,'BaseLine'),'LineStyle',':')
    set(h,'MarkerFaceColor','red')
    for i=1:length(MaxTabSF)
        text(MaxTabSF(i,1),MaxTabSF(i,2),num2str(i))
    end
    figure(11)
    h = stem(MaxTabSF(:,1)*n_hop,MaxTabSF(:,2)/max(MaxTabSF(:,2)),'fill','--','color',red2);
    set(get(h,'BaseLine'),'LineStyle',':')
    set(h,'MarkerFaceColor','red')
    hold on
    plot(abs(x/max(x)),'.-')
    plot(NSFx,SFx/max(SFx),'r')
    plot(NSFx,SFx_filt/max(SFx_filt),'g')
    axis tight
    for i=1:length(MaxTabSF)
        text(MaxTabSF(i,1)*n_hop,MaxTabSF(i,2)/max(MaxTabSF(:,2)),num2str(i))
    end
end

n_Pmax     = round(1.20/hop);
n_Pmin     = round(0.24/hop);

MAX_AGENTS         = 30;
MAX_OUTER          = 8;
REDUNDANCY_P_MAX   = 11.6e-3*(fs/n_hop);
REDUNDANCY_PHI_MAX = 23.2e-3*(fs/n_hop);

[agents] = tracking(agents,MaxTabSF,n_Pmax,n_Pmin,hop,MAX_AGENTS,MAX_OUTER,opt,REDUNDANCY_P_MAX,REDUNDANCY_PHI_MAX);
% agents = tracking(agents,MaxTabSF,n_Pmax,n_Pmin,hop,MAX_AGENTS,MAX_OUTER,opt);

%% Referee

S = zeros(length(agents),1);
P = zeros(length(agents),1);
for j=1:length(agents)
    S(j)=agents(j).S(end);
    P(j)=agents(j).Pm(end);
end

[a,b]=max(S);

beats_m = agents(b).Phi';
beats_t = beats_m/(fs/n_hop);

y = zeros(size(x));
while beats_m(end)*n_hop>size(x,1)
    beats_m(end)=[];
end
y(round(beats_m*n_hop),1) = 1;

% y=zeros(size(beats_m));
% y(round(beats_m*n_hop)) = 1;
% if length(y)<length(x)
%     y=[y;zeros(length(x)-length(y),1)];
% else
%     y=y(1:length(x));
% end

% load click sound
[click,fs_click] = wavread(click_path);
if fs_click ~= fs
    fprintf('===========================================\nOJOOOO!!! FRECUENCIAS DE MUESTREO DISTINTAS\n===========================================\n')
end

% wav write
tracked_beats = conv(y,click);
tracked_beats = tracked_beats(1:length(x))/max(abs(tracked_beats(1:length(x))));
signal_out    = (x+tracked_beats)/max(abs(x+tracked_beats)+.0001);

if opt.wav_write
    fprintf('Salvando %s_tracked.wav\n-----------------------------------------------\n',file)
    if opt.compu_mati
        wavwrite(signal_out,fs,['./proyecto/' file '_tracked.wav'])
    else
        wavwrite(signal_out,fs,[file '_tracked.wav'])
    end
end

% text
if opt.txt_write
    fprintf('Salvando %s_tracked.txt\n-----------------------------------------------\n',file)
    % crear un archivo .txt con datos
    fileID = fopen (['./proyecto/' file '_tracked.txt'],'w');
    fprintf (fileID,'%6.2f',beats_t');
    fclose (fileID);
end

%% Test

fprintf('Comparando con el groundTruth de %s_tracked.txt\n-----------------------------------------------\n',file)

% if opt.sintetica
%     fprintf('Valor esperado: %s\n',file);
%     beat_ground_truth(['./proyecto/' file '_tracked.txt']);
% else
%     beat_ground_truth([file '.txt']);
%     beat_ground_truth(['./proyecto/' file '_tracked.txt']);
% end

%% Plots

if opt.show_plots >= 1
    figure;
    plot(x,'color',blue1);
    hold on;
    plot(tracked_beats/2,'color',green1)
    lines = find(y);
    for i=1:length(lines);
        line([lines(i) lines(i)],[-1 1],'linewidth',2.2,'color',red2);
    end
    %     h = stem(ejex(y~=0),y(y~=0),'color',red2,'markersize',0,'linewidth',2);
    axis tight
end

[Psorted,PsortInd]=sort(P);
for i=1:length(agents)
    sorted_agents(i) = agents(PsortInd(i));
    Ssorted(i) = sorted_agents(i).S(end);
end

if opt.show_plots >= 1
    figure;
    subplot(2,1,1)
    h=plot(Psorted,Ssorted,'o','color',red2,'markersize',6);
    set(h,'MarkerFaceColor',red2)
    xlabel('Periodos [muestras]')
    subplot(2,1,2)
    %  plot_with_colormap(1:length(agents),S,'Score','Numero de agente','Score',2,'hot')
    stem(1:length(agents),S,'o','color',blue1,'markersize',6);
end
 
% startup

close all
clear all
home
audio_colors

%% Params

% read file

path = '';
file = 'train13';
% path = './proyecto/';
% file = 'orishas_wav';

[x,fs] = wavread([path file '.wav']);
if size(x,2)>1, x = x(:,1); end
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

opt.show_plots = 2;
opt.save_plots = 0;
opt.log        = 1;
opt.wav_write  = 1;

%% Feature detection: Spectral Flux

SFx = SF(x,n_win,n_hop,n_bins,'hamming',opt);

% Filtrado
[B,A]   = butter(10,.4,'low');
SFx_filt = filtfilt(B,A,SFx);

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

if opt.show_plots >= 1
    figure(10)
    plot(SFx_filt,'color',blue1)
    hold on
    h = plot(MaxTabSF(:,1),MaxTabSF(:,2),'o','color',red2);
    %set(get(h,'BaseLine'),'LineStyle',':')
    set(h,'MarkerFaceColor','red')
    axis tight
    for i=1:length(MaxTabSF)
        text(MaxTabSF(i,1),MaxTabSF(i,2),num2str(i))
    end
end

n_Pmax     = round(1.2/hop);

MAX_AGENTS         = 30;
MAX_OUTER          = 8;
REDUNDANCY_P_MAX   = 11.6e-3*(fs/n_hop);
REDUNDANCY_PHI_MAX = 23.2e-3*(fs/n_hop);

agents = tracking(agents,MaxTabSF,n_Pmax,hop,MAX_AGENTS,MAX_OUTER,opt,REDUNDANCY_P_MAX,REDUNDANCY_PHI_MAX);
% agents = tracking(agents,MaxTabSF,n_Pmax,hop,MAX_AGENTS,MAX_OUTER,opt);

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

y=zeros(size(beats_m));
y(round(beats_m*n_hop)) = 1;
if length(y)<length(x)
    y=[y;zeros(length(x)-length(y),1)];
else
    y=y(1:length(x));
end

% load click sound
[click,fs_click] = wavread('../../../../matlab/audio/beatroot/audio/31-sticks.wav');
if fs_click ~= fs
    fprintf('===========================================\nOJOOOO!!! FRECUENCIAS DE MUESTREO DISTINTAS\n===========================================\n')
end

% wav write
tracked_beats = conv(y,click,'same');
tracked_beats = tracked_beats/max(abs(tracked_beats));
signal_out    = (x+tracked_beats)/max(abs(x+tracked_beats)+.0001);

if opt.wav_write
    fprintf('|---------------------------|\n| Salvando %s_tracked.wav |\n|---------------------------| \n',file)
    wavwrite(signal_out,fs,['./proyecto/' file '_tracked.wav'])
end

%% Plots

if opt.show_plots >= 1
    figure;
    plot(x,'color',blue1);
    hold on;
    lines = find(y);
    for i=1:length(lines);
        line([lines(i) lines(i)],[-1 1],'linewidth',1.8,'color',red2);
    end
%     h = stem(ejex(y~=0),y(y~=0),'color',red2,'markersize',0,'linewidth',2);
    axis tight
end

[Psorted,PsortInd]=sort(P);
for i=1:length(agents)
    sorted_agents(i) = agents(PsortInd(i));
    Ssorted(i) = sorted_agents(i).S(end);
end

figure;
subplot(2,1,1)
h=plot(Psorted,Ssorted,'o','color',red2,'markersize',6);
set(h,'MarkerFaceColor',red2)
xlabel('Periodos [muestras]')
subplot(2,1,2)
%  plot_with_colormap(1:length(agents),S,'Score','Numero de agente','Score',2,'hot')
stem(1:length(agents),S,'o','color',blue1,'markersize',6);
 
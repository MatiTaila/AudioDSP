close all
clear all
clc
audio_colors

%% Stuff

% read file
path = '';
file = 'train13';
[x,fs] = wavread([path file '.wav']); if size(x,2)>1, x = x(:,1); end
L      = length(x);

% time windows for spectral flux
n_win  = 1024;
n_hop  = 100;
win    = n_win/fs;
hop    = n_hop/fs;
t      = win/2:hop:L/fs;
n_bins = 4096;

opt.show_plots = 2;
opt.save_plots = 0;
opt.log        = 1;
opt.wav_write  = 1;
opt.txt_write  = 1;
opt.compu_mati = 1;
opt.sintetica  = 1;

%% Majuga

SFx      = SF(x,n_win,n_hop,n_bins,'hamming',opt);
[B,A]    = butter(2,0.28,'low');
SFx_filt = filtfilt(B,A,SFx);

[MaxTabSF, MinTabSF] = peakdet(SFx_filt, 1e-4);

%%

MaxTabSF2   = peak_filt(SFx_filt);

dists = zeros(size(MaxTabSF,1),1);
for i=1:size(MaxTabSF,1)
    min_l = find(MinTabSF(:,1)<MaxTabSF(i,1),1,'last');
    min_r = find(MinTabSF(:,1)>MaxTabSF(i,1),1,'first');
    if ~isempty(min_l) & ~isempty(min_r)
        dist = min(MaxTabSF(i,2)-MinTabSF(min_l,2),MaxTabSF(i,2)-MinTabSF(min_r,2));
    elseif isempty(min_l) & ~isempty(min_r)
        dist = MaxTabSF(i,2)-MinTabSF(min_r,2);
    elseif ~isempty(min_l) & isempty(min_r)
        dist = MaxTabSF(i,2)-MinTabSF(min_l,2);
    else
        disp('Error al eliminar los maximos espureos de la autocorrelacion');
        break
    end
    dists(i) = dist;
end

th = 0.20;
thr = min(dists)+(max(dists)-min(dists))*th;
MaxTabSF3 = MaxTabSF(dists>thr,:);

median_thr  = zeros(length(SFx_filt),1);
mean_thr    = zeros(length(SFx_filt),1);
median_size = 5;
for i=median_size+1:length(SFx_filt)-median_size
    frame = SFx_filt(i-median_size:i+median_size);
    median_thr(i) = median(frame);
    mean_thr(i) = mean(frame);
end
    
median_thr(1:median_size) = ones(median_size,1)*median_thr(median_size+1);
median_thr(end-median_size+1:end) = ones(median_size,1)*median_thr(end-median_size);

mean_thr(1:median_size) = ones(median_size,1)*mean_thr(median_size+1);
mean_thr(end-median_size+1:end) = ones(median_size,1)*mean_thr(end-median_size);

%%

figure;
hold on
plot(SFx_filt,'color',blue1)
% h1 = plot(MaxTabSF(:,1),MaxTabSF(:,2),'o','color',red3);
% set(h1,'MarkerFaceColor',red3)
h2 = plot(MaxTabSF2(:,1),MaxTabSF2(:,2),'o','color',red3);
set(h2,'MarkerFaceColor',red3)
h3 = plot(MaxTabSF3(:,1),MaxTabSF3(:,2),'o','color',red2);
set(h3,'MarkerFaceColor',red2)
plot(median_thr,'--','color',green1,'linewidth',2)
plot(mean_thr,'--','color',orange1,'linewidth',2)

figure;
hold on
plot(SFx,'color',blue1)
plot(SFx_filt,'color',red2)
plot(median_thr,'color',green1)
plot(mean_thr,'color',orange1)


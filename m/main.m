% startup
close all
clear all
home
audio_colors

%%

[x,fs] = wavread('train13.wav');
L      = length(x);

n_win  = 1024;
n_hop  = n_win/2;
n_bins = 2048;

win    = n_win/fs;
hop    = n_hop/fs;
n      = win/2:hop:L/fs;           % frame times [s]
N      = length(n);                % frames quantity

ind_w  = 5;                        % induction window: 5 seconds
nind_w = length(win/2:hop:ind_w);  % frames quantity in induction window

SFx    = SF(x,n_win,n_hop,n_bins,'hamming');

%% Pre-Tracking

K      = 250;
[acf,lags,bounds] = autocorr(SFx(1:nind_w),K-1);

[MaxTab, MinTab] = peakdet(acf, 1e-4);
% [pks,locs] = findpeaks(acf);



figure;
plot(acf,'color',green2)
hold on
% plot(locs,pks,'.','MarkerSize',20,'color',orange1);
plot(MaxTab(:,1),MaxTab(:,2),'.','MarkerSize',20,'color',blue1);
plot(MinTab(:,1),MinTab(:,2),'.','MarkerSize',15,'color',red2);

dists = zeros(size(MaxTab,1),1);
for i=1:size(MaxTab,1)
    min_l = find(MinTab(:,1)<MaxTab(i,1),1,'last');
    min_r = find(MinTab(:,1)>MaxTab(i,1),1,'first');
    if ~isempty(min_l) & ~isempty(min_r)
        dist = (2*MaxTab(i,2)-MinTab(min_l,2)-MinTab(min_r,2))/2;
    elseif isempty(min_l) & ~isempty(min_r)
        dist = MaxTab(i,2)-MinTab(min_r,2);
    elseif ~isempty(min_l) & isempty(min_r)
        dist = MaxTab(i,2)-MinTab(min_l,2);
    else
        disp('Error al eliminar los maximos espureos de la autocorrelacion');
        break
    end
    dists(i) = dist;
end

th = 0.15;
thr = min(dists)+(max(dists)-min(dists))*th;

MaxTab2 = MaxTab(dists>thr,:);

figure;
plot((0:length(acf)-1)*hop,acf,'color',green2)
hold on
plot((MaxTab2(:,1)-1)*hop,MaxTab2(:,2),'.','MarkerSize',20,'color',blue1);

n_win_rms = 15;
n_overlap_rms = 14;
tt = round(.24/hop)+n_win_rms/2:n_win_rms-n_overlap_rms:length(acf)-n_win_rms/2;
AcfRms = rms(acf(round(.24/hop):end),n_win_rms,n_overlap_rms,0);    

AcfRms2 = interp1(tt,AcfRms,1:length(acf),'linear','extrap');

M = 0.24:(1.2-0.24)/K:1.2-(1.2-0.24)/K;
delta = 0.75;

plot((0:length(acf)-1)*hop,delta*AcfRms2./M,'r')
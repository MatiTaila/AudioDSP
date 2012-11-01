function [P,Phi,S] = pre_tracking(S,n_win,n_hop,fs)

audio_colors
K      = 250;                      % Bins for autocorrelation
ind_w  = 5;                        % induction window: 5 seconds
win    = n_win/fs;                 % induction window length [s]
hop    = n_hop/fs;                 % hop length [s]
n      = win/2:hop:ind_w;          % frame times [s] in induction window
N      = length(n);                % frames quantity in induction window

[acf,lags,bounds] = autocorr(S(1:N),K-1);

[MaxTab, MinTab] = peakdet(acf, 1e-4);

figure;
    plot(acf,'color',green2)
    hold on
    % plot(locs,pks,'.','MarkerSize',20,'color',orange1);
    plot(MaxTab(:,1),MaxTab(:,2),'.','MarkerSize',20,'color',blue1);
    plot(MinTab(:,1),MinTab(:,2),'.','MarkerSize',15,'color',red2);
    title('\fontsize{16}Autocorrelacion')
    legend('Autocorr','Maximos locales','Minimos locales')

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
    title('\fontsize{16}Maximos locales filtrados')

%% Period P    
    
n_win_rms = 16;
n_overlap_rms = 7;

% tt = round(.24/hop)+n_win_rms/2:n_win_rms-n_overlap_rms:length(acf)-n_win_rms/2;
tt = n_win_rms/2:n_win_rms-n_overlap_rms:length(acf)-n_win_rms/2;

AcfRms  = rms(acf,n_win_rms,n_overlap_rms,0);
AcfRms2 = sqrt(mean(acf));

AcfRmsInterp = interp1(tt,AcfRms,1:length(acf),'linear','extrap');

M = 0.24:(1.2-0.24)/K:1.2-(1.2-0.24)/K;
delta = 0.75;

peak_thr = delta*AcfRms2./M;
plot((0:length(acf)-1)*hop,delta*AcfRmsInterp./M,'color',red2,'linewidth',1.9)
plot((0:length(acf)-1)*hop,peak_thr,'color',orange1,'linewidth',2)

peak_thr_interp = interp1(1:K,peak_thr,MaxTab2(:,1),'linear','extrap');
P = MaxTab2(MaxTab2(:,2)>peak_thr_interp,:);

%%

Phi = 0;
S = 0;
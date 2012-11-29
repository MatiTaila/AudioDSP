function [agents] = pre_tracking(S,n_win,n_hop,fs,opt)

audio_colors

win    = n_win/fs;                 % window length [s]
hop    = n_hop/fs;                 % hop length [s]
% ind_w  = n_ind_win/fs;             % induction window: 5 seconds
t = (0:length(S)-1)*hop+win/2;     % frame times [s] in induction window
% t      = win/2:hop:ind_w;          % frame times [s] in induction window
% m      = length(t);                % frames quantity in induction window

L      = length(S);                % induction window size [muestras]

%% Period P    

BPS_min = 0.24;
BPS_max = 1.2;

% Autocorrelation
K         = round(BPS_max*fs/n_hop);     % Bins for autocorrelation
[acf,~,~] = autocorr(S,K-1);

MaxTabAcf = peak_filt(acf);

% Delete first peak (correspond to autocorrlation in 0)
if MaxTabAcf(1,1) == 1
    MaxTabAcf(1,:) = [];
end

% OJOOO TODO no estoy seguro q este bien, pero creo q si
MaxTabAcf(:,1)=MaxTabAcf(:,1)-1;

M       = BPS_min:(BPS_max-BPS_min)/K:BPS_max-(BPS_max-BPS_min)/K;
delta   = 0.75;
AcfRms  = sqrt(mean(acf));
if mean(acf) < 0
    fprintf('=====================================================\nEN EL PRE TRACKING MEAN(ACF) DIO NEGATIVO!!!\n=====================================================\n')
    AcfRms = 0;
end

peak_thr = inf*ones(size(M));
while sum(peak_thr(MaxTabAcf(:,1)) < MaxTabAcf(:,2)')<size(MaxTabAcf,1)/3
    peak_thr = delta*AcfRms./M;
    delta = delta/2;
end

% peak_thr_interp = interp1(1:K,peak_thr,MaxTabAcf(:,1),'linear','extrap');
% P = MaxTabAcf(MaxTabAcf(:,2)>peak_thr_interp,:);

P =  MaxTabAcf(MaxTabAcf(:,2)>peak_thr(MaxTabAcf(:,1))',:); % Periodo en muestras

% init agents
agents(size(P,1)).Pm = 0;
agents(size(P,1)).Phi = 0;
agents(size(P,1)).Sraw = 0;
agents(size(P,1)).Srel = 0;
agents(size(P,1)).S = 0;

% agents(length(P)).Pt = 0;
% agents(length(P)).trains = 0;
% for i = 1:length(agents)
%     agents(i).Pt = P(i,1)/(fs/n_hop);
% end

for i = 1:length(agents)
    agents(i).Pm = P(i,1);
end

if opt.show_plots >= 2
    figure;
    plot(0:K-1,acf,'color',green2)
    hold on
    plot(P(:,1),P(:,2),'.','MarkerSize',20,'color',blue1);
    plot(peak_thr,'color',orange1,'linewidth',2)
    title('\fontsize{16}AUTOCORRELACION DEL SF - Maximos locales filtrados')
    xlabel('\fontsize{16}Muestras')
    legend('Autocorrelacion','Picos','Umbral')
end

%% Phase

MaxTabSF = peak_filt(S);

for i=1:length(agents)
    [score,phi] = S_raw(MaxTabSF,L,agents(i).Pm,n_hop,fs);
    
    agents(i).Phi = phi;
    agents(i).Sraw = score;
    
    if opt.show_plots && i==1
        figure
        h = stem(MaxTabSF(:,1),MaxTabSF(:,2),'fill','--','color',red2);
        set(get(h,'BaseLine'),'LineStyle',':')
        set(h,'MarkerFaceColor','red')
        hold on
        [bp] = beat_train_template(agents(i).Pm,L,phi);
        stem(2000*bp,'color',blue1)
        legend('\fontsize{15}SF','\fontsize{15}Predicted Beats')
    end
end

%% Score

for i = 1:length(agents)
    tmp = 0;
    for j = 1:length(agents)
        if j~=i
            tmp = tmp+r(agents(i).Pm(1),agents(j).Pm(1),n_hop,fs)*agents(j).Sraw;
        end
    end
    agents(i).Srel = 10*agents(i).Sraw+tmp;
end

Sraw = zeros(length(agents),1);
Srel = zeros(length(agents),1);
for i = 1:length(agents)
    Sraw(i) = agents(i).Sraw;
    Srel(i) = agents(i).Srel;
end

maxSraw = max(Sraw);
maxSrel = max(Srel);

S = zeros(length(agents),1);
for i = 1:length(agents)
    agents(i).S = maxSraw*agents(i).Srel/maxSrel;
    S(i)    = agents(i).S;
end

if opt.show_plots > 2
    figure
    plot(Sraw,'-*','color',blue1)
    hold on
    plot(S,'-*','color',red2)
    legend('S raw','S','location','southeast')
end

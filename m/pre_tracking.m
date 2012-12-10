function [agents,BPM_estimado] = pre_tracking(S,n_win,n_hop,fs,opt)

if opt.show_plots
    audio_colors;
end

win = n_win/fs;                  % window length [s]
hop = n_hop/fs;                  % hop length [s]
t   = (0:length(S)-1)*hop+win/2; % frame times [s] in induction window
L   = length(S);                 % induction window size [muestras]

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

MaxTabAcf(:,1)=MaxTabAcf(:,1)-1;

M       = BPS_min:(BPS_max-BPS_min)/K:BPS_max-(BPS_max-BPS_min)/K;
delta   = 0.75;

AcfRms   = sqrt(sum(acf.^2)/length(acf));
peak_thr = delta*AcfRms./M;
step     = .2*peak_thr(end);

peak_thr = inf*ones(size(M));
j=0;
while sum(peak_thr(MaxTabAcf(:,1)) < MaxTabAcf(:,2)')<size(MaxTabAcf,1)/3
    peak_thr = delta*AcfRms./M-j*step;
    if opt.show_plots>=1
        figure(1)
        hold on
        plot(peak_thr,'color',orange1,'linewidth',2)
        if opt.show_plots>=2
            figure(2)
            hold on
            plot(60./((1:length(peak_thr))*n_hop/fs), peak_thr,'color',orange1,'linewidth',2)
            legend('Umbral')
        end
    end
    j=j+1;
end

P =  MaxTabAcf(MaxTabAcf(:,2)>peak_thr(MaxTabAcf(:,1))',:); % Periodo en muestras

% init agents
agents(size(P,1)).pre = 0;
agents(size(P,1)).Pm = 0;
agents(size(P,1)).Phi = 0;
agents(size(P,1)).Sraw = 0;
agents(size(P,1)).Srel = 0;
agents(size(P,1)).S = 0;

if opt.show_plots >= 1
    figure(1)
    plot(0:K-1,acf,'color',green2)
    plot(P(:,1),P(:,2),'.','MarkerSize',20,'color',blue1);
    title('\fontsize{16}AUTOCORRELACION DEL SF - Maximos locales filtrados')
    xlabel('\fontsize{16}Muestras')
    legend('Autocorrelacion','Umbral','Picos')
    hold off
    
    if opt.show_plots >= 2
        figure(2)
        plot(60./((0:K-1)*n_hop/fs),acf,'color',green2)
        plot(60./(P(:,1)*n_hop/fs),P(:,2),'.','MarkerSize',20,'color',blue1);
        title('\fontsize{16}AUTOCORRELACION DEL SF - Maximos locales filtrados')
        xlabel('\fontsize{16}BPMs')
        legend('Autocorrelacion','Picos','location','southeast')
        axis([60/BPS_max-10 60/BPS_min+10 min(acf) max(P(:,2))*1.2])
        hold off
    end
end

%% Phase

MaxTabSF = peak_filt(S);

for i=1:length(agents)
    [score,phi] = S_raw(MaxTabSF,L,P(i,1),n_hop,fs);
    
    agents(i).pre  = phi;
    agents(i).Pm   = P(i,1);
    agents(i).Phi  = 0;
    agents(i).Sraw = score;
    
    if opt.show_plots && i==1
        figure(3)
        h = stem(MaxTabSF(:,1),MaxTabSF(:,2),'fill','--','color',red2);
        set(get(h,'BaseLine'),'LineStyle',':')
        set(h,'MarkerFaceColor','red')
        hold on
        [bp] = beat_train_template(agents(i).Pm,L,phi);
        stem(2000*bp,'color',blue1)
        legend('\fontsize{15}SF','\fontsize{15}Predicted Beats')
        hold off
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
    figure(4)
    plot(Sraw,'-*','color',blue1)
    hold on
    plot(S,'-*','color',red2)
    legend('S raw','S','location','southeast')
    hold off
end

string = num2str(frames2bpm(agents(1).Pm,fs,n_hop));
for i=2:length(agents)
    string = [string '\t' num2str(frames2bpm(agents(i).Pm,fs,n_hop))];
end
if opt.log
    fprintf(['Salida del Pre-Tracking: Periodos:\n' string '\n'])
    fprintf('-----------------------------------------------\n')
end


%% BPM estimado

if size(P,1)> 2
    for i=1:length(P)
        BPM_estimado=60/(P(i,1)/(fs/n_hop));
        if (BPM_estimado > 60/BPS_min) | (BPM_estimado < 60/BPS_max)
            P(i,2)=0;
        end
    end
    [unUsed,I]=max(P);
    if P(I(2),2)~= 0
        BPM_estimado=60/(P(I(2))/(fs/n_hop));
    else
        BPM_estimado=0;
    end
else
    BPM_estimado=0;
end
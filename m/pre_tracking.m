function [agents] = pre_tracking(S,n_win,n_hop,n_ind_win,fs,opt)

audio_colors

win    = n_win/fs;                 % window length [s]
hop    = n_hop/fs;                 % hop length [s]
% ind_w  = n_ind_win/fs;             % induction window: 5 seconds
t = (0:length(S)-1)*hop+win/2;     % frame times [s] in induction window
% t      = win/2:hop:ind_w;          % frame times [s] in induction window
% m      = length(t);                % frames quantity in induction window

L      = length(S);                % induction window size [muestras]

%% Period P    

% Autocorrelation
K         = 250;                      % Bins for autocorrelation
[acf,~,~] = autocorr(S,K-1);

MaxTabAcf = peak_filt(acf);

% Delete first peak (correspond to autocorrlation in 0)
if MaxTabAcf(1,1) == 1
    MaxTabAcf(1,:) = [];
end

BPS_min = 0.24;
BPS_max = 1.2;

M       = BPS_min:(BPS_max-BPS_min)/K:BPS_max-(BPS_max-BPS_min)/K;
delta   = 0.75;
AcfRms  = sqrt(mean(acf));

peak_thr = delta*AcfRms./M;

% peak_thr_interp = interp1(1:K,peak_thr,MaxTabAcf(:,1),'linear','extrap');
% P = MaxTabAcf(MaxTabAcf(:,2)>peak_thr_interp,:);

P =  MaxTabAcf(MaxTabAcf(:,2)>peak_thr(MaxTabAcf(:,1))',:); % Periodo en muestras

% init agents
agents(length(P)).Pm = 0;
agents(length(P)).Pt = 0;
agents(length(P)).trains = 0;

for i = 1:length(agents)
    agents(i).Pt = P(i,1)/fs;
end
for i = 1:length(agents)
    agents(i).Pm = P(i,:);
end

if opt.show_plots >= 2
    figure;
    plot(acf,'color',green2)
    hold on
    plot(P(:,1),P(:,2),'.','MarkerSize',20,'color',blue1);
    plot(peak_thr,'color',orange1,'linewidth',2)
    title('\fontsize{16}AUTOCORRELACION DEL SF - Maximos locales filtrados')
    xlabel('\fontsize{16}Muestras')
    legend('Autocorrelacion','Picos','Umbral')
end

%% Phase

MaxTabS = peak_filt(S);
n_trains = 20;

for i = 1:length(P)
    agents(i).trains = beat_train_template(P(i,1),n_trains,L);
end

if opt.show_plots >= 2
    figure;
    plot(t,S,'color',green2)
    hold on
    plot((MaxTabS(:,1)-1)*hop+win/2,MaxTabS(:,2),'.','MarkerSize',20,'color',blue1);
    title('\fontsize{16}Maximos locales del SF filtrados')
    xlabel('\fontsize{16}Tiempo [s]')
    
    figure;
    plot(S,'color',green2)
    hold on
    plot(MaxTabS(:,1),MaxTabS(:,2),'.','MarkerSize',20,'color',blue1);
    title('\fontsize{16}Maximos locales del SF filtrados')
    xlabel('\fontsize{16}Muestras')
    
    j=1;
    for k = 1:5
        stem(max(S)*agents(j).trains(:,k),'color',colors{k})
    end
    axis tight
end

% Limits to search
lims = zeros(1,size(MaxTabS,1));
for i=1:size(MaxTabS,1)-1
    lims(i)=(MaxTabS(i)+MaxTabS(i+1))/2;
end
lims = [1 lims length(S)];

if opt.show_plots >= 2
    figure;
    stem(MaxTabS(:,1),ones(size(MaxTabS(:,1))),'color',blue1)
    hold on
    stem(lims,ones(length(lims)),'color',green1)
    stem(.5*agents(1).trains(:,10),'color',red2)
    legend('posicion de los maximos del SF','limites de los rangos de busqueda','train template de ejemplo')
end

% keyboard

% select best phase hipothesis
for i = 1:length(agents)
    for j=1:n_trains
        error = 0;
        for k = 1:size(MaxTabS,1)
            % hay q tener en cuenta q puede no haber ningun beat del train
            % template en el rango y q puede haber mas de uno...
            % rango ya esta definido por lims
            indexes = find(agents(i).trains(:,j)==1);
            cands = indexes(find(indexes>=lims(k) & indexes<lims(k+1)));
            if ~isempty(cands)
                difs = abs(MaxTabS(k,1)-cands);
                error = error + min(difs);
            end
            % pararse en cada uno de los MaxTabS(k,1) con el mas cercano de los agents(i).train(:,j)
            % sumar errores en todos los k
        end
        agents(i).error(j) = error;
    end
    [~,ind] = min(agents(i).error);
    agents(i).phi = find(agents(i).trains(:,ind)==1,1,'first');
    % me quedo con el train(:,indice) donde indice es el que tiene menos error acumulado
    % y lo guardo en agents(i).phi
end


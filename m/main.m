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
opt.log        = 2;

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

MAX_AGENTS = 30;
MAX_OUTER  = 8;

MaxTabSF = peak_filt(SFx_filt);
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

for i=1:size(MaxTabSF,1)-1
%     if i == 25, keyboard, end
    
    delete = zeros(length(agents),1);
    inner_count = 0;
    outer_count = 0;
    waiting_agents = [];
    ind = 1;
    for j=1:length(agents)
        
        Tout_R = 0.4*agents(j).Pm(end);
        Tout_L = 0.2*agents(j).Pm(end);
        Tin    = round(46.4e-3/hop);
        m      = MaxTabSF(i,1);
        bp     = agents(j).Phi(end);
        error  = m-bp;
        
        while bp+Tout_R < m
            bp=bp+agents(j).Pm(end);
        end
        
        if(abs(MaxTabSF(i+1,1)-bp)>=abs(error))
            
            if ( error > Tout_R) || ( error < -Tout_L )
                delete(j) = 1;
            else % cae dentro de alguno de los intevalos
                if abs(error)<Tin % inner region
                    inner_count = inner_count+1;
                    agents(j).loss = 0;
                    delta_s = (1-abs(error)/Tout_R)*MaxTabSF(i,2);
%                     agents(j).Pm  = [agents(j).Pm agents(j).Pm(end)+0.25*error];
%                     agents(j).Phi = [agents(j).Phi agents(j).Phi(end)+agents(j).Pm(end)];
                    agents(j).Pm  = agents(j).Pm(end)+0.25*error;
                    agents(j).Phi = agents(j).Phi(end)+agents(j).Pm(end);
                else % outer region
                    % creo hijos
                    outer_count = outer_count+1;
                    agents(j).loss = agents(j).loss + 1;
                    delta_s = -(abs(error)/Tout_R)*MaxTabSF(i,2);
                    P_hijos = {0,error,0.5*error};
                    Phi_hijos = {error,error,0.5*error};
                    for k=1:3
%                         waiting_agents(ind).Pm = [agents(j).Pm agents(j).Pm(end)+P_hijos{k}];
%                         waiting_agents(ind).Phi = [agents(j).Phi agents(j).Phi(end)+Phi_hijos{k}+waiting_agents(ind).Pm(end)];
                        waiting_agents(ind).Pm = agents(j).Pm(end)+P_hijos{k};
                        waiting_agents(ind).Phi = agents(j).Phi(end)+Phi_hijos{k}+waiting_agents(ind).Pm(end);
                        
                        waiting_agents(ind).Sraw = 0;
                        waiting_agents(ind).Srel = 0;
%                         waiting_agents(ind).S = [agents(j).S 0.9*agents(j).S(end)];
                        waiting_agents(ind).S = 0.9*agents(j).S(end);
                        waiting_agents(ind).T_ = 0;
                        waiting_agents(ind).T = 0;
                        waiting_agents(ind).age = 0;
                        waiting_agents(ind).loss = 0;
                        
                        ind = ind+1;
                    end
                end
%                 agents(j).S = [agents(j).S agents(j).S+delta_s];
                agents(j).S = agents(j).S+delta_s;
            end    
        end
    end
    
    if opt.log >= 2
        fprintf('pico numero: %d\n',i)
        fprintf('cantidad de agentes en el inner: %d\n',inner_count)
        fprintf('cantidad de agentes en el outer: %d\n',outer_count)
        fprintf('---------------------------------------\n')
    end
    
    % delete agents / capaz q es mejor no borrarlos sino guardarlos para
    % conservar toda la historia e ir viendo cuando mueren
    j = 1;
    while j<=length(agents)
        if delete(j)
            agents(j) = [];
            delete(j) = [];
        else
            j=j+1;
        end
    end
    
    % redundancy
    
    % loss
    j=1;
    while j<=length(agents)
        if agents(j).loss > MAX_OUTER
            agents(j) = [];
        else
            j=j+1;
        end
    end
    
    % obsolescence
    aux = zeros(length(agents),1);
    for j=1:length(agents)
        aux(j)=agents(j).S(end);
    end
    [maxS,index] = max(aux);
    j=1;
    while j<=length(agents)
        if abs(agents(j).S(end)-maxS) > 0.8*maxS
            agents(j) = [];
        else
            j=j+1;
        end
    end
    
    % replacement
    while length(agents)>MAX_AGENTS
        aux = zeros(length(agents),1);
        for j=1:length(agents)
            aux(j)=agents(j).S(end);
        end
        [a,index] = min(aux);
        agents(index)=[];
    end
    
    
    % add waiting agents
    L = length(agents);
    for j=1:length(waiting_agents)
        agents(L+j)=waiting_agents(j);
    end
    
    % Debug
    P = zeros(length(agents),1);
    S = zeros(length(agents),1);
    for j=1:length(agents)
        P(j)=agents(j).Pm(end);
        S(j)=agents(j).S(end);
    end
    
    
    [a,b]=max(S);
    periodos(i) = agents(b).Pm(end);
    fases(i) = agents(b).Phi(end);
    
    if opt.log >= 1
        if length(agents) == 0
            fprintf('mierda: %d\n=======================================\n',i)
        end
    end
    
end



% y=zeros(size(fases));
% y(round(fases*n_hop))=10;
% [B,A]   = butter(1,0.01,'low');
% y_filt = filtfilt(B,A,y);
% 
% len = min(length(x),length(y));
% wavwrite(x(1:len)+y(1:len)',fs,16,'prueba.wav')


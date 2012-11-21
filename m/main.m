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

%% Pre-Tracking

t_ind_win = t(1:N);
SFx_filt_ind_win = SFx_filt(1:N);
agents =  pre_tracking(SFx_filt_ind_win,n_win,n_hop,fs,opt);

%% Tracking

MaxTabSF = peak_filt(SFx_filt);
if opt.show_plots >= 1
    figure(10)
    plot(SFx_filt,'color',blue1)
    hold on
    h = plot(MaxTabSF(:,1),MaxTabSF(:,2),'o','color',red2);
    %set(get(h,'BaseLine'),'LineStyle',':')
    set(h,'MarkerFaceColor','red')
    axis tight
end

for i=1:size(MaxTabSF,1)
    delete = zeros(length(agents),1);
    inner_count = 0;
    outter_count = 0;
    for j=1:length(agents)
        
        Tout_R = 0.4*agents(j).Pm(end);
        Tout_L = 0.2*agents(j).Pm(end);
        Tin    = round(46.4e-3/hop);
        
        if ( (MaxTabSF(i,1)-agents(j).Phi(end)) > Tout_R) || ( (MaxTabSF(i,1)-agents(j).Phi(end)) < -Tout_L )
            delete(j) = 1;
        else % cae dentro de alguno de los intevalos
            error = MaxTabSF(i,1)-agents(j).T(end);
            if abs(error)<Tin % inner region
                inner_count = inner_count+1;
                agents(j).Pm  = [agents(j).Pm agents(j).Pm+0.25*error];
                agents(j).Phi = [agents(j).Phi agents(j).Phi(end)+0.25*error+agents(j).Pm(end)];
            else % outter region
                % creo hijos
                outter_count = outter_count+1;
            end
        end
    end
    
    fprintf('pico numero: %d\n',i)
    fprintf('cantidad de agentes en el inner: %d\n',inner_count)
    fprintf('cantidad de agentes en el outter: %d\n',outter_count)
    fprintf('---------------------------------------\n')
    
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
    
    if length(agents) == 0
        fprintf('mierda: %d\n=======================================\n',i)
    end
    
end





function beats = beat_track(wavfilename)
% -------------------------------------------------------------------------
% function beats = beat_track(wavfilename);
% -------------------------------------------------------------------------
% Computes the beat tracking of the file given by ./path/file.wav.
% Original implementation: "Beat tracking for multiple applications: A
% multi-agent system architecture with state recovery" Joao Lobato
% Oliveira, Matthew E. P. Davies, Fabien Gouyon and Luis Paulo Reis
% -------------------------------------------------------------------------
% Inputs
%   path    : path to the wav file
%   file    : name of the file
%   opt     : struct of options
% Outputs
%   beats_t : list of beats time in seconds
% -------------------------------------------------------------------------

opt.sintetica  = 0;
opt.show_plots = 0;
opt.save_plots = 0;
opt.log        = 0;
opt.wav_write  = 0;
opt.txt_write  = 0;
opt.compu_mati = 1;
opt.referee    = 1;
opt.cmp_gt     = 0;

if opt.show_plots
    audio_colors;
end

[x,fs] = wavread(wavfilename);
if size(x,2)>1, x = x(:,1); end
L      = length(x);

% time windows for spectral flux
n_win = 1024;
% n_hop = n_win/2;
n_hop = 100;
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

if opt.compu_mati
    click_path = '/home/mat/Documentos/matlab/audio/beatroot/audio/31-sticks.wav';
else
    click_path = '..\audio\31-sticks.wav';
end

%% Feature detection: Spectral Flux

SFx = SF(x,n_win,n_hop,n_bins,'hamming',opt);
W_c=0.28;
% Filtering
[B,A]    = butter(2,W_c,'low');
SFx_filt = filtfilt(B,A,SFx);

%% Pre-Tracking

t_ind_win = t(1:N);
SFx_filt_ind_win = SFx_filt(1:N);
[agents,BPM_estimado] =  pre_tracking(SFx_filt_ind_win,n_win,n_hop,fs,opt);
for i=1:length(agents)
    agents(i).age = 0;
    agents(i).loss = 0;
    agents(i).outer_loss = 0;
    agents(i).wins = 0;
    agents(i).id = i;
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
    hold off
    figure(11)
    h = stem(MaxTabSF(:,1)*n_hop+n_win/2-n_hop,MaxTabSF(:,2)/max(MaxTabSF(:,2)),'fill','--','color',red2);
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
    hold off
end

n_Pmax     = round(1.20/hop);
n_Pmin     = round(0.24/hop);

MAX_AGENTS         = 50;
MAX_OUTER          = 8;
REDUNDANCY_P_MAX   = 11.6e-3*(fs/n_hop);
REDUNDANCY_PHI_MAX = 23.2e-3*(fs/n_hop);

agents = tracking(agents,MaxTabSF,n_Pmax,n_Pmin,hop,MAX_AGENTS,MAX_OUTER,opt,REDUNDANCY_P_MAX,REDUNDANCY_PHI_MAX);
% agents = tracking(agents,MaxTabSF,n_Pmax,n_Pmin,hop,MAX_AGENTS,MAX_OUTER,opt);

%% Referee

S = [agents.S]';
P = [agents.Pm]';
W = [agents.wins]';

[unUsed,Windex] = max(W);
[unUsed,b]=max(S);

if opt.referee == 2
    b = Windex;
end

if opt.log>=1
    fprintf('MaxS: %0.0f\n',b);
    fprintf('MaxW: %0.0f\n',Windex);
end

if opt.referee == 3    
    BPM_max=0;
    BPM_candidatos=0;
    for i=1:length(agents)
        beats_m = agents(i).Phi';
        beats = (beats_m*n_hop+n_win/2-n_hop)/fs;
        bpm_tracker=frames2bpm(fs*mean(diff(beats)),fs,1);
        agentes_bpm(i)=bpm_tracker;
    end;
    
    for i=1:length(agentes_bpm)-1
        for j=i+1:length(agentes_bpm)
            n=multiplo(agentes_bpm(i),agentes_bpm(j));
            if n>1
                BPM_max=[BPM_max max([agentes_bpm(i) agentes_bpm(j)])];
                BPM_candidatos=[BPM_candidatos agentes_bpm(i) agentes_bpm(j)];
            end
        end
    end
    
    % %% Metodo 1
    BPM_ganador_1=max(BPM_max);
    if BPM_ganador_1 ~= 0
        [unUsed,I]=min(abs(agentes_bpm-BPM_ganador_1));
        b=I;
        fprintf('MaxGTR: %0.0f\n',I);
    end;
end

if opt.log>=1
    fprintf('-----------------------------------------------\n');
end

% %% Metodo 2
% %viendo si alguno de los multiplos es parecido al BPM estimado en el
% %pretracking
% [unUsed,I]=min(abs(BPM_candidatos-BPM_estimado));
% BPM_ganador=BPM_candidatos(I);
% if abs(BPM_ganador-BPM_estimado)<10
%     [unUsed,I]=min(abs(agentes_bpm-BPM_ganador));%para sacar cual agente es
%     if abs(S(b)-S(I))<6000  b=I;% comparo el puntaje con el de Matias
%     else %segundo metodo, este elige el maximo BPM que es multiplo, sino usa el b de Matias
%         BPM_ganador_1=max(BPM_max);
%         if BPM_ganador_1 ~= 0
%             [unUsed,I]=min(abs(agentes_bpm-BPM_ganador_1));
%             if abs(S(b)-S(I))<6000  b=I;% comparo el puntaje con el de Matias
%             end;
%         end;
%     end;
% end;

if opt.log>=1
    if Windex ~= b
        fprintf('=====================================================\nOJOOOO!!! DISTINTOS REFEREES DAN DISTINTOS RESULTADOS\n=====================================================\n')
    end
end

beats_m = agents(b).Phi';
beats   = (beats_m*n_hop+n_win/2-n_hop)/fs;

y = zeros(size(x));
while beats_m(end)*n_hop>size(x,1)
    beats_m(end)=[];
end
y(round(beats_m*n_hop+n_win/2-n_hop),1) = 1;

if opt.wav_write>=1
    % load click sound
    [click,fs_click] = wavread(click_path);
    if fs_click ~= fs
        fprintf('===========================================\nOJOOOO!!! FRECUENCIAS DE MUESTREO DISTINTAS\n===========================================\n')
    end
    % wav write
    tracked_beats = conv(y,click);
    tracked_beats = tracked_beats(1:length(x))/max(abs(tracked_beats(1:length(x))));
    signal_out    = (x+tracked_beats)/max(abs(x+tracked_beats)+.0001);
    if opt.log>=1
        fprintf('Salvando %s_tracked.wav\n-----------------------------------------------\n',wavfilename(1:end-4))
    end
    if opt.compu_mati
        wavwrite(signal_out,fs,['./proyecto/results/' wavfilename(1:end-4) '_tracked.wav'])
    else
        wavwrite(signal_out,fs,['..\train\Tracked_Songs\' wavfilename(19:end-4) '_tracked.wav'])
    end
end

% text
if opt.txt_write
    fprintf('Salvando %s_tracked.txt\n-----------------------------------------------\n',wavfilename(1:end-4))
    % crear un archivo .txt con datos
    fileID = fopen (['./proyecto/results/' wavfilename(1:end-4) '_tracked.txt'],'w');
    fprintf (fileID,'%6.2f',beats');
    fclose (fileID);
end

%% Plots

if opt.show_plots >= 1
    figure(50);
    plot(x,'color',blue1);
    hold on;
    lines = find(y);
    for i=1:length(lines);
        line([lines(i) lines(i)],[-1 1],'linewidth',2.2,'color',red2);
    end
    axis tight
    hold off
end

[Psorted,PsortInd]=sort(P);
for i=1:length(agents)
    sorted_agents(i) = agents(PsortInd(i));
    Ssorted(i) = sorted_agents(i).S(end);
end

aux = zeros(length(agents),1);
for i=1:length(agents)
    aux(i) = frames2bpm(Psorted(i),fs,n_hop);
end
if opt.show_plots >= 1
    figure(100);
    subplot(2,1,1)
    h=plot(aux,Ssorted,'o','color',red2,'markersize',6);
    set(h,'MarkerFaceColor',red2)
    for k=1:length(agents)
        text(aux(k),Ssorted(k),num2str(k))
    end
    xlabel('Periodos [BPM]')
    subplot(2,1,2)
    %  plot_with_colormap(1:length(agents),S,'Score','Numero de agente','Score',2,'hot')
    stem(1:length(agents),S,'o','color',blue1,'markersize',6);
    hold off
end

%% Performance

if opt.cmp_gt
    gt = load([wavfilename(1:end-4) '.lab']);
    tracked = load([wavfilename(1:end-4) '_tracked.txt']);
    [cmlC1,cmlT1,amlC1,amlT1] = be_continuityBased(gt,beats);
    [f1,p1,r1,a1] = be_fMeasure(gt,tracked);
    if opt.log>=1
        fprintf('Performance:\n\tCont-Based:\t cC:\t%0.2f\tcT:\t%0.2f\taC:\t%0.2f\taT:\t%0.2f\n\tF-Mesure:\t f:\t%0.2f\tp:\t%0.2f\tr:\t%0.2f\ta:\t%0.2f\n-------------------------------------------------------------------------------------\n',cmlC1,cmlT1,amlC1,amlT1,f1,p1,r1,a1)
    end
end
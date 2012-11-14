function [score] = S_raw(MaxTabSF,N,Pi,phi,n_win,n_hop,fs,L,n_ind_win,train)

%Pi tiene que estar en muestras
%phi tiene que estar en muestras




global bp
global bp1
global I
global n_error
global aux
global n_Cluster_Width
global n_Tin
global n_Tout
global Pii
global Pm
global Flux
aux=MaxTabSF(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure  
h = stem(MaxTabSF(:,1),MaxTabSF(:,2),'fill','--');
set(get(h,'BaseLine'),'LineStyle',':')
set(h,'MarkerFaceColor','red')
hold on
[bp] = beat_train_template(Pi,1,N);
%se cuenta a partir de la primer muestra y termina en la muestra
%phi hay que verificar esto
bp=[zeros(phi,1); bp(1:length(bp)-phi)];
stem(2000*bp,'b')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



hop    = n_hop/fs; 
n_Tin=fix((46.4*10^(-3))/hop);% en muestras
%% si le pasamos Pi y phi_i hay que armar el bp y calcular Pm

[bpi]=find(bp); %lugares indice de la prediccion de beat
n_Tout=fix(0.4*Pi); % en muestras
Pm=length(bpi);%fix((n_ind_win/n_hop)/n_Pi);
score=0;
for i=1:length(bpi),
    [n_error,I_MaxTabSF]=min(abs(MaxTabSF(:,1)-bpi(i)));
    if n_error <= n_Tout
        if n_error <= n_Tin
            delta_s=(1-n_error/n_Tout)*MaxTabSF(I_MaxTabSF,2);
        else  delta_s=-n_error/n_Tout*MaxTabSF(I_MaxTabSF,2);
        end;
    else delta_s=0;
     end;

score= score + delta_s*Pi/Pm ;  
end;


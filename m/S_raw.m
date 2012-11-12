function [score] = S_raw(MaxTabSF,n_win,n_hop,fs,L,n_ind_win,train)


N =length(round(n_win/2):n_hop:L); % frames totales


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure  
h = stem(MaxTabSF(:,1),MaxTabSF(:,2),'fill','--');
set(get(h,'BaseLine'),'LineStyle',':')
set(h,'MarkerFaceColor','red')
hold on
stem(2000*train,'b')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global bp
global I
global n_error
global aux
global n_Cluster_Width
global n_Tin
global n_Tout
global n_Pi
global Pm
aux=MaxTabSF(:,1);



n_Tin=fix((46.4*10^(-3))/hop);% en muestras
%% si le pasamos Pi y phi_i hay que armar el bp y calcular Pm
[bp]=find(train); %lugares indice de la prediccion de beat
n_Pi= bp(2)-bp(1);% en muestras
n_Tout=fix(0.4*n_Pi); % en muestras
Pm=length(bp);%fix((n_ind_win/n_hop)/n_Pi);
score=0;
for i=1:length(bp),
    [n_error,I_MaxTabSF]=min(abs(MaxTabSF(:,1)-bp(i)));
    if n_error <= n_Tout
        if n_error <= n_Tin
            delta_s=(1-n_error/n_Tout)*MaxTabSF(I_MaxTabSF,2);
        else  delta_s=-n_error/n_Tout*MaxTabSF(I_MaxTabSF,2);
        end;
    else delta_s=0;
     end;

score= score + delta_s*n_Pi/Pm ;  
end;


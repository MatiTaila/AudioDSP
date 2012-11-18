function [score,phase] = S_raw(MaxTabSF,N,Pi,n_hop,fs)

%Pi tiene que estar en muestras







hop    = n_hop/fs; 
n_Pmax=round(1.2/hop);
M=103.3594;
% M=103.3594/2;
n_step_phi=ceil(n_Pmax/M);

n_Tin=fix((46.4*10^(-3))/hop);% en muestras
n_Tout_R=fix(0.4*Pi); % en muestras
n_Tout_L=fix(0.2*Pi);
cantidad_phase=fix(Pi/n_step_phi);
score=zeros(cantidad_phase,1);
aux=0;

for j=1:cantidad_phase
    [bp] = beat_train_template(Pi,1,N,(j-1)*n_step_phi); % genero el tren de beat variable en phi
    [bpi]=find(bp); %lugares indice de la prediccion de beat
    for i=1:length(bpi),
        [NoSeUsa,I_MaxTabSF]=min(abs(MaxTabSF(:,1)-bpi(i)));
        n_error_rel=MaxTabSF(I_MaxTabSF,1)-bpi(i);
        inner_window= ((bpi(i)-n_Tout_L <= n_error_rel)&(n_error_rel <= bpi(i)-n_Tin)) | ((bpi(i)+n_Tin<=n_error_rel)&(n_error_rel<=bpi(i)+n_Tout_R)) | abs(n_error_rel)<n_Tin;
        
        if inner_window
            if abs(n_error_rel) <= n_Tin
                delta_s=(1-abs(n_error_rel)/n_Tout_R)*MaxTabSF(I_MaxTabSF,2);
            else  delta_s=(n_error_rel/n_Tout_R)*MaxTabSF(I_MaxTabSF,2);
            end;
        else delta_s=0;
        end;
        aux= aux + delta_s*Pi/n_Pmax ;
    end;
     score(j,1)=aux;
     aux=0;
   
end;

[score,I]=max(score);
phase=(I-1)*n_step_phi; %es I-1 porque arranca en fase 0

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure  
h = stem(MaxTabSF(:,1),MaxTabSF(:,2),'fill','--');
set(get(h,'BaseLine'),'LineStyle',':')
set(h,'MarkerFaceColor','red')
hold on
[bp] = beat_train_template(Pi,1,N,phase);
stem(2000*bp,'b')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

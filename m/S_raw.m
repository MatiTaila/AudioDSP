function [score,phi] = S_raw(MaxTabSF,L,Pi,n_hop,fs)
% -------------------------------------------------------------------------
% function [score,phi] = S_raw(MaxTabSF,L,Pi,n_hop,fs)
% -------------------------------------------------------------------------
% MaxTabSF : Position of spectral flux maxima
% L        : Number of frames in induction window
% Pi       : Periodo en muestras
% n_hop    : Hop in muestras
% fs       : Sample rate
% -------------------------------------------------------------------------

BPS_max    = 1.2;
hop        = n_hop/fs;
n_Pmax     = round(BPS_max/hop);
M          = round(BPS_max*fs/n_hop);
% M          = 103.3594;                  % FIXME
n_step_phi = ceil(n_Pmax/M);

n_Tin      = round(46.4e-3/hop);        % en muestras
n_Tout_R   = round(0.4*Pi);             % en muestras
n_Tout_L   = round(0.2*Pi);
N_phi    = fix(Pi/n_step_phi);
scores      = zeros(N_phi,1);

for j=1:N_phi
    cum_err = 0;
    bp  = beat_train_template(Pi,L,(j)*n_step_phi); % genero el tren de beat variable en phi
    bpi = find(bp);                                     % lugares indice de la prediccion de beat
    for i=1:length(bpi),
        [NoSeUsa,I_MaxTabSF] = min(abs(MaxTabSF(:,1)-bpi(i)));
        n_error_rel          = MaxTabSF(I_MaxTabSF,1)-bpi(i);
        inner_window         = ( (bpi(i)-n_Tout_L <= n_error_rel) & (n_error_rel <= bpi(i)-n_Tin) ) ...
            | ( (bpi(i)+n_Tin <= n_error_rel) & (n_error_rel <= bpi(i)+n_Tout_R) ) ...
            | abs(n_error_rel) < n_Tin;
        
        if inner_window
            if abs(n_error_rel) <= n_Tin
                delta_s = (1-abs(n_error_rel)/n_Tout_R)*MaxTabSF(I_MaxTabSF,2);
            else
                delta_s = -(abs(n_error_rel)/n_Tout_R)*MaxTabSF(I_MaxTabSF,2);
            end
        else
            delta_s = -2*MaxTabSF(I_MaxTabSF,2);
        end
        cum_err = cum_err + delta_s*Pi/n_Pmax;
    end
    scores(j,1) = cum_err;
end

[score,I] = max(scores);
phi       = (I)*n_step_phi; % es I-1 porque arranca en fase 0

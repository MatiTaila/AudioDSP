function [train] = beat_train_template(P,N,n_win,phi)
% P: period [muestras]
% N: number of templates
% n_win: windows size
% phi: phase of train

M = zeros(n_win,N);

lags = 0:P/N:P-P/N;

for i = 1:N
    M(round(1+lags(i)):P:n_win,i)=1;
end

train=[zeros(phi,1);M(1:length(M)-phi)];
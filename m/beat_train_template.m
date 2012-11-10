function [M] = beat_train_template(P,N,n_win)
% P: period [muestras]
% N: number of templates
% n_win: windows size

M = zeros(n_win,N);

lags = 0:P/N:P-P/N;

for i = 1:N
    M(round(1+lags(i)):P:n_win,i)=1;
end
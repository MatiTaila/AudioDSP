function [train] = beat_train_template(P,n_win,phi)
% P: period [muestras]
% N: number of templates
% n_win: windows size
% phi: phase of train

if ~exist('phi','var')
    phi = 1;
end
% keyboard
% M = zeros(n_win,N);
% for i = 1:N
%     M(phi:P:n_win,i)=1;
% end

% train=[zeros(phi,1);M(1:length(M)-phi)];
train = zeros(n_win,1);
train(phi:P:n_win,1)=1;
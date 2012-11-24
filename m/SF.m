function [SF] = SF(x,n_win,n_hop,n_bins,win,opt)

% path = 'train13';
% [x,fs] = wavread(['train13' '.wav']);
L      = length(x);
% gt     = load([path '.txt']);
% t      = (0:L-1)/fs;
% n_win  = 1024;
% n_hop  = n_win/2;
% win    = n_win/fs;
% hop    = n_hop/fs;
% n_bins = 2*1024;

switch win
    case 'hanning'
        window = hann(n_win);
    case 'hamming'
        window = hamming(n_win);
end

N    = length(round(n_win/2):n_hop:L);

% window = hamming(n_win);
fx = zeros(n_bins,N);
for i = 1:N
    t_st  = (i-1)*n_hop;
    if length(x)>=t_st+n_win
        frame = x(t_st+1:t_st+n_win).*window;
    else
        frame = x(t_st+1:end).*hamming(length(t_st+1:L));
    end
    fx(:,i) = fft(frame,n_bins);
end

diff = abs([fx(:,2:end) zeros(n_bins,1)]) - abs(fx);

h = @(x) (x+abs(x))/2;

SF = sum(h(diff));

if opt.show_plots >= 3
    figure; plot(SF)
end

% [B,A] = butter(7,2500/(fs/2),'low');
% SF_filt = filter(B,A,SF);

% figure; plot(SF_filt)

%%

% window    = 2^10;
% noverlap  = window/2;
% nfft      = 2^12;
% [S F T P] = spectrogram(x,window,noverlap,nfft,fs);
% figure
% plot_spect(T,F,P)

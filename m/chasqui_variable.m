function chasqui_variable(fs,t,BPMini,BPMfin,step,path,opt)
% -------------------------------------------------------------------------
% chasqui_variable(fs,t,BPMini,BPMfin,step,path,opt)
% -------------------------------------------------------------------------
% Example: 
%   chasqui_variable(44100,10,90,180,5,'./proyecto/90a180step5.wav',1);
% -------------------------------------------------------------------------
% Inputs:
%   fs     : sample rate
%   t      : total time of de generated synthetic signal
%   BPMini : initial beat frequency in BPMs
%   BPMfin : final beat frequency in BPMs
%   step   : step for changing veat frequency
%   path   : path for saving .wav
%   opt    : pathes for different computers: 1=Mat, 0=GTR
% -------------------------------------------------------------------------

% load click sound
if opt
    [click,fs_click] = wavread('../../../../matlab/audio/beatroot/audio/53-cymbell.wav');
else
    [click,fs_click] = wavread('..\audio\31-sticks.wav');
end
if fs_click ~= fs
    fprintf('===========================================\nOJOOOO!!! FRECUENCIAS DE MUESTREO DISTINTAS\n===========================================\n')
end

BPMs = min(BPMini,BPMfin):step:max(BPMini,BPMfin);

P = 60*fs./BPMs;
L = length(P);
n_win = fs*t/L;

train = zeros(L*n_win,1); train(1)=1;

for i=1:L
    train(fix(find(train==1,1,'last'):P(i):n_win*i))=1;
end

tracked_beats = conv(train,click);
tracked_beats = tracked_beats/max(abs(tracked_beats)+.00001);
wavwrite(tracked_beats,fs,path);


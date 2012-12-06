close all
clear all
clc
audio_colors

N = 20;

%% Params

% read file
% path = './proyecto/';
% file = '100a120step5';

%% Set 1

% for i=1:20
%     wavfile = sprintf('train%0.0d.wav',i);
%     beats   = beat_track(wavfile);
% end

%% Set 2

cmlC2 = zeros(N,1); cmlT2 = zeros(N,1); amlC2 = zeros(N,1); amlT2 = zeros(N,1);
f2    = zeros(N,1); p2    = zeros(N,1); r2    = zeros(N,1); a2    = zeros(N,1);

for i=1:N
    wavfile = sprintf('datos2_%0.0d_A.wav',i);
    beats   = beat_track(wavfile);
    gt      = load([wavfile(1:end-4) '.lab']);
    tracked = load([wavfile(1:end-4) '_tracked.txt']);
    [cmlC2(i),cmlT2(i),amlC2(i),amlT2(i)] = be_continuityBased(gt,tracked);
    [f2(i),p2(i),r2(i),a2(i)] = be_fMeasure(gt,tracked);
end

beContinuityBased = [(1:20)' cmlC2 cmlT2 amlC2 amlT2];
disp(beContinuityBased2);
beFMesure2 = [(1:20)' f2 p2 r2 a2];
disp(beFMesure2);

%% Set 3

cmlC3 = zeros(N,1); cmlT3 = zeros(N,1); amlC3 = zeros(N,1); amlT3 = zeros(N,1);
f3    = zeros(N,1); p3    = zeros(N,1); r3    = zeros(N,1); a3    = zeros(N,1);

for i=1:N
    wavfile = sprintf('datos3_%0.0d_A.wav',i);
    beats   = beat_track(wavfile);
    gt      = load([wavfile(1:end-4) '.lab']);
    tracked = load([wavfile(1:end-4) '_tracked.txt']);
    [cmlC3(i),cmlT3(i),amlC3(i),amlT3(i)] = be_continuityBased(gt,tracked);
    [f3(i),p3(i),r3(i),a3(i)] = be_fMeasure(gt,tracked);
end

beContinuityBased3 = [(1:20)' cmlC3 cmlT3 amlC3 amlT3];
disp(beContinuityBased3);
beFMesure3 = [(1:20)' f3 p3 r3 a3];
disp(beFMesure3);

close all
clear all
clc
audio_colors

N  = 20;
fs = 44100;

%% Generate train ground truth

% for i=1:N
%     wavfile = sprintf('train%0.0d.txt',i);
%     read_ground_truth(wavfile);
% end

%% Synth

% read file
% path = './proyecto/';
% file = '100a120step5';

%% Set 1

cmlC1 = zeros(N,1); cmlT1 = zeros(N,1); amlC1 = zeros(N,1); amlT1 = zeros(N,1);
f1    = zeros(N,1); p1    = zeros(N,1); r1    = zeros(N,1); a1    = zeros(N,1);
BPM1  = zeros(N,3);

for i=1:N
    wavfile = sprintf('train%0.0d.wav',i);
    beats   = beat_track(wavfile);
    gt      = load([wavfile(1:end-4) '.lab']);
    tracked = load([wavfile(1:end-4) '_tracked.txt']);
    [cmlC1(i),cmlT1(i),amlC1(i),amlT1(i)] = be_continuityBased(gt,tracked);
    [f1(i),p1(i),r1(i),a1(i)] = be_fMeasure(gt,tracked);
    [~,bpm]=beat_ground_truth([wavfile(1:end-4) '.lab'],0);
    BPM1(i,:) = [i bpm frames2bpm(fs*mean(diff(beats)),fs,1)];
end

beContinuityBased1 = [(1:20)' cmlC1 cmlT1 amlC1 amlT1];
disp(beContinuityBased1);
beFMesure1 = [(1:20)' f1 p1 r1 a1];
disp(beFMesure1);

%% Set 2

cmlC2 = zeros(N,1); cmlT2 = zeros(N,1); amlC2 = zeros(N,1); amlT2 = zeros(N,1);
f2    = zeros(N,1); p2    = zeros(N,1); r2    = zeros(N,1); a2    = zeros(N,1);
BPM2  = zeros(N,3);

for i=1:N
    wavfile = sprintf('datos2_%0.0d_A.wav',i);
    beats   = beat_track(wavfile);
    gt      = load([wavfile(1:end-4) '.lab']);
    tracked = load([wavfile(1:end-4) '_tracked.txt']);
    [cmlC2(i),cmlT2(i),amlC2(i),amlT2(i)] = be_continuityBased(gt,tracked);
    [f2(i),p2(i),r2(i),a2(i)] = be_fMeasure(gt,tracked);
    BPM2(i,:) = [i frames2bpm(mean(diff(gt))*fs,fs,1) frames2bpm(fs*mean(diff(beats)),fs,1)];
end

beContinuityBased2 = [(1:20)' cmlC2 cmlT2 amlC2 amlT2];
disp(beContinuityBased2);
beFMesure2 = [(1:20)' f2 p2 r2 a2];
disp(beFMesure2);

%% Set 3

cmlC3 = zeros(N,1); cmlT3 = zeros(N,1); amlC3 = zeros(N,1); amlT3 = zeros(N,1);
f3    = zeros(N,1); p3    = zeros(N,1); r3    = zeros(N,1); a3    = zeros(N,1);
BPM3  = zeros(N,3);

for i=1:N
    wavfile = sprintf('datos3_%0.0d_A.wav',i);
    beats   = beat_track(wavfile);
    gt      = load([wavfile(1:end-4) '.lab']);
    tracked = load([wavfile(1:end-4) '_tracked.txt']);
    [cmlC3(i),cmlT3(i),amlC3(i),amlT3(i)] = be_continuityBased(gt,tracked);
    [f3(i),p3(i),r3(i),a3(i)] = be_fMeasure(gt,tracked);
    BPM3(i,:) = [i frames2bpm(mean(diff(gt))*fs,fs,1) frames2bpm(fs*mean(diff(beats)),fs,1)];
end

beContinuityBased3 = [(1:20)' cmlC3 cmlT3 amlC3 amlT3];
disp(beContinuityBased3);
beFMesure3 = [(1:20)' f3 p3 r3 a3];
disp(beFMesure3);

%% Save results

save('resultados_2012.12.09_W28_sin_limites','beContinuityBased1','beFMesure1','BPM1','beContinuityBased2','beFMesure2','BPM2','beContinuityBased3','beFMesure3','BPM3')

%% Load Results

% r1 = load('resultados');
% r2 = load('resultados_referee2');
% r3 = load('resultados_sin_Pi_Pmax');
% r4 = load('resultados_desfavoreciendo_multiplos');
% r5 = load('resultados_desfavoreciendo_multiplos_referee_2');
% r6 = load('resultados_desfavoreciendo_multiplos_pero_no_tanto');
% 
% beContinuityBased = [...
%     mean([r1.beContinuityBased1;r1.beContinuityBased2;r1.beContinuityBased3]);
%     mean([r2.beContinuityBased1;r2.beContinuityBased2;r2.beContinuityBased3]);
%     mean([r3.beContinuityBased1;r3.beContinuityBased2;r3.beContinuityBased3]);
%     mean([r4.beContinuityBased1;r4.beContinuityBased2;r4.beContinuityBased3]);
%     mean([r5.beContinuityBased1;r5.beContinuityBased2;r5.beContinuityBased3]);
%     mean([r6.beContinuityBased1;r6.beContinuityBased2;r6.beContinuityBased3]) ]
% 
% beFMesure = [...
%     mean([r1.beFMesure1;r1.beFMesure2;r1.beFMesure3]);
%     mean([r2.beFMesure1;r2.beFMesure2;r2.beFMesure3]);
%     mean([r3.beFMesure1;r3.beFMesure2;r3.beFMesure3]);
%     mean([r4.beFMesure1;r4.beFMesure2;r4.beFMesure3]);
%     mean([r5.beFMesure1;r5.beFMesure2;r5.beFMesure3]);
%     mean([r6.beFMesure1;r6.beFMesure2;r6.beFMesure3]) ]
% 
% cocientes = [...
%     [r1.BPM1(:,2)./r1.BPM1(:,3);r1.BPM2(:,2)./r1.BPM2(:,3);r1.BPM3(:,2)./r1.BPM3(:,3)]';
%     [r2.BPM1(:,2)./r2.BPM1(:,3);r2.BPM2(:,2)./r2.BPM2(:,3);r2.BPM3(:,2)./r2.BPM3(:,3)]';
%     [r3.BPM1(:,2)./r3.BPM1(:,3);r3.BPM2(:,2)./r3.BPM2(:,3);r3.BPM3(:,2)./r3.BPM3(:,3)]';
%     [r4.BPM1(:,2)./r4.BPM1(:,3);r4.BPM2(:,2)./r4.BPM2(:,3);r4.BPM3(:,2)./r4.BPM3(:,3)]';
%     [r5.BPM1(:,2)./r5.BPM1(:,3);r5.BPM2(:,2)./r5.BPM2(:,3);r5.BPM3(:,2)./r5.BPM3(:,3)]';
%     [r6.BPM1(:,2)./r6.BPM1(:,3);r6.BPM2(:,2)./r6.BPM2(:,3);r6.BPM3(:,2)./r6.BPM3(:,3)]' ]

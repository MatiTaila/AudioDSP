close all
clear all
clc

%% 90 BPM
beats = beat_track('90BPM.wav');
figure(50)
xlabel('\fontsize{15}Tiempo [muestras]')
ylabel('\fontsize{15}Amplitud')
legend('\fontsize{13}Senal','\fontsize{13}Beats')
print('-depsc','proyecto/3_entrega_FINAL/pics/bpm90beats');

%% BPM variable de 90 a 120
startup
beats = beat_track('90a100step5.wav');
figure(50)
xlabel('\fontsize{15}Tiempo [muestras]')
ylabel('\fontsize{15}Amplitud')
legend('\fontsize{13}Senal','\fontsize{13}Beats')
print('-depsc','proyecto/3_entrega_FINAL/pics/bpm90a100beats');

%%
startup
beats = beat_track('train13.wav');
figure(50)
xlabel('\fontsize{15}Tiempo [muestras]')
ylabel('\fontsize{15}Amplitud')
legend('\fontsize{13}Senal','\fontsize{13}Beats')

figure(1)
xlabel('\fontsize{15}Tiempo [Muestras]')
ylabel('\fontsize{15}Amplitud')
legend('\fontsize{15}Umbral','\fontsize{15}Autocorrelacion del SF','\fontsize{15}Picos filtrados')
axis tight
title('')

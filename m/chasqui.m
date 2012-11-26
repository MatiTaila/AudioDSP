function chasqui(fs,t,BPM,texto)

% 
% [train] = beat_train_template(60*fs/BPM,1,fs*t,10);
% [B,A]   = butter(1,0.01,'low');
% audio = filtfilt(B,A,train);
% wavwrite(abs(audio*(1/0.02)),fs,texto);
% 
% plot(abs(audio*(1/0.02)))

[train] = beat_train_template(60*fs/BPM,1,fs*t,10);
% load click sound
[click,fs_click] = wavread('..\audio\31-sticks.wav');
if fs_click ~= fs
    fprintf('===========================================\nOJOOOO!!! FRECUENCIAS DE MUESTREO DISTINTAS\n===========================================\n')
end

tracked_beats = conv(train,click,'same');
tracked_beats = tracked_beats/max(abs(tracked_beats));
wavwrite(tracked_beats,fs,texto);

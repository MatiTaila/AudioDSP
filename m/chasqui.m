function chasqui(fs,t,BPM,texto,opt)
% -------------------------------------------------------------------------
% [train] = beat_train_template(60*fs/BPM,1,fs*t,10);
% -------------------------------------------------------------------------
% Example: chasqui(44100,10,180,'180BPM',opt);
% -------------------------------------------------------------------------
% [B,A]   = butter(1,0.01,'low');
% audio = filtfilt(B,A,train);
% wavwrite(abs(audio*(1/0.02)),fs,texto);
% plot(abs(audio*(1/0.02)))
% -------------------------------------------------------------------------

[train] = beat_train_template(60*fs/BPM,fs*t,10);
% load click sound
if opt
    [click,fs_click] = wavread('../../../../matlab/audio/beatroot/audio/53-cymbell.wav');
else
    [click,fs_click] = wavread('..\audio\31-sticks.wav');
end
if fs_click ~= fs
    fprintf('===========================================\nOJOOOO!!! FRECUENCIAS DE MUESTREO DISTINTAS\n===========================================\n')
end

tracked_beats = conv(train,click,'same');
tracked_beats = tracked_beats/max(abs(tracked_beats));
wavwrite(tracked_beats,fs,texto);

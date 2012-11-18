function chasqui(fs,t,BPM,texto)


[train] = beat_train_template(60*fs/BPM,1,fs*t,10);
[B,A]   = butter(1,0.01,'low');
audio = filtfilt(B,A,train);
wavwrite(abs(audio*(1/0.02)),fs,texto);

plot(abs(audio*(1/0.02)))
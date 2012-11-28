function bpm = frames2bpm(input,fs,n_hop)
bpm = 60/(input/(fs/n_hop));
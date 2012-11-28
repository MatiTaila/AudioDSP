function frames = bpm2frames(input,fs,n_hop)
frames = fs/n_hop/(input/60);
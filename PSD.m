function [pxx,f] = PSD(dati,fs,percentage_overlap,segment_duration)

%Numero di punti di una finestra
n_window = segment_duration*fs;
%Finestra di Hann
window = hann(n_window);
%Sovrapposizione tra segmenti
noverlap = round(n_window * percentage_overlap / 100);
nfft = n_window;
[pxx, f] = pwelch(double(dati), window, noverlap, nfft, fs);
end
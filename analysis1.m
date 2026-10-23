
%SOLO SEGNALE IN  uV SENZA FILTRI

%clear;clc;close all;

%PER L'USO DELLE FUNZIONI

h = heatmaps();
p = plots();

%VISUALIZZAZIONE DEL FILE H5

nomeFile = '../Dati/file0109.h5';
%Dico quale dataset sto guardando
date = h5readatt(nomeFile,'/Data','Date');
fprintf("Dati del: %s\n",date);

%Informazioni generali
%h5disp(nomeFile);
%fprintf("\n");

%Infromazioni sul mapping canale-posizione elettrodo nel MEA
info = h5read(nomeFile,'/Data/Recording_0/AnalogStream/Stream_0/InfoChannel');
T = table(info.ChannelID,info.RowIndex, string(info.Label), 'VariableNames', {'ChannelID','RowIndex', 'Label'});
%disp(T)
%fprintf("\n");

%Apertura del dataset
pathDati = '/Data/Recording_0/AnalogStream/Stream_0/ChannelData';
dati = h5read(nomeFile,pathDati);


%CONTROLLO DEI DATI

%Check sulle dimensioni della matrice di dati e del corretto orientamento
[n_samples,n_channels] = size(dati);
if(n_samples==60 && n_channels>60)
    fprintf("Riorentamento corretto dellla matrice di dati\n");
    dati = dati.';
    [n_samples,n_channels] = size(dati);
end
fprintf("Numero di campioni presi: %d\n",n_samples);
fprintf("Numero di canali totali: %d\n",n_channels);

%Check sulla validità dei dati
assert(all(isfinite(dati(:))),'Presenti Not Number e/o dati non finiti, trattarli prima di proseguire');

%Check sulla durata della registrazione
%Dal laboratorio è stata data la frequenza di campionamento, nei metadati %non si trova (?)
fs = 10000;
duration = (n_samples/fs)/60;
fprintf("Durata della registrazione dell'attività elettrica dell' organoide: %f minuti\n",duration);

%Vettore con gli istanti di campionamento (timestamps)
timestamps = (0:n_samples-1) / fs;


dati_noOff = dati - int32(round(mean(dati,1)));
time = timestamps>=0 & timestamps<0.1;
%segnale in microVolt
dati_uV = dati_noOff * 8670 * 1e-6;

clear dati;
clear dati_noOff;


%ANALISI IN FREQUENZA DE SEGNALE IN uV

s = FFT(dati_uV,n_samples);
% Asse delle frequenze
frequencies = (0:floor(n_samples/2))*fs/n_samples;
%p.frequenze(frequencies,s,'Ampiezza Spettrale','Spettro delle frequenze');
%h.frequenze(frequencies,s,'Ampiezza Spettrale (uV)','Heatmap della FFT del segnale')

%PSD DEL SEGNALE IN uV

%psd calcolata, vettore freq corrispondenti
[psd_uV,frequencies_uVPsd] = PSD(dati_uV,fs,50,1);
%p.psd(frequencies_uVPsd,psd_uV,'PSD (dB re 1 uV^2/Hz)','Densità spettrale di potenza del segnale in uV');
%index = frequencies_uVPsd<=150;
%h.frequenze(frequencies_uVPsd(index),10*log10(psd_uV(index,:)),'PSD (dB re 1 uV^2/Hz)','Heatmap della PSD del segnale')

%-------------------------- PLOT DEI SEGNALI IN VARI CANALI NON SECUTIVI ------------------------------------------
time = timestamps>= 0 & timestamps<0.1;
canali = [24,26,54,56];
p.tempi(time,dati_uV,timestamps,canali);

%{
canale = 30;
segnale = double(dati_uV(:,canale));
figure;
spectrogram(segnale, hann(10000), 5000, 10000, fs, 'yaxis');
title(sprintf('Spettrogramma - Canale %d',canale));
colorbar;
%}


%{
% Calcolo dello spettrogramma
[S,F,T,P] = spectrogram(segnale, ...
    hann(10000), 5000, 10000, fs);
% Visualizzazione della PSD senza conversione in dB
figure;
imagesc(T, F, P);
axis xy;
xlabel('Tempo (s)');
ylabel('Frequenza (Hz)');
title(sprintf('Spettrogramma - Canale %d',canale));
cb = colorbar;
cb.Label.String = 'PSD (\muV^2/Hz)';
%}


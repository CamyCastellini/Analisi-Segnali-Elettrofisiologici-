clear;clc;close all;

%PER L'USO DELLE FUNZIONI
h = heatmaps();
p = plots();

%VISUALIZZAZIONE DEL FILE H5
nomeFile = '../Dati/file0109.h5';
%Dico quale dataset sto guardando
date = h5readatt(nomeFile,'/Data','Date');
fprintf("Dati del: %s\n",date);

%Informazioni generali
h5disp(nomeFile);
%fprintf("\n");

%Infromazioni sul mapping canale-posizione elettrodo nel MEA
info = h5read(nomeFile,'/Data/Recording_0/AnalogStream/Stream_0/InfoChannel');
T = table(info.ChannelID,info.RowIndex, string(info.Label), 'VariableNames', {'ChannelID','RowIndex', 'Label'});
%disp(T)
%fprintf("\n");
%Infomrazio sui parametri di conversione 
% Visualizza i parametri
fprintf('ADZero:');
disp(info.ADZero(1));
fprintf('ConversionFactor:');
disp(info.ConversionFactor(1));
fprintf('Exponent:');
disp(info.Exponent(1));
fprintf('Unit:')
disp(info.Unit(1));

%Apertura del dataset
pathDati = '/Data/Recording_0/AnalogStream/Stream_0/ChannelData';
dati = h5read(nomeFile,pathDati);


%STEP 1 ---> CONTROLLO DEI DATI

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


%STEP 2 -->CONTROLLO DI QUALITà DEI SEGNALI

%Calcolo un vettore riga con n_channels colonne che contiene le mediane delle colonne
%Le mediane fungono da offset--> togliendolo rialliniamo la matrice "intorno allo zero"
dati_noOff = dati - int32(round(mean(dati,1)));
time = timestamps>=0 & timestamps<0.02;
%FIGURA 1
%p.tempo(time,dati_noOff,timestamps,'Ampiezza','Tracce grezze');

%Visualizzazione del segnale in microVolt
dati_uV = dati_noOff * 8670 * 1e-6;
%FIGURA 2
p.tempo(time,dati_uV,timestamps,'Ampiezza (uV)','Tracce in micro Volt');
h.tempo(timestamps(time),dati_uV(time,:),'Heatmap temporale del segnale');

%ORA CHE HO IL SEGNALE CON L'UNITà DI MISURA GIUSTO CANCELLO LE COPIE PRIMA
clear dati;
clear dati_noOff;

%Calcolo del RMS per ogni canale
rms = sqrt(mean(dati_uV.^2,1)); 
%FIGURA 3
%{
figure;
bar(rms);
xlabel('Canale');
ylabel('RMS '); 
%Vediamo gli elettrodi 15 (min RSM) e il 30 (max RSM)
%FIGURA 4
figure
plot(timestamps(time),dati_uV(time,15));
ylabel(sprintf('Canale %d',15));
%FIGURA 5
figure
plot(timestamps(time),dati_uV(time,30));
ylabel(sprintf('Canale %d',30));
%}

%{
%AMPIEZZE MASSIME E MINIME
valMax = max(abs(dati_uV), [], 'all');
[rMax, cMax] = find(abs(dati_uV) == valMax);
valMin = min(abs(dati_uV), [], 'all');
[rMin, cMin] = find(abs(dati_uV) == valMin);

%Stampa Massimi--->TUTTI RAGGIUNGONO L'AMPIEZZA MINIMA
maxPerCanale = max(abs(dati_uV), [], 1);
valMax = max(maxPerCanale);
canaliMax = find(maxPerCanale == valMax);
fprintf('Ampiezza massima: %.3f\n', valMax);
fprintf('Canali: ')
for i = 1:numel(canaliMax)
    fprintf('%d ', canaliMax(i));
end
fprintf('\n');

%Stampa Minimi
minPerCanale = min(abs(dati_uV), [], 1);   
valMin = min(minPerCanale);                  
canaliMin = find(minPerCanale == valMin);     
fprintf('Ampiezza minima: %.3f\n', valMin);
fprintf('Canali: ')
for i = 1:numel(canaliMin)
    fprintf('%d ',canaliMin(i));
end
fprintf('\n');
%}

%ANALISI IN FREQUENZA DE SEGNALE IN uV
s = FFT(dati_uV,n_samples);
% Asse delle frequenze
frequencies = (0:floor(n_samples/2))*fs/n_samples;
%FIGURA 6
%p.frequenze(frequencies,s,'Ampiezza Spettrale','Spettro delle frequenze');
%h.frequenze(frequencies,s,'Ampiezza Spettrale (uV)','Heatmap della FFT del segnale')

%PSD DEL SEGNALE IN uV
%psd calcolata, vettore freq corrispondenti
[psd_uV,frequencies_uVPsd] = PSD(dati_uV,fs,50,1);
%FIGURA 7
%p.psd(frequencies_uVPsd,psd_uV,'PSD (dB re 1 uV^2/Hz)','Densità spettrale di potenza del segnale in uV');
%index = frequencies_uVPsd<=150;
%h.frequenze(frequencies_uVPsd(index),10*log10(psd_uV(index,:)),'PSD (dB re 1 uV^2/Hz)','Heatmap della PSD del segnale')

%STEP 3--> FILTRAGGIO 

%SPIKING EXTRA-CELLULARE ---> FILTRO PASSA-BANDA 300-3000
%4--> lv di attenuazione delle frequenze fuori dal range
%normalizzazione delle frequenze rispetto alla frequenza di Nyquist
[b,a] = butter(4, [300 3000]/(fs/2), 'bandpass');
dati_uV = double(dati_uV);
dati_filtrati_se = filtfilt(b, a, dati_uV);
%FIGURA 8
p.tempo(time,dati_filtrati_se,timestamps,'Ampiezza (uV)','Spiking Extracellulare');
h.tempo(timestamps(time),dati_filtrati_se(time,:),'Heatmap temporale dello spiking Extra-Cellulare');
%ANALISI IN FREQUENZA 
sf_se = FFT(dati_filtrati_se,n_samples);
%FIGURA 9
%p.frequenze(frequencies,sf_se,'Ampiezza Spettrale','Spettro delle frequenze dello Spiking Extra-cellulare');
%index = frequencies>=300 & frequencies<=3000;
%h.frequenze(frequencies(index),sf_se(index,:),'Ampiezza Spettrale','Heatmap dello spettro delle frequenze dello Spiking Extra-Cellulare');

%LOCAL FIELD POTENTIAL ---> FILTRO PASSA BASSO <1000
%frequenza di taglio
f_cut = 1000;
[b, a] = butter(4, f_cut/(fs/2), 'low');
dati_filtrati_lfp = filtfilt(b, a, double(dati_uV));
%FIGURA 10
%p.tempo(time,dati_filtrati_lfp,timestamps,'Ampiezza (uV)','Local Field Potential');
%h.tempo(timestamps(time),dati_filtrati_lfp(time,:),'Heatmap temporale del Local Field Potential')
%ANALISI IN FREQUENZA 
sf_lfp = FFT(dati_filtrati_lfp,n_samples);
%FIGURA 11
%p.frequenze(frequencies,sf_lfp,'Ampiezza Spettrale','Spettro delle frequenze del Local Field Potential');
index = frequencies<=1000;
%h.frequenze(frequencies(index),sf_lfp(index,:),'Ampiezza Spettrale','Heatmap delle frequenze del Local Field Potential');

%PSD al LFP --> capire a quali frequenze è distribuita la potenza del segnale elettrico, 
% di individuare eventuali oscillazioni dominanti e quali componenti in frequenza contribuiscono maggiormente alla varianza del segnale.
%La sinusoide con ampiezza maggiore ha potenza maggiore (P=A^2/2)
%l risultato ti dice quanta potenza per unità di frequenza è presente a ciascuna frequenza.
%Dopo per distingure il fondo aperiodico dai vari picchi --> FOOOF

%psd calcolata, vettore freq corrispondenti
[psd_lfp, frequenciesPsd] = PSD(dati_filtrati_lfp,fs,50,2);
%FIGURA 12
%p.psd(frequenciesPsd,psd_lfp,'PSD (dB re 1 uV^2/Hz)','Densità spettrale di potenza LFP');
index = frequenciesPsd<=1000;
%h.frequenze(frequenciesPsd(index),10*log10(psd_lfp(index,:)),'PSD (dB re 1 uV^2/Hz)','Heatmap della PSD del Local Field Potential');

%{
fc = 100;
[b, a] = butter(4, fc/(fs/2), 'low');
dati_sonno = filtfilt(b, a, double(dati_uV));
% Isolamento della componente a bassa frequenza associata al sonno
timeSonno = timestamps >= 0 & timestamps < 60;
p.tempo(timeSonno,dati_sonno,timestamps,'Ampiezza (uV)','Segnale neel tempo filtrato passa-basso 100Hz');
h.tempo(timestamps(timeSonno), dati_sonno(timeSonno,:),'Heatmap temporale del segnale filtrato passa-basso 100Hz');
% Analisi spettrale della componente associata al sonno
spettroSonno = FFT(dati_sonno,n_samples);
indexSonno = frequencies <= fc;
frequenze100 = frequencies(indexSonno);
h.frequenze(frequencies(indexSonno),spettroSonno(indexSonno,:),'Ampiezza Spettrale (uV)','Heatmap dello spettro durante il sonno');
p.frequenze(frequencies,spettroSonno,'Ampiezza Spettrale','FFT del segnale filtrato passa-basso 100Hz');
% PSD della componente a bassa frequenza associata al sonno
[psdSonno, frequenciesPsdSonno] = PSD(dati_sonno, fs, 50, 2);
indexSonnoPsd = frequenciesPsdSonno <= fc;
h.frequenze(frequenciesPsdSonno(indexSonnoPsd),10*log10(psdSonno(indexSonnoPsd,:)),  'PSD (dB re 1 uV^2/Hz)', 'Heatmap PSD segnale filtrato passa-basso 100Hz');
%}


%% =========================================================
% STEP 4 - SPIKE DETECTION SU TUTTI I 60 CANALI
% =========================================================

% Variabili gia' disponibili nel tuo script:
% dati_filtrati_se -> matrice campioni x canali, in uV
% timestamps       -> vettore dei tempi in secondi
% fs               -> frequenza di campionamento
% n_samples        -> numero di campioni
% n_channels       -> numero di canali

% Parametri della spike detection
fattore_soglia = 5;
distanza_minima = round(0.001 * fs);  % 1 ms

% Parametri per estrarre le forme d'onda
pre = round(0.001 * fs);    % 1 ms prima del picco
post = round(0.002 * fs);   % 2 ms dopo il picco

% Seleziona i canali per il confronto delle forme d'onda
% e per il secondo raster plot
canali_scelti = [24 26 27 38];    % MODIFICA QUI

assert(all(canali_scelti >= 1 & ...
           canali_scelti <= n_channels & ...
           mod(canali_scelti,1) == 0), ...
       'Inserisci numeri di canale validi.');

% Risultati della spike detection
tempi_spike = cell(1, n_channels);
indici_spike = cell(1, n_channels);
numero_spike = zeros(1, n_channels);
rumore_spike = zeros(1, n_channels);
soglie_spike = zeros(1, n_channels);

%% RILEVAZIONE DEGLI SPIKE

for ch = 1:n_channels

    % Segnale gia' filtrato tra 300 e 3000 Hz
    segnale_spike = double(dati_filtrati_se(:,ch));

    % Stima robusta del rumore
    sigma_spike = median(abs(segnale_spike - ...
                  median(segnale_spike))) / 0.6745;

    % Soglia negativa
    soglia_spike = -fattore_soglia * sigma_spike;

    rumore_spike(ch) = sigma_spike;
    soglie_spike(ch) = soglia_spike;

    % Inversione per trovare i minimi negativi
    segnale_invertito = -segnale_spike;

    % Evita warning se nessun campione supera la soglia
    if isfinite(sigma_spike) && sigma_spike > 0 && ...
            any(segnale_invertito > abs(soglia_spike))

        [~, locs_spike] = findpeaks( ...
            segnale_invertito, ...
            'MinPeakHeight', abs(soglia_spike), ...
            'MinPeakDistance', distanza_minima);

    else
        locs_spike = [];
    end

    % Conserva gli indici dei campioni e i tempi
    indici_spike{ch} = locs_spike(:);
    tempi_spike{ch} = (locs_spike(:) - 1) / fs;

    numero_spike(ch) = numel(locs_spike);

end

%% TABELLA RIASSUNTIVA

durata_secondi = n_samples / fs;

frequenza_spike = numero_spike / durata_secondi;

tabella_spike = table( ...
    (1:n_channels)', ...
    numero_spike', ...
    frequenza_spike', ...
    rumore_spike', ...
    soglie_spike', ...
    'VariableNames', ...
    {'Canale','NumeroEventi','EventiAlSecondo', ...
     'Rumore_uV','Soglia_uV'});

disp(tabella_spike);

%% =========================================================
% GRAFICO 1 - RASTER PLOT DI TUTTI I CANALI
% =========================================================

figure('Color','w', ...
       'Position',[100 100 1200 700]);

ax = gca;
set(ax,'Color','w','XColor','k','YColor','k');

hold on;

for ch = 1:n_channels

    t = tempi_spike{ch};

    if ~isempty(t)
        plot(t, ch * ones(size(t)), ...
            'k.', 'MarkerSize', 6);
    end

end

xlabel('Tempo (s)');
ylabel('Canale di acquisizione');
title('Raster plot - Tutti i canali');

xlim([0 durata_secondi]);
ylim([0.5 n_channels + 0.5]);

yticks(1:2:n_channels);

grid on;
box on;

exportgraphics(gcf, 'raster_tutti_canali.png', ...
    'Resolution', 300);

%% =========================================================
% GRAFICO 2 - RASTER PLOT DEI CANALI SELEZIONATI
% =========================================================

figure('Color','w', ...
       'Position',[100 100 1200 400]);

ax = gca;
set(ax,'Color','w','XColor','k','YColor','k');

hold on;

colori_raster = lines(numel(canali_scelti));

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);
    t = tempi_spike{ch};

    if ~isempty(t)
        plot(t, i * ones(size(t)), '|', ...
            'Color', colori_raster(i,:), ...
            'MarkerSize', 12, ...
            'LineWidth', 1.2);
    end

end

xlabel('Tempo (s)');
ylabel('Canale di acquisizione');
title('Raster plot - Canali selezionati');

yticks(1:numel(canali_scelti));
yticklabels(string(canali_scelti));

ylim([0.5 numel(canali_scelti) + 0.5]);
xlim([0 durata_secondi]);

grid on;
box on;

exportgraphics(gcf, 'raster_canali_scelti.png', ...
    'Resolution', 300);

%% =========================================================
% GRAFICO 3 - ZOOM TEMPORALE DEI CANALI SELEZIONATI
% =========================================================

% Modifica l'intervallo temporale a piacere
t_inizio = 56,5;
t_fine = 59;

assert(t_inizio >= 0 && t_fine > t_inizio && ...
       t_fine <= durata_secondi, ...
       'Intervallo temporale non valido.');

figure('Color','w', ...
       'Position',[100 100 1200 400]);

ax = gca;
set(ax,'Color','w','XColor','k','YColor','k');

hold on;

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);
    t = tempi_spike{ch};

    % Seleziona gli spike nell'intervallo scelto
    t = t(t >= t_inizio & t <= t_fine);

    if ~isempty(t)
        plot(t, i * ones(size(t)), '|', ...
            'Color', colori_raster(i,:), ...
            'MarkerSize', 12, ...
            'LineWidth', 1.2);
    end

end

xlabel('Tempo (s)');
ylabel('Canale di acquisizione');

title(sprintf('Raster plot - Intervallo %.1f-%.1f s', ...
    t_inizio, t_fine));

yticks(1:numel(canali_scelti));
yticklabels(string(canali_scelti));

ylim([0.5 numel(canali_scelti) + 0.5]);
xlim([t_inizio t_fine]);

grid on;
box on;

exportgraphics(gcf, 'raster_zoom.png', ...
    'Resolution', 300);

%% =========================================================
% GRAFICO 4 - FORME D'ONDA DEI CANALI SELEZIONATI
% =========================================================

% Tempo relativo al picco, espresso in millisecondi
tempo_forma = (-pre:post)' / fs * 1000;

% Conserva le forme d'onda per ogni canale scelto
forme_onda = cell(1, numel(canali_scelti));

figure('Color','w', ...
       'Position',[100 100 1200 650]);

n_selezionati = numel(canali_scelti);

for i = 1:n_selezionati

    ch = canali_scelti(i);

    segnale_spike = double(dati_filtrati_se(:,ch));

    locs_spike = indici_spike{ch};

    % Esclude eventi troppo vicini ai bordi
    locs_validi = locs_spike( ...
        locs_spike > pre & ...
        locs_spike + post <= n_samples);

    n_eventi = numel(locs_validi);

    onde = zeros(pre + post + 1, n_eventi);

    for k = 1:n_eventi

        indici = (locs_validi(k)-pre): ...
                 (locs_validi(k)+post);

        onde(:,k) = segnale_spike(indici);

    end

    forme_onda{i} = onde;

    subplot(ceil(sqrt(n_selezionati)), ...
            ceil(n_selezionati / ceil(sqrt(n_selezionati))), i);

    hold on;

    if ~isempty(onde)

        % Mostra fino a 200 forme d'onda
        n_mostrate = min(size(onde,2), 200);

        indici_mostrati = round(linspace( ...
            1, size(onde,2), n_mostrate));

        plot(tempo_forma, onde(:,indici_mostrati), ...
            'Color',[0.65 0.75 0.9]);

        % Forma d'onda mediana
        plot(tempo_forma, median(onde,2), ...
            'r','LineWidth',2);

    end

    xline(0,'--k');

    xlabel('Tempo rispetto al picco (ms)');
    ylabel('Ampiezza (\muV)');

    title(sprintf('Canale %d - %d eventi', ...
        ch, n_eventi));

    grid on;
    box on;

end

sgtitle('Forme d''onda dei possibili spike extracellulari');

exportgraphics(gcf, 'forme_onda_canali_scelti.png', ...
    'Resolution', 300);

%% =========================================================
% RISULTATI DEI CANALI SELEZIONATI
% =========================================================

fprintf('\n===== CANALI SELEZIONATI =====\n');

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    fprintf('\nCanale %d\n', ch);
    fprintf('Rumore: %.3f uV\n', rumore_spike(ch));
    fprintf('Soglia: %.3f uV\n', soglie_spike(ch));
    fprintf('Eventi rilevati: %d\n', numero_spike(ch));
    fprintf('Frequenza: %.4f eventi/s\n', ...
        frequenza_spike(ch));

end

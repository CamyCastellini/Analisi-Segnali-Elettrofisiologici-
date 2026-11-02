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
%h5disp(nomeFile);
%fprintf("\n");

%Infromazioni sul mapping canale-posizione elettrodo nel MEA
info = h5read(nomeFile,'/Data/Recording_0/AnalogStream/Stream_0/InfoChannel');
T = table(info.ChannelID,info.RowIndex, string(info.Label), 'VariableNames', {'ChannelID','RowIndex', 'Label'});
%disp(T)
%fprintf("\n");
%Infomrazio sui parametri di conversione 
% Visualizza i parametri
%{
fprintf('ADZero:');
disp(info.ADZero(1));
fprintf('ConversionFactor:');
disp(info.ConversionFactor(1));
fprintf('Exponent:');
disp(info.Exponent(1));
fprintf('Unit:')
disp(info.Unit(1));
%}
%Apertura del dataset
pathDati = '/Data/Recording_0/AnalogStream/Stream_0/ChannelData';
dati = h5read(nomeFile,pathDati);

%================================
%STEP 1 ---> CONTROLLO DEI DATI
%================================

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

%================================
%STEP 2 -->CONTROLLO DI QUALITà DEI SEGNALI
%================================

%Calcolo un vettore riga con n_channels colonne che contiene le mediane delle colonne
%Le mediane fungono da offset--> togliendolo rialliniamo la matrice "intorno allo zero"
dati_noOff = dati - int32(round(mean(dati,1)));
time = timestamps>=0 & timestamps<0.02;
%FIGURA 1
%p.tempo(time,dati_noOff,timestamps,'Ampiezza','Tracce grezze');

%Visualizzazione del segnale in microVolt
dati_uV = dati_noOff * 8670 * 1e-6;
%FIGURA 2
%p.tempo(time,dati_uV,timestamps,'Ampiezza (uV)','Tracce in micro Volt');
%h.tempo(timestamps(time),dati_uV(time,:),'Heatmap temporale del segnale');

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

%================================
%STEP 3--> FILTRAGGIO 
%================================

%SPIKING EXTRA-CELLULARE ---> FILTRO PASSA-BANDA 300-3000
%4--> lv di attenuazione delle frequenze fuori dal range
%normalizzazione delle frequenze rispetto alla frequenza di Nyquist
[b,a] = butter(4, [300 3000]/(fs/2), 'bandpass');
dati_uV = double(dati_uV);
dati_filtrati_se = filtfilt(b, a, dati_uV);
%FIGURA 8
%p.tempo(time,dati_filtrati_se,timestamps,'Ampiezza (uV)','Spiking Extracellulare');
%h.tempo(timestamps(time),dati_filtrati_se(time,:),'Heatmap temporale dello spiking Extra-Cellulare');
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


%% ========================================================================
% STEP 4 - RILEVAMENTO BIPOLARE DEI CANDIDATI SPIKE EXTRACELLULARI
% ========================================================================
% Il segnale e' gia' stato filtrato nella banda 300-3000 Hz.
%
% Si cercano sia picchi negativi sia picchi positivi, utilizzando
% il valore assoluto del segnale e una soglia pari a 5 sigma.
%
% La stima sigma e' ricavata dalla MAD:
% sigma = median(abs(x - median(x))) / 0.6745
%
% Gli eventi sono considerati candidati extracellulari, senza
% attribuire automaticamente loro un'origine neuronale.
%
% NOTA:
% La distanza minima di 1 ms limita il doppio conteggio delle
% componenti bifasiche, ma non lo elimina necessariamente.

% -------------------- PARAMETRI -----------------------------------------

fattore_soglia = 5;

distanza_minima = round(0.001 * fs);  % 1 ms

pre = round(0.001 * fs);             % 1 ms prima del picco
post = round(0.002 * fs);            % 2 ms dopo il picco

canali_scelti = [27 38];

durata_secondi = n_samples / fs;

tempo_forma = (-pre:post)' / fs * 1000;

% -------------------- INIZIALIZZAZIONE ----------------------------------

indici_spike = cell(n_channels,1);
tempi_spike = cell(n_channels,1);
forme_onda = cell(n_channels,1);

polarita_spike = cell(n_channels,1);

numero_spike = zeros(n_channels,1);
numero_spike_positivi = zeros(n_channels,1);
numero_spike_negativi = zeros(n_channels,1);

rumore_spike = zeros(n_channels,1);
soglie_spike = zeros(n_channels,1);

% -------------------- RILEVAMENTO ---------------------------------------

for ch = 1:n_channels

    x = double(dati_filtrati_se(:,ch));
    x = x(:);

    % Stima robusta del livello di fondo
    sigma = median(abs(x - median(x))) / 0.6745;

    rumore_spike(ch) = sigma;
    soglie_spike(ch) = fattore_soglia * sigma;

    if ~isfinite(sigma) || sigma <= 0
        continue;
    end

    % Rilevamento sul valore assoluto:
    % permette di identificare eventi di entrambe le polarita'.
    [~,locs] = findpeaks(abs(x), ...
        'MinPeakHeight',fattore_soglia*sigma, ...
        'MinPeakDistance',distanza_minima);

    locs = locs(:);

    indici_spike{ch} = locs;
    tempi_spike{ch} = (locs-1)/fs;

    % Polarita' del picco dominante
    pol = sign(x(locs));
    polarita_spike{ch} = pol;

    numero_spike(ch) = numel(locs);
    numero_spike_positivi(ch) = sum(pol > 0);
    numero_spike_negativi(ch) = sum(pol < 0);

    % Estrazione delle forme d'onda attorno al picco dominante
    locs_validi = locs(locs > pre & locs+post <= n_samples);

    onde = nan(pre+post+1,numel(locs_validi));

    for k = 1:numel(locs_validi)

        idx = locs_validi(k);

        onde(:,k) = x(idx-pre:idx+post);

    end

    forme_onda{ch} = onde;

end

% -------------------- TABELLA RISULTATI ---------------------------------

frequenza_spike = numero_spike / durata_secondi;

tabella_spike = table( ...
    (1:n_channels)', ...
    numero_spike, ...
    numero_spike_positivi, ...
    numero_spike_negativi, ...
    frequenza_spike, ...
    rumore_spike, ...
    soglie_spike, ...
    'VariableNames',{ ...
    'Canale', ...
    'CandidatiTotali', ...
    'CandidatiPositivi', ...
    'CandidatiNegativi', ...
    'EventiAlSecondo', ...
    'Rumore_uV', ...
    'SogliaAssoluta_uV'});

fprintf('\n===== RILEVAMENTO BIPOLARE DEI CANDIDATI =====\n');
fprintf('Soglia: +/- %.1f sigma\n',fattore_soglia);
fprintf('Distanza minima: %.2f ms\n',distanza_minima/fs*1000);

disp(tabella_spike);

fprintf('\n===== CANALI SELEZIONATI =====\n');

disp(tabella_spike(canali_scelti,:));


%% ========================================================================
% STEP 4 - RASTER PLOT DI TUTTI I CANALI
% ========================================================================
% Ogni punto rappresenta un candidato rilevato.
% Rosso: picco dominante negativo.
% Blu: picco dominante positivo.

figure('Color','w', ...
    'Name','Distribuzione temporale dei candidati spike');

hold on;

for ch = 1:n_channels

    t = tempi_spike{ch};
    p = polarita_spike{ch};

    if isempty(t)
        continue;
    end

    plot(t(p<0),ch*ones(sum(p<0),1), ...
        '.','Color',[0.85 0.20 0.20], ...
        'MarkerSize',7);

    plot(t(p>0),ch*ones(sum(p>0),1), ...
        '.','Color',[0.15 0.40 0.85], ...
        'MarkerSize',7);

end

xlabel('Tempo (s)');
ylabel('Canale di acquisizione');

title('Distribuzione temporale dei candidati extracellulari');

xlim([0 durata_secondi]);
ylim([0.5 n_channels+0.5]);

grid on;

set(gca,'Color','w', ...
    'XColor','k','YColor','k', ...
    'GridColor',[0.8 0.8 0.8], ...
    'GridAlpha',0.25);

hold off;


%% ========================================================================
% STEP 4 - RASTER PLOT DEI CANALI SELEZIONATI
% ========================================================================
% Si visualizza l'intera registrazione, senza zoom temporale.

figure('Color','w', ...
    'Name','Candidati extracellulari nei canali selezionati');

hold on;

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    t = tempi_spike{ch};
    p = polarita_spike{ch};

    for k = 1:numel(t)

        if p(k) < 0
            colore = [0.85 0.20 0.20];
        else
            colore = [0.15 0.40 0.85];
        end

        line([t(k) t(k)], ...
            [i-0.3 i+0.3], ...
            'Color',colore, ...
            'LineWidth',1.2);

    end

end

xlabel('Tempo (s)');
ylabel('Canale di acquisizione');

title('Distribuzione temporale dei candidati - Canali selezionati');

yticks(1:numel(canali_scelti));
yticklabels(string(canali_scelti));

ylim([0.5 numel(canali_scelti)+0.5]);
xlim([0 durata_secondi]);

grid on;

set(gca,'Color','w', ...
    'XColor','k','YColor','k', ...
    'GridColor',[0.8 0.8 0.8], ...
    'GridAlpha',0.25);

hold off;


%% ========================================================================
% STEP 4 - MORFOLOGIA DEI CANDIDATI
% ========================================================================
% Le forme d'onda vengono separate per polarita'.
% In questo modo non si mescolano direttamente picchi
% positivi e negativi nella stessa forma mediana.

figure('Color','w', ...
    'Name','Morfologia dei candidati extracellulari');

tiledlayout(2,2);

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    nexttile;
    hold on;

    x = double(dati_filtrati_se(:,ch));
    x = x(:);

    locs = indici_spike{ch};

    locs = locs(locs > pre & locs+post <= n_samples);

    onde_pos = [];
    onde_neg = [];

    for k = 1:numel(locs)

        idx = locs(k);

        onda = x(idx-pre:idx+post);

        if x(idx) > 0
            onde_pos(:,end+1) = onda;
        else
            onde_neg(:,end+1) = onda;
        end

    end

    if ~isempty(onde_neg)

        plot(tempo_forma,onde_neg, ...
            'Color',[1.00 0.75 0.75]);

        plot(tempo_forma,median(onde_neg,2), ...
            'r','LineWidth',2);

    end

    if ~isempty(onde_pos)

        plot(tempo_forma,onde_pos, ...
            'Color',[0.70 0.80 1.00]);

        plot(tempo_forma,median(onde_pos,2), ...
            'b','LineWidth',2);

    end

    xline(0,'k:');

    xlabel('Tempo rispetto al picco (ms)');
    ylabel('Ampiezza (\muV)');

    title(sprintf('Canale %d - %d candidati', ...
        ch,numero_spike(ch)));

    grid on;

    set(gca,'Color','w', ...
        'XColor','k','YColor','k', ...
        'GridColor',[0.8 0.8 0.8], ...
        'GridAlpha',0.25);

    hold off;

end

sgtitle('Morfologia dei candidati extracellulari per polarita''');


%% ========================================================================
% STEP 4B - CARATTERIZZAZIONE DEI CANDIDATI
% ========================================================================
% Per ogni evento si calcolano:
%
% 1. Ampiezza assoluta del picco dominante
% 2. Rapporto segnale-rumore
% 3. Larghezza a meta' ampiezza
% 4. Intervallo tra eventi consecutivi (ISI)
%
% La larghezza viene calcolata rispetto alla polarita'
% del picco dominante, usando l'interpolazione lineare
% degli attraversamenti della meta' ampiezza.

ampiezze_4b = cell(n_channels,1);
snr_4b = cell(n_channels,1);
larghezze_4b = cell(n_channels,1);
isi_4b = cell(n_channels,1);

ampiezza_mediana_4b = nan(n_channels,1);
snr_mediano_4b = nan(n_channels,1);
larghezza_mediana_4b = nan(n_channels,1);
isi_mediano_4b = nan(n_channels,1);

larghezze_valide_4b = zeros(n_channels,1);

for ch = 1:n_channels

    x = double(dati_filtrati_se(:,ch));
    x = x(:);

    locs = indici_spike{ch};

    n = numel(locs);

    amp = nan(n,1);
    snr = nan(n,1);
    lar = nan(n,1);

    for k = 1:n

        idx = locs(k);

        amp(k) = abs(x(idx));

        if rumore_spike(ch) > 0
            snr(k) = amp(k)/rumore_spike(ch);
        end

        if idx <= pre || idx+post > n_samples
            continue;
        end

        onda = x(idx-pre:idx+post);

        % Si orienta la forma d'onda affinche' il picco
        % dominante sia sempre positivo per il calcolo.
        segno = sign(x(idx));

        onda_orientata = segno * onda;

        centro = pre+1;

        % Baseline locale robusta
        baseline = median(onda_orientata);

        altezza = onda_orientata(centro)-baseline;

        if altezza <= 0
            continue;
        end

        meta = baseline + altezza/2;

        % Attraversamento sinistro
        sx = find(onda_orientata(1:centro-1) < meta,1,'last');

        % Attraversamento destro
        dx_rel = find(onda_orientata(centro+1:end) < meta,1,'first');

        if isempty(sx) || isempty(dx_rel)
            continue;
        end

        dx = centro+dx_rel;

        % Interpolazione lineare a sinistra
        y1 = onda_orientata(sx);
        y2 = onda_orientata(sx+1);

        if y2 == y1
            continue;
        end

        t_sx = sx + (meta-y1)/(y2-y1);

        % Interpolazione lineare a destra
        y1 = onda_orientata(dx-1);
        y2 = onda_orientata(dx);

        if y2 == y1
            continue;
        end

        t_dx = (dx-1) + (meta-y1)/(y2-y1);

        lar(k) = (t_dx-t_sx)/fs*1000;

    end

    ampiezze_4b{ch} = amp;
    snr_4b{ch} = snr;
    larghezze_4b{ch} = lar;

    isi = diff(tempi_spike{ch})*1000;

    isi_4b{ch} = isi;

    if any(isfinite(amp))
        ampiezza_mediana_4b(ch) = median(amp(isfinite(amp)));
    end

    if any(isfinite(snr))
        snr_mediano_4b(ch) = median(snr(isfinite(snr)));
    end

    if any(isfinite(lar))
        larghezza_mediana_4b(ch) = median(lar(isfinite(lar)));
    end

    if ~isempty(isi)
        isi_mediano_4b(ch) = median(isi);
    end

    larghezze_valide_4b(ch) = sum(isfinite(lar));

end

tabella_spike_4b = table( ...
    canali_scelti(:), ...
    numero_spike(canali_scelti), ...
    ampiezza_mediana_4b(canali_scelti), ...
    snr_mediano_4b(canali_scelti), ...
    larghezza_mediana_4b(canali_scelti), ...
    larghezze_valide_4b(canali_scelti), ...
    isi_mediano_4b(canali_scelti), ...
    'VariableNames',{ ...
    'Canale', ...
    'NumeroEventi', ...
    'AmpiezzaMediana_uV', ...
    'SNRmediano', ...
    'LarghezzaMediana_ms', ...
    'LarghezzeValide', ...
    'ISImediano_ms'});

fprintf('\n===== CARATTERIZZAZIONE DEI CANDIDATI =====\n');

disp(tabella_spike_4b);


%% ========================================================================
% STEP 4B - DISTRIBUZIONE DELLE LARGHEZZE
% ========================================================================

figure('Color','w', ...
    'Name','Distribuzione delle larghezze dei candidati');

tiledlayout(2,2);

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    nexttile;
    hold on;

    lar = larghezze_4b{ch};
    lar = lar(isfinite(lar));

    if ~isempty(lar)

        histogram(lar, ...
            'FaceColor',[0.30 0.50 0.75]);

        xline(median(lar),'r', ...
            'LineWidth',2);

    end

    xlabel('Larghezza a meta'' ampiezza (ms)');
    ylabel('Numero di candidati');

    title(sprintf('Canale %d',ch));

    grid on;

    set(gca,'Color','w', ...
        'XColor','k','YColor','k');

    hold off;

end

sgtitle('Distribuzione delle larghezze a meta'' ampiezza');


%% ========================================================================
% STEP 4B - FREQUENZA DEI CANDIDATI NEL TEMPO
% ========================================================================
% Conteggio degli eventi in finestre consecutive da 10 s.

ampiezza_bin = 10;

bordi = 0:ampiezza_bin:durata_secondi;

if bordi(end) < durata_secondi
    bordi(end+1) = durata_secondi;
end

figure('Color','w', ...
    'Name','Frequenza temporale dei candidati');

tiledlayout(2,2);

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    conteggi = histcounts(tempi_spike{ch},bordi);

    durate_bin = diff(bordi);

    frequenze = conteggi ./ durate_bin;

    nexttile;

    stairs(bordi(1:end-1),frequenze, ...
        'LineWidth',1.5, ...
        'Color',[0.15 0.45 0.80]);

    xlabel('Tempo (s)');
    ylabel('Candidati/s');

    title(sprintf('Canale %d',ch));

    xlim([0 durata_secondi]);

    grid on;

    set(gca,'Color','w', ...
        'XColor','k','YColor','k');

end

sgtitle('Frequenza dei candidati in finestre di 10 s');


%% ========================================================================
% STEP 4C - CONTROLLO DI QUALITA' MORFOLOGICO
% ========================================================================
% Criteri esplorativi:
%
% Larghezza: 0.10-1.00 ms
% Correlazione con il template mediano: >= 0.80
%
% I template vengono calcolati separatamente per
% candidati positivi e negativi.
%
% La correlazione e' descrittiva: il template deriva
% dagli stessi eventi che vengono valutati.

s4c_corr = cell(n_channels,1);

s4c_n_valutabili = zeros(n_channels,1);
s4c_n_ok = zeros(n_channels,1);
s4c_n_da_rivedere = zeros(n_channels,1);

s4c_corr_med = nan(n_channels,1);
s4c_isi_brevi = zeros(n_channels,1);

s4c_ok_larg = cell(n_channels,1);
s4c_ok_corr = cell(n_channels,1);
s4c_ok_combinato = cell(n_channels,1);

for ch = 1:n_channels

    x = double(dati_filtrati_se(:,ch));
    x = x(:);

    locs = indici_spike{ch};
    pol = polarita_spike{ch};

    n = numel(locs);

    corr_eventi = nan(n,1);

    % Valutazione della larghezza
    lar = larghezze_4b{ch};

    ok_larg = isfinite(lar) & ...
        lar >= 0.10 & lar <= 1.00;

    ok_corr = false(n,1);

    % Si analizzano separatamente le due polarita'
    for segno = [-1 1]

        indici_pol = find(pol == segno);

        indici_validi = indici_pol( ...
            locs(indici_pol) > pre & ...
            locs(indici_pol)+post <= n_samples);

        % Servono almeno 3 eventi per costruire il template
        if numel(indici_validi) < 3
            continue;
        end

        onde = nan(pre+post+1,numel(indici_validi));

        for j = 1:numel(indici_validi)

            k = indici_validi(j);
            idx = locs(k);

            onde(:,j) = x(idx-pre:idx+post);

        end

        template = median(onde,2);

        for j = 1:numel(indici_validi)

            k = indici_validi(j);

            C = corrcoef(onde(:,j),template);

            if numel(C) == 4
                corr_eventi(k) = C(1,2);
            end

        end

    end

    ok_corr = isfinite(corr_eventi) & ...
        corr_eventi >= 0.80;

    ok_combinato = ok_larg & ok_corr;

    valutabili = isfinite(corr_eventi);

    s4c_corr{ch} = corr_eventi;

    s4c_ok_larg{ch} = ok_larg;
    s4c_ok_corr{ch} = ok_corr;
    s4c_ok_combinato{ch} = ok_combinato;

    s4c_n_valutabili(ch) = sum(valutabili);

    s4c_n_ok(ch) = sum(ok_combinato);

    s4c_n_da_rivedere(ch) = ...
        sum(valutabili & ~ok_combinato);

    if any(valutabili)

        s4c_corr_med(ch) = ...
            median(corr_eventi(valutabili));

    end

    s4c_isi_brevi(ch) = ...
        sum(isi_4b{ch} < 1.5);

end

s4c_tabella = table( ...
    canali_scelti(:), ...
    numero_spike(canali_scelti), ...
    s4c_n_valutabili(canali_scelti), ...
    s4c_n_ok(canali_scelti), ...
    s4c_n_da_rivedere(canali_scelti), ...
    s4c_corr_med(canali_scelti), ...
    s4c_isi_brevi(canali_scelti), ...
    'VariableNames',{ ...
    'Canale', ...
    'CandidatiTotali', ...
    'Valutabili', ...
    'CompatibiliCriteri', ...
    'DaRivedere', ...
    'CorrelazioneMediana', ...
    'ISIbrevi'});

fprintf('\n===== CONTROLLO DI QUALITA'' DEI CANDIDATI =====\n');

disp(s4c_tabella);


%% ========================================================================
% STEP 4C - DISTRIBUZIONE DELLE CORRELAZIONI
% ========================================================================

figure('Color','w', ...
    'Name','Somiglianza morfologica dei candidati');

tiledlayout(2,2);

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    nexttile;
    hold on;

    valori = s4c_corr{ch};
    valori = valori(isfinite(valori));

    if ~isempty(valori)

        histogram(valori, ...
            'BinEdges',-1:0.05:1, ...
            'FaceColor',[0.30 0.50 0.75]);

        xline(0.80,'r--','LineWidth',2);

    else

        text(0,0.5,'Dati insufficienti', ...
            'HorizontalAlignment','center');

    end

    xlabel('Correlazione con il template mediano');
    ylabel('Numero di candidati');

    title(sprintf('Canale %d',ch));

    xlim([-1 1]);

    grid on;

    set(gca,'Color','w', ...
        'XColor','k','YColor','k');

    hold off;

end

sgtitle('Somiglianza morfologica dei candidati spike');


%% ========================================================================
% STEP 4C - STABILITA' TEMPORALE DELLA MORFOLOGIA
% ========================================================================
% La registrazione viene suddivisa in tre periodi.
% Per confrontare correttamente le forme positive e negative,
% le onde vengono orientate rispetto al picco dominante.
%
% Questa rappresentazione e' utile per valutare la stabilita'
% della forma, ma non conserva l'informazione sulla polarita'.

bordi_periodi = linspace(0,durata_secondi,4);

colori_periodi = [
    0.15 0.45 0.85
    0.95 0.45 0.10
    0.85 0.70 0.10
];

figure('Color','w', ...
    'Name','Evoluzione temporale della morfologia');

tiledlayout(2,2);

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    x = double(dati_filtrati_se(:,ch));
    x = x(:);

    locs = indici_spike{ch};
    tempi = tempi_spike{ch};

    nexttile;
    hold on;

    for periodo = 1:3

        if periodo < 3

            selezione = ...
                tempi >= bordi_periodi(periodo) & ...
                tempi < bordi_periodi(periodo+1);

        else

            selezione = ...
                tempi >= bordi_periodi(periodo) & ...
                tempi <= bordi_periodi(periodo+1);

        end

        locs_periodo = locs(selezione);

        locs_periodo = locs_periodo( ...
            locs_periodo > pre & ...
            locs_periodo+post <= n_samples);

        if isempty(locs_periodo)
            continue;
        end

        onde = nan(pre+post+1,numel(locs_periodo));

        for k = 1:numel(locs_periodo)

            idx = locs_periodo(k);

            % Orientamento comune per il confronto
            onde(:,k) = ...
                sign(x(idx))*x(idx-pre:idx+post);

        end

        mediana_periodo = median(onde,2);

        plot(tempo_forma,mediana_periodo, ...
            'Color',colori_periodi(periodo,:), ...
            'LineWidth',1.8, ...
            'DisplayName',sprintf( ...
            'Periodo %d (n=%d)', ...
            periodo,numel(locs_periodo)));

    end

    xlabel('Tempo rispetto al picco (ms)');
    ylabel('Ampiezza orientata (\muV)');

    title(sprintf('Canale %d',ch));

    legend('show','Location','best');

    grid on;

    set(gca,'Color','w', ...
        'XColor','k','YColor','k');

    hold off;

end

sgtitle('Evoluzione temporale della morfologia dei candidati');


%% ========================================================================
% STEP 4C - STABILITA' TEMPORALE SEPARATA PER POLARITA'
% ========================================================================
% Analisi esplorativa della morfologia dei candidati positivi e negativi
% in tre periodi consecutivi della registrazione.
%
% Le forme d'onda mantengono il segno originale.
% I template mediani sono calcolati separatamente per polarita'.
%
% Se in un periodo sono presenti meno di 3 candidati, la mediana
% non viene rappresentata, per evitare confronti poco informativi.
%
% Questa analisi non modifica i candidati individuati nello STEP 4.

s4c_bordi_periodi = linspace(0,durata_secondi,4);
s4c_min_eventi_periodo = 3;

s4c_colori_periodi = [
    0.15 0.45 0.85
    0.95 0.45 0.10
    0.30 0.65 0.35
];

% Memorizzazione dei conteggi per ciascun canale e polarita'
s4c_conteggi_temporali = zeros(numel(canali_scelti),2,3);

figure('Color','w', ...
    'Name','Stabilita temporale separata per polarita');

tiledlayout(numel(canali_scelti),2, ...
    'TileSpacing','compact','Padding','compact');

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    x = double(dati_filtrati_se(:,ch));
    x = x(:);

    locs = indici_spike{ch};
    tempi = tempi_spike{ch};
    polarita = polarita_spike{ch};

    for ip = 1:2

        if ip == 1
            segno = -1;
            nome_polarita = 'Negativa';
        else
            segno = 1;
            nome_polarita = 'Positiva';
        end

        nexttile;
        hold on;

        for periodo = 1:3

            t_inizio = s4c_bordi_periodi(periodo);
            t_fine = s4c_bordi_periodi(periodo+1);

            if periodo < 3
                nel_periodo = ...
                    tempi >= t_inizio & tempi < t_fine;
            else
                nel_periodo = ...
                    tempi >= t_inizio & tempi <= t_fine;
            end

            selezione = ...
                nel_periodo & polarita == segno;

            locs_periodo = locs(selezione);

            % Esclusione delle forme incomplete ai bordi
            locs_periodo = locs_periodo( ...
                locs_periodo > pre & ...
                locs_periodo + post <= n_samples);

            n_eventi = numel(locs_periodo);

            s4c_conteggi_temporali(i,ip,periodo) = n_eventi;

            if n_eventi < s4c_min_eventi_periodo
                continue;
            end

            onde = nan(pre+post+1,n_eventi);

            for k = 1:n_eventi

                idx = locs_periodo(k);

                % Si mantiene la polarita' originale
                onde(:,k) = x(idx-pre:idx+post);

            end

            forma_mediana = median(onde,2);

            plot(tempo_forma,forma_mediana, ...
                'Color',s4c_colori_periodi(periodo,:), ...
                'LineWidth',1.8, ...
                'DisplayName',sprintf( ...
                'Periodo %d (n=%d)',periodo,n_eventi));

        end

        xline(0,'k:','HandleVisibility','off');
        yline(0,'k:','HandleVisibility','off');

        xlabel('Tempo rispetto al picco (ms)');
        ylabel('Ampiezza (\muV)');

        title(sprintf('Canale %d - Polarita %s', ...
            ch,nome_polarita));

        grid on;

        set(gca, ...
            'Color','w', ...
            'XColor','k', ...
            'YColor','k', ...
            'GridColor',[0.8 0.8 0.8], ...
            'GridAlpha',0.25);

        if ~isempty(findobj(gca,'Type','line', ...
                '-not','LineStyle',':'))
            legend('show','Location','best');
        end

        hold off;

    end

end

sgtitle('Stabilita temporale delle forme d''onda per polarita');


%% ========================================================================
% STEP 4C - TABELLA DEI CANDIDATI PER POLARITA' E PERIODO
% ========================================================================

fprintf('\n===== STABILITA'' TEMPORALE PER POLARITA'' =====\n');

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    negativi = reshape( ...
        s4c_conteggi_temporali(i,1,:),1,3);

    positivi = reshape( ...
        s4c_conteggi_temporali(i,2,:),1,3);

    tabella_periodi = table( ...
        (1:3)', ...
        negativi(:), ...
        positivi(:), ...
        'VariableNames',{ ...
        'Periodo', ...
        'CandidatiNegativi', ...
        'CandidatiPositivi'});

    fprintf('\nCanale %d:\n',ch);

    disp(tabella_periodi);

end



%% ========================================================================
% STEP 4D - CARATTERIZZAZIONE DEL FONDO AD ALTA FREQUENZA
% ========================================================================
% Si confrontano i candidati bipolari con estremi locali
% sub-soglia di entrambe le polarita'.
%
% Sono escluse le regioni entro +/- 5 ms dai candidati.
% Per ogni canale si considerano al massimo 500 estremi
% sub-soglia.
%
% Per il confronto morfologico tutte le onde vengono
% orientate in modo che il picco dominante sia positivo.

s4d_mad = nan(numel(canali_scelti),1);
s4d_rms_fondo = nan(numel(canali_scelti),1);

s4d_n_candidati = zeros(numel(canali_scelti),1);
s4d_n_minimi = zeros(numel(canali_scelti),1);

s4d_amp_med_candidati = nan(numel(canali_scelti),1);
s4d_amp_med_rumore = nan(numel(canali_scelti),1);

s4d_wave_candidati = cell(numel(canali_scelti),1);
s4d_wave_rumore = cell(numel(canali_scelti),1);

s4d_amp_candidati = cell(numel(canali_scelti),1);
s4d_amp_rumore = cell(numel(canali_scelti),1);

s4d_psd_completa = cell(numel(canali_scelti),1);
s4d_psd_fondo = cell(numel(canali_scelti),1);
s4d_freq = cell(numel(canali_scelti),1);

s4d_n_finestre = zeros(numel(canali_scelti),1);

esclusione = round(0.005*fs);

lunghezza_psd = round(0.100*fs);

for i = 1:numel(canali_scelti)

    ch = canali_scelti(i);

    x = double(dati_filtrati_se(:,ch));
    x = x(:);

    locs = indici_spike{ch};

    sigma = rumore_spike(ch);

    s4d_mad(i) = sigma;

    % ---------------- ESCLUSIONE DEI CANDIDATI ---------------------------

    maschera_fondo = true(n_samples,1);

    for k = 1:numel(locs)

        a = max(1,locs(k)-esclusione);
        b = min(n_samples,locs(k)+esclusione);

        maschera_fondo(a:b) = false;

    end

    fondo = x(maschera_fondo);

    if ~isempty(fondo)
        s4d_rms_fondo(i) = sqrt(mean(fondo.^2));
    end

    % ---------------- FORME DEI CANDIDATI -------------------------------

    locs_validi = locs( ...
        locs > pre & locs+post <= n_samples);

    onde_candidati = [];

    for k = 1:numel(locs_validi)

        idx = locs_validi(k);

        onda = sign(x(idx))*x(idx-pre:idx+post);

        onde_candidati(:,end+1) = onda;

    end

    s4d_wave_candidati{i} = onde_candidati;

    s4d_amp_candidati{i} = abs(x(locs_validi));

    s4d_n_candidati(i) = numel(locs);

    if ~isempty(locs_validi)

        s4d_amp_med_candidati(i) = ...
            median(abs(x(locs_validi)));

    end

    % ---------------- ESTREMI SUB-SOGLIA -------------------------------

    [~,locs_estremi] = findpeaks(abs(x), ...
        'MinPeakDistance',distanza_minima);

    locs_estremi = locs_estremi(:);

    validi = ...
        locs_estremi > pre & ...
        locs_estremi+post <= n_samples & ...
        abs(x(locs_estremi)) < fattore_soglia*sigma & ...
        abs(x(locs_estremi)) > 0;

    locs_estremi = locs_estremi(validi);

    % Si mantengono solo le finestre completamente esterne
    % alle regioni escluse attorno ai candidati.
    validi_fondo = false(size(locs_estremi));

    for k = 1:numel(locs_estremi)

        idx = locs_estremi(k);

        validi_fondo(k) = ...
            all(maschera_fondo(idx-pre:idx+post));

    end

    locs_estremi = locs_estremi(validi_fondo);

    % Campionamento uniforme lungo la registrazione
    if numel(locs_estremi) > 500

        posizioni = round(linspace( ...
            1,numel(locs_estremi),500));

        locs_estremi = locs_estremi(posizioni);

    end

    onde_fondo = [];

    for k = 1:numel(locs_estremi)

        idx = locs_estremi(k);

        onda = sign(x(idx))*x(idx-pre:idx+post);

        onde_fondo(:,end+1) = onda;

    end

    s4d_wave_rumore{i} = onde_fondo;

    s4d_amp_rumore{i} = abs(x(locs_estremi));

    s4d_n_minimi(i) = numel(locs_estremi);

    if ~isempty(locs_estremi)

        s4d_amp_med_rumore(i) = ...
            median(abs(x(locs_estremi)));

    end

    % ---------------- PSD DEL SEGNALE COMPLETO -------------------------

    finestra = hann(lunghezza_psd);

    [Pcompleta,f] = pwelch( ...
        x, ...
        finestra, ...
        round(lunghezza_psd/2), ...
        lunghezza_psd, ...
        fs);

    s4d_psd_completa{i} = Pcompleta;
    s4d_freq{i} = f;

    % ---------------- PSD DELLE FINESTRE DI FONDO ----------------------

    n_blocchi = floor(n_samples/lunghezza_psd);

    psd_blocchi = [];

    for k = 1:n_blocchi

        a = (k-1)*lunghezza_psd+1;
        b = k*lunghezza_psd;

        if all(maschera_fondo(a:b))

            segmento = x(a:b);

            [Pseg,~] = pwelch( ...
                segmento, ...
                finestra, ...
                0, ...
                lunghezza_psd, ...
                fs);

            psd_blocchi(:,end+1) = Pseg;

        end

    end

    s4d_n_finestre(i) = size(psd_blocchi,2);

    if ~isempty(psd_blocchi)

        s4d_psd_fondo{i} = mean(psd_blocchi,2);

    else

        s4d_psd_fondo{i} = nan(size(Pcompleta));

    end

end

% -------------------- TABELLA DEL FONDO ---------------------------------

s4d_tabella = table( ...
    canali_scelti(:), ...
    s4d_n_candidati, ...
    s4d_n_minimi, ...
    s4d_mad, ...
    s4d_rms_fondo, ...
    s4d_n_finestre, ...
    s4d_amp_med_candidati, ...
    s4d_amp_med_rumore, ...
    'VariableNames',{ ...
    'Canale', ...
    'Candidati', ...
    'EstremiSubSoglia', ...
    'MAD_uV', ...
    'RMSfondo_uV', ...
    'FinestrePSD', ...
    'AmpMedCandidati_uV', ...
    'AmpMedSubSoglia_uV'});

fprintf('\n===== FONDO EXTRACELLULARE AD ALTA FREQUENZA =====\n');

disp(s4d_tabella);


%% ========================================================================
% STEP 4D - FORME D'ONDA MEDIANE: CANDIDATI E SUB-SOGLIA
% ========================================================================

figure('Color','w', ...
    'Name','Candidati e fluttuazioni sub-soglia');

tiledlayout(2,2);

for i = 1:numel(canali_scelti)

    nexttile;
    hold on;

    onde_candidati = s4d_wave_candidati{i};
    onde_fondo = s4d_wave_rumore{i};

    if ~isempty(onde_candidati)

        plot(tempo_forma,median(onde_candidati,2), ...
            'r','LineWidth',2, ...
            'DisplayName','Candidati');

    end

    if ~isempty(onde_fondo)

        plot(tempo_forma,median(onde_fondo,2), ...
            'b','LineWidth',2, ...
            'DisplayName','Estremi sub-soglia');

    end

    xline(0,'k:');

    xlabel('Tempo rispetto al picco (ms)');
    ylabel('Ampiezza orientata (\muV)');

    title(sprintf('Canale %d',canali_scelti(i)));

    legend('show','Location','best');

    grid on;

    set(gca,'Color','w', ...
        'XColor','k','YColor','k');

    hold off;

end

sgtitle('Forme d''onda mediane: candidati ed estremi sub-soglia');


%% ========================================================================
% STEP 4D - PSD DEL SEGNALE COMPLETO E DEL FONDO
% ========================================================================

figure('Color','w', ...
    'Name','PSD alta frequenza: segnale completo e fondo');

tiledlayout(2,2);

for i = 1:numel(canali_scelti)

    nexttile;
    hold on;

    f = s4d_freq{i};

    Pcompleta = s4d_psd_completa{i};
    Pfondo = s4d_psd_fondo{i};

    selezione = f >= 300 & f <= 3000;

    plot(f(selezione),Pcompleta(selezione), ...
        'r','LineWidth',1.2, ...
        'DisplayName','Segnale completo');

    plot(f(selezione),Pfondo(selezione), ...
        'b','LineWidth',1.2, ...
        'DisplayName','Finestre di fondo');

    xlabel('Frequenza (Hz)');
    ylabel('PSD (\muV^2/Hz)');

    title(sprintf('Canale %d',canali_scelti(i)));

    legend('show','Location','best');

    grid on;

    set(gca,'Color','w', ...
        'XColor','k','YColor','k');

    hold off;

end

sgtitle('PSD 300-3000 Hz: segnale completo e finestre di fondo');


%% ========================================================================
% FINE DELL'ANALISI BIPOLARE DEI CANDIDATI
% ========================================================================

fprintf('\nAnalisi bipolare dei candidati completata.\n');
fprintf('Sono stati considerati picchi positivi e negativi.\n');
fprintf('I risultati devono essere verificati prima dell''interpretazione.\n');


%{

%% ============================================================
% STEP 6 - ANALISI SPETTRALE
% PSD, FONDO APERIODICO E POTENZA OSCILLATORIA
%
% Inserire alla fine di analysis(5).m, dopo lo STEP 5.
%
% Richiede:
%   dati_uV, fs, n_samples, n_channels
%   funzione PSD.m gia' presente nel progetto
%
% Non modifica lo spike detection.
% Non salva automaticamente le figure.
%% ============================================================


%% 6.1 - PARAMETRI

canali_psd_scelti = [24 26 27 38];

% Intervallo del fit aperiodico
f_fit_min = 0.5;
f_fit_max = 20;

% Banda per la potenza oscillatoria
banda_osc = [2.5 4.5];

% Frequenza massima mostrata nel primo grafico
f_plot_max = 100;

assert(size(dati_uV,1)==n_samples && ...
       size(dati_uV,2)==n_channels, ...
       'Dimensioni di dati_uV non coerenti.');

assert(fs > 2*f_plot_max, ...
       'Frequenza di campionamento insufficiente.');


%% 6.2 - CALCOLO DELLA PSD SU TUTTI I CANALI

% La funzione PSD.m utilizza pwelch con finestra di Hann.
%
% Durata finestra: 2 secondi
% Sovrapposizione: 50%, cioe' 1 secondo
% Risoluzione in frequenza: 0.5 Hz
%
% Segnale: dati_uV, NON filtrato a 300-3000 Hz.

[psd_step6, frequenze_step6] = PSD(dati_uV,fs,50,2);

psd_step6 = double(psd_step6);
frequenze_step6 = double(frequenze_step6(:));

% Evita log10(0)
psd_step6 = max(psd_step6,realmin('double'));

% Frequenze utilizzate per il fit
idx_fit = frequenze_step6 >= f_fit_min & ...
          frequenze_step6 <= f_fit_max;

% Frequenze della banda oscillatoria
idx_osc = frequenze_step6 >= banda_osc(1) & ...
          frequenze_step6 <= banda_osc(2);

% Frequenze da visualizzare
idx_plot = frequenze_step6 > 0 & ...
           frequenze_step6 <= f_plot_max;

assert(nnz(idx_fit)>=3 && nnz(idx_osc)>=2, ...
       'Risoluzione spettrale insufficiente.');

% Asse delle frequenze in scala logaritmica
log_f_fit = log10(frequenze_step6(idx_fit));

log_f_all = nan(size(frequenze_step6));

idx_positivi = frequenze_step6 > 0;

log_f_all(idx_positivi) = ...
    log10(frequenze_step6(idx_positivi));


%% 6.3 - INIZIALIZZAZIONE DEI RISULTATI

pendenza_aperiodica = nan(n_channels,1);
intercetta_aperiodica = nan(n_channels,1);

potenza_osc_log10 = nan(n_channels,1);
potenza_osc_dB = nan(n_channels,1);

rapporto_osc_fondo = nan(n_channels,1);

potenza_banda_uV2 = nan(n_channels,1);

% Valori del fit per ogni frequenza e canale
fit_log_psd = nan(numel(frequenze_step6),n_channels);


%% 6.4 - FIT LINEARE E POTENZA OSCILLATORIA

for ch = 1:n_channels

    % PSD in scala logaritmica
    log_psd = log10(psd_step6(:,ch));

    % PSD nell'intervallo 0.5-20 Hz
    y_fit = log_psd(idx_fit);

    % --------------------------------------------------------
    % FIT ROBUSTO IN SCALA LOG-LOG
    %
    % log10(PSD) = intercetta + pendenza*log10(f)
    % --------------------------------------------------------

    if exist('robustfit','file') == 2

        % Richiede Statistics and Machine Learning Toolbox
        coefficienti = robustfit( ...
            log_f_fit,y_fit,'bisquare');

    else

        % Alternativa se robustfit non e' disponibile
        warning(['robustfit non disponibile: ' ...
                 'utilizzo polyfit non robusto.']);

        p = polyfit(log_f_fit,y_fit,1);

        coefficienti = [p(2);p(1)];

    end

    intercetta_aperiodica(ch) = coefficienti(1);
    pendenza_aperiodica(ch) = coefficienti(2);

    % Ricostruzione del fondo aperiodico stimato
    fit_log_psd(idx_positivi,ch) = ...
        coefficienti(1) + ...
        coefficienti(2)*log_f_all(idx_positivi);

    % --------------------------------------------------------
    % POTENZA OSCILLATORIA RELATIVA AL FONDO
    % Banda: 2.5-4.5 Hz
    % --------------------------------------------------------

    residuo = log_psd(idx_osc) - ...
              fit_log_psd(idx_osc,ch);

    % Eccesso medio in unita' log10
    potenza_osc_log10(ch) = mean(residuo);

    % Eccesso medio in dB
    potenza_osc_dB(ch) = 10*mean(residuo);

    % Rapporto geometrico PSD/fondo
    rapporto_osc_fondo(ch) = ...
        10^potenza_osc_log10(ch);

    % --------------------------------------------------------
    % POTENZA TOTALE DELLA BANDA
    % Non corretta per il fondo aperiodico
    % --------------------------------------------------------

    potenza_banda_uV2(ch) = trapz( ...
        frequenze_step6(idx_osc), ...
        psd_step6(idx_osc,ch));

end


%% 6.5 - TABELLA GENERALE DEI 60 CANALI

tabella_psd = table( ...
    (1:n_channels)', ...
    pendenza_aperiodica, ...
    potenza_osc_log10, ...
    potenza_osc_dB, ...
    rapporto_osc_fondo, ...
    potenza_banda_uV2, ...
    'VariableNames', ...
    {'Canale', ...
     'PendenzaLogLog', ...
     'EccessoLog10', ...
     'Eccesso_dB', ...
     'RapportoOscFondo', ...
     'PotenzaTotaleBanda_uV2'});

fprintf('\n============================================\n');
fprintf('STEP 6 - ANALISI SPETTRALE\n');
fprintf('============================================\n');

fprintf('Canali analizzati: %d\n',n_channels);

fprintf('Welch: finestre 2 s, overlap 1 s\n');

fprintf('Risoluzione: %.2f Hz\n', ...
    frequenze_step6(2)-frequenze_step6(1));

fprintf('Fit aperiodico: %.1f-%.1f Hz\n', ...
    f_fit_min,f_fit_max);

fprintf('Banda oscillatoria: %.1f-%.1f Hz\n', ...
    banda_osc(1),banda_osc(2));

fprintf('\n===== TABELLA PSD - TUTTI I CANALI =====\n');

disp(tabella_psd);


%% 6.6 - RIEPILOGO DEI CANALI SELEZIONATI

fprintf('\n===== CANALI SELEZIONATI - PSD =====\n');

for ch = canali_psd_scelti

    fprintf('\nCanale %d\n',ch);

    fprintf('Pendenza aperiodica: %.3f\n', ...
        pendenza_aperiodica(ch));

    fprintf('Eccesso oscillatorio: %.4f log10\n', ...
        potenza_osc_log10(ch));

    fprintf('Eccesso oscillatorio: %.3f dB\n', ...
        potenza_osc_dB(ch));

    fprintf('Rapporto oscillazione/fondo: %.3f\n', ...
        rapporto_osc_fondo(ch));

    fprintf('Potenza totale %.1f-%.1f Hz: %.4g uV^2\n', ...
        banda_osc(1),banda_osc(2), ...
        potenza_banda_uV2(ch));

end


%% ============================================================
% 6.7 - GRAFICO 1
% PSD 0.5-100 Hz DEI QUATTRO CANALI
%% ============================================================

figure('Color','w', ...
       'Name','STEP 6 - PSD 0.5-100 Hz', ...
       'Position',[100 100 1200 700]);

for i = 1:numel(canali_psd_scelti)

    ch = canali_psd_scelti(i);

    subplot(2,2,i);

    semilogx( ...
        frequenze_step6(idx_plot), ...
        10*log10(psd_step6(idx_plot,ch)), ...
        'Color',[0.15 0.40 0.75], ...
        'LineWidth',1.2);

    xlabel('Frequenza (Hz)','Color','k');

    ylabel('PSD (dB re 1 \muV^2/Hz)','Color','k');

    title(sprintf('Canale %d',ch),'Color','k');

    xlim([0.5 f_plot_max]);

    grid on;
    box on;

    set(gca, ...
        'Color','w', ...
        'XColor','k', ...
        'YColor','k', ...
        'GridColor',[0.75 0.75 0.75]);

end

sgtitle('PSD del segnale originale - Metodo di Welch', ...
        'Color','k');


%% ============================================================
% 6.8 - GRAFICO 2
% PSD E FIT APERIODICO TRA 0.5 E 20 Hz
%% ============================================================

figure('Color','w', ...
       'Name','STEP 6 - PSD e fit aperiodico', ...
       'Position',[100 100 1200 700]);

for i = 1:numel(canali_psd_scelti)

    ch = canali_psd_scelti(i);

    subplot(2,2,i);
    hold on;

    % PSD sperimentale
    h_psd = plot( ...
        log10(frequenze_step6(idx_fit)), ...
        log10(psd_step6(idx_fit,ch)), ...
        'Color',[0.15 0.40 0.75], ...
        'LineWidth',1.3);

    % Fit del fondo aperiodico
    h_fit = plot( ...
        log10(frequenze_step6(idx_fit)), ...
        fit_log_psd(idx_fit,ch), ...
        '--', ...
        'Color',[0.85 0.30 0.15], ...
        'LineWidth',1.8);

    % Limiti della banda oscillatoria
    xline(log10(banda_osc(1)),':k');
    xline(log10(banda_osc(2)),':k');

    xlabel('log_{10}(Frequenza / Hz)','Color','k');

    ylabel('log_{10}(PSD / (\muV^2/Hz))', ...
           'Color','k');

    title(sprintf( ...
        'Canale %d - Eccesso %.2f dB', ...
        ch,potenza_osc_dB(ch)), ...
        'Color','k');

    legend([h_psd h_fit], ...
        {'PSD','Fit aperiodico'}, ...
        'Location','best');

    grid on;
    box on;

    set(gca, ...
        'Color','w', ...
        'XColor','k', ...
        'YColor','k', ...
        'GridColor',[0.75 0.75 0.75]);

    hold off;

end

sgtitle('PSD e fit lineare in scala log-log (0.5-20 Hz)', ...
        'Color','k');


%% ============================================================
% 6.9 - GRAFICO 3
% ECCESSO OSCILLATORIO DEI 60 CANALI
%% ============================================================

figure('Color','w', ...
       'Name','STEP 6 - Potenza oscillatoria 60 canali', ...
       'Position',[100 100 1200 550]);

bar(1:n_channels, ...
    potenza_osc_dB, ...
    'FaceColor',[0.30 0.50 0.75]);

hold on;

yline(0,'--k','Fondo aperiodico');

xlabel('Canale MATLAB','Color','k');

ylabel('Eccesso medio 2.5-4.5 Hz (dB)','Color','k');

title('Potenza oscillatoria relativa al fondo - 60 canali', ...
      'Color','k');

xlim([0.5 n_channels+0.5]);

xticks(1:2:n_channels);

grid on;
box on;

set(gca, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k', ...
    'GridColor',[0.75 0.75 0.75]);

hold off;


%% ============================================================
% STEP 6.10 - ANALISI DI CONTROLLO DELLE BASSE FREQUENZE
%
% Obiettivi:
% 1) Confrontare le PSD normalizzate
% 2) Calcolare la coerenza spettrale tra i canali
% 3) Visualizzare il segnale originale e la banda 1-4 Hz
%
% Canali selezionati: 24, 26, 27, 38
%
% Utilizza:
%   dati_uV
%   fs
%   n_samples
%   psd_step6
%   frequenze_step6
%
% Richiede Signal Processing Toolbox.
%
% Nessuna figura viene salvata automaticamente.
%% ============================================================


%% 6.10.1 - PARAMETRI GENERALI

canali_controllo = [24 26 27 38];

n_canali_controllo = numel(canali_controllo);

colori_controllo = lines(n_canali_controllo);

% Intervallo di frequenze da studiare
f_min_controllo = 0.5;
f_max_controllo = 20;

% Banda di interesse
banda_controllo = [2.5 4.5];

% Parametri Welch
durata_finestra_coh = 2;  % secondi

campioni_finestra_coh = ...
    round(durata_finestra_coh*fs);

campioni_overlap_coh = ...
    round(campioni_finestra_coh/2);

finestra_coh = hann(campioni_finestra_coh);

% Controlli
assert(all(canali_controllo <= size(dati_uV,2)), ...
    'Uno o piu'' canali non sono presenti nei dati.');

assert(size(psd_step6,2) >= max(canali_controllo), ...
    'La matrice PSD non contiene tutti i canali.');

fprintf('\n============================================\n');
fprintf('STEP 6.10 - CONTROLLI BASSE FREQUENZE\n');
fprintf('============================================\n');


%% ============================================================
% 6.10.2 - PSD NORMALIZZATE SOVRAPPOSTE
%
% Normalizzazione rispetto alla potenza media
% nell'intervallo 0.5-20 Hz.
%
% Questo grafico confronta la FORMA degli spettri,
% non la loro potenza assoluta.
%% ============================================================

idx_psd_controllo = ...
    frequenze_step6 >= f_min_controllo & ...
    frequenze_step6 <= f_max_controllo;

f_psd_controllo = ...
    frequenze_step6(idx_psd_controllo);

figure( ...
    'Color','w', ...
    'Name','STEP 6.10 - PSD normalizzate', ...
    'Position',[100 100 1100 600]);

hold on;

for c = 1:n_canali_controllo

    ch = canali_controllo(c);

    % PSD del canale
    P = double(psd_step6(idx_psd_controllo,ch));

    % Normalizzazione
    P_media = mean(P);

    P_normalizzata = P/P_media;

    % Conversione in dB
    P_normalizzata_dB = ...
        10*log10(P_normalizzata);

    % Grafico
    semilogx( ...
        f_psd_controllo, ...
        P_normalizzata_dB, ...
        'Color',colori_controllo(c,:), ...
        'LineWidth',1.8, ...
        'DisplayName',sprintf('Canale %d',ch));

end

% Limiti della banda di interesse
xline(banda_controllo(1),'--k', ...
      'HandleVisibility','off');

xline(banda_controllo(2),'--k', ...
      'HandleVisibility','off');

xlabel('Frequenza (Hz)','Color','k');

ylabel('PSD normalizzata (dB)','Color','k');

title('Confronto delle PSD normalizzate', ...
      'Color','k');

xlim([f_min_controllo f_max_controllo]);

legend('Location','best');

grid on;
box on;

set(gca, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k', ...
    'GridColor',[0.75 0.75 0.75]);

hold off;


%% ============================================================
% 6.10.3 - COERENZA SPETTRALE
%
% Calcolo della magnitude-squared coherence
% tra tutte le coppie dei quattro canali.
%
% La coerenza e' compresa tra 0 e 1.
%
% Un valore alto indica una relazione spettrale
% consistente, ma NON dimostra origine neuronale.
%% ============================================================

% Coppie di canali
coppie_canali = nchoosek(canali_controllo,2);

n_coppie = size(coppie_canali,1);

% Risultati
coerenza_2Hz = nan(n_coppie,1);

coerenza_media_banda = nan(n_coppie,1);

% Matrice contenente le curve di coerenza
coerenze_calcolate = [];

frequenze_coerenza = [];

figure( ...
    'Color','w', ...
    'Name','STEP 6.10 - Coerenza spettrale', ...
    'Position',[100 100 1200 850]);

for p = 1:n_coppie

    ch1 = coppie_canali(p,1);
    ch2 = coppie_canali(p,2);

    % Segnali originali
    x1 = double(dati_uV(:,ch1));
    x2 = double(dati_uV(:,ch2));

    x1 = x1(:);
    x2 = x2(:);

    % Rimozione offset
    x1 = x1-mean(x1);
    x2 = x2-mean(x2);

    % Coerenza spettrale
    [Cxy,f_coh] = mscohere( ...
        x1,x2, ...
        finestra_coh, ...
        campioni_overlap_coh, ...
        campioni_finestra_coh, ...
        fs);

    % Conservazione dei risultati
    if p == 1

        frequenze_coerenza = f_coh;

        coerenze_calcolate = ...
            nan(numel(f_coh),n_coppie);

    end

    coerenze_calcolate(:,p) = Cxy;

    % Coerenza al bin piu' vicino a 2 Hz
    [~,idx_2Hz] = min(abs(f_coh-2));

    coerenza_2Hz(p) = Cxy(idx_2Hz);

    % Coerenza media nella banda 2.5-4.5 Hz
    idx_banda = ...
        f_coh >= banda_controllo(1) & ...
        f_coh <= banda_controllo(2);

    coerenza_media_banda(p) = ...
        mean(Cxy(idx_banda));

    % Grafico della coppia
    subplot(3,2,p);

    idx_plot_coh = ...
        f_coh >= f_min_controllo & ...
        f_coh <= f_max_controllo;

    plot( ...
        f_coh(idx_plot_coh), ...
        Cxy(idx_plot_coh), ...
        'Color',[0.15 0.40 0.75], ...
        'LineWidth',1.5);

    hold on;

    xline(banda_controllo(1),'--k');
    xline(banda_controllo(2),'--k');

    xlabel('Frequenza (Hz)','Color','k');
    ylabel('Coerenza','Color','k');

    title(sprintf('Canali %d - %d',ch1,ch2), ...
          'Color','k');

    xlim([f_min_controllo f_max_controllo]);
    ylim([0 1]);

    grid on;
    box on;

    set(gca, ...
        'Color','w', ...
        'XColor','k', ...
        'YColor','k', ...
        'GridColor',[0.75 0.75 0.75]);

    hold off;

end

sgtitle('Coerenza spettrale tra i canali selezionati', ...
        'Color','k');


%% ============================================================
% 6.10.4 - TABELLA DELLA COERENZA
%% ============================================================

tabella_coerenza = table( ...
    coppie_canali(:,1), ...
    coppie_canali(:,2), ...
    coerenza_2Hz, ...
    coerenza_media_banda, ...
    'VariableNames', ...
    {'Canale1', ...
     'Canale2', ...
     'Coerenza_2Hz', ...
     'CoerenzaMedia_2p5_4p5Hz'});

fprintf('\n===== COERENZA TRA I CANALI =====\n');

disp(tabella_coerenza);


%% ============================================================
% 6.10.5 - SEGNALE ORIGINALE E FILTRATO 1-4 Hz
%
% Procedura:
% 1) Riduzione da fs a circa 200 Hz tramite resample
%    (con filtro antialiasing)
% 2) Butterworth 1-4 Hz con sezioni del secondo ordine
% 3) Filtraggio a fase zero
%
% Il filtro e' usato SOLO per visualizzare il segnale.
% Non modifica la PSD calcolata nello STEP 6.
%% ============================================================

% Finestra temporale da mostrare
t_inizio = 50;    % secondi
t_fine = 60;      % secondi

% Frequenza di campionamento desiderata
fs_ridotta_desiderata = 200;

% Fattori interi di ricampionamento
[p_resample,q_resample] = rat( ...
    fs_ridotta_desiderata/fs,1e-12);

% Frequenza effettiva dopo ricampionamento
fs_ridotta = fs*p_resample/q_resample;

% Filtro Butterworth ordine 4, passa-banda 1-4 Hz
% La rappresentazione SOS e' numericamente stabile.
[z_lf,p_lf,k_lf] = butter( ...
    4,[1 4]/(fs_ridotta/2),'bandpass');

[sos_lf,g_lf] = zp2sos(z_lf,p_lf,k_lf);

% Usiamo un oggetto filtro per il filtraggio a fase zero.
% In questo modo evitiamo di passare direttamente
% coefficienti b,a mal condizionati a 10 kHz.
filtro_lf = designfilt( ...
    'bandpassiir', ...
    'FilterOrder',8, ...
    'HalfPowerFrequency1',1, ...
    'HalfPowerFrequency2',4, ...
    'SampleRate',fs_ridotta);

figure( ...
    'Color','w', ...
    'Name','STEP 6.10 - Segnali 1-4 Hz', ...
    'Position',[100 100 1250 900]);

for c = 1:n_canali_controllo

    ch = canali_controllo(c);

    % Segnale originale in microvolt
    x = double(dati_uV(:,ch));
    x = x(:);

    % Rimozione dell'offset
    x = x-mean(x);

    % Ricampionamento con antialiasing
    x_ridotto = resample( ...
        x,p_resample,q_resample);

    % Tempo del segnale ricampionato
    t_ridotto = ...
        (0:numel(x_ridotto)-1)'/fs_ridotta;

    % Filtraggio passa-banda 1-4 Hz a fase zero
    x_filtrato = filtfilt(filtro_lf,x_ridotto);

    % Indici della finestra temporale
    idx_originale = ...
        (floor(t_inizio*fs)+1): ...
        min(floor(t_fine*fs),numel(x));

    idx_ridotto = ...
        t_ridotto >= t_inizio & ...
        t_ridotto <= t_fine;

    t_originale = (idx_originale-1)/fs;

    % --------------------------------------------------------
    % Grafico del segnale originale
    % --------------------------------------------------------

    subplot(4,2,2*c-1);

    plot( ...
        t_originale, ...
        x(idx_originale), ...
        'Color',[0.35 0.35 0.35], ...
        'LineWidth',0.7);

    xlabel('Tempo (s)','Color','k');
    ylabel('Ampiezza (\muV)','Color','k');

    title(sprintf( ...
        'Canale %d - Segnale originale',ch), ...
        'Color','k');

    xlim([t_inizio t_fine]);

    grid on;
    box on;

    set(gca, ...
        'Color','w', ...
        'XColor','k', ...
        'YColor','k', ...
        'GridColor',[0.75 0.75 0.75]);

    % --------------------------------------------------------
    % Grafico del segnale filtrato
    % --------------------------------------------------------

    subplot(4,2,2*c);

    plot( ...
        t_ridotto(idx_ridotto), ...
        x_filtrato(idx_ridotto), ...
        'Color',colori_controllo(c,:), ...
        'LineWidth',1.3);

    xlabel('Tempo (s)','Color','k');
    ylabel('Ampiezza (\muV)','Color','k');

    title(sprintf( ...
        'Canale %d - Filtrato 1-4 Hz',ch), ...
        'Color','k');

    xlim([t_inizio t_fine]);

    grid on;
    box on;

    set(gca, ...
        'Color','w', ...
        'XColor','k', ...
        'YColor','k', ...
        'GridColor',[0.75 0.75 0.75]);

end

sgtitle( ...
    'Segnale originale e componente a bassa frequenza', ...
    'Color','k');

%}
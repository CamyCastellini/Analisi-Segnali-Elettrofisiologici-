function p = plots()

p.tempo = @plots_tempo;
p.frequenze = @plots_frequenze;
p.psd = @plots_PSD;
p.tempi = @plots_stacked;


end

function plots_tempo(time, dati, timestamps, testoY, titolo)

figure;

% Tre canali, tre grafici uno sotto l'altro
for channel = 1:3

    subplot(3,1,channel);

    % Visualizzazione del canale selezionato
    plot(timestamps(time), dati(time,channel));

    ylabel(testoY);
    % Titolo del singolo canale
    title(sprintf('Canale %d', channel));
    grid on;

end

xlabel('Tempo (s)');
sgtitle(titolo);

end

function plots_frequenze(frequencies, s, testoY, titolo)

figure;

for channel = 1:3

    subplot(3,1,channel);

    plot(frequencies, s(:,channel));

    %PER VEDERE LE LINEE VERTICALI PIù CHIARAMENTE
    %stem(frequencies, s(:,channel), 'Marker', 'none');

    xlim([0 200]);

    ylabel(testoY);

    % Titolo del singolo canale
    title(sprintf('Canale %d', channel));

    grid on;

end

xlabel('Frequenza (Hz)');
sgtitle(titolo);

end

function plots_PSD(f, pxx, testoY, titolo)

figure;

for channel = 1:3
    subplot(3,1,channel);
    %unità di misura: uV^2/Hz--> densità di potenza
    %passo alla scala log. in decibel
    plot(f, 10*log10(pxx(:,channel)));

    %xlim([0 1000]);
    ylabel(testoY);
    title(sprintf('Canale %d',channel));
    grid on;
end

xlabel('Frequenza (Hz)');
sgtitle(titolo);

end


%VISUALIZZAZIONE DI PIù CANALI IN CONTEMPORANEA
function plots_stacked(time, dati, timestamps, canali)

% Selezione dell'intervallo temporale
t = timestamps(time);
segnali = double(dati(time, canali));

% Controllo dei dati
if isempty(t) || isempty(segnali)
    error('Intervallo temporale vuoto');
end

% Rimozione della media per ogni canale
segnali = segnali - mean(segnali, 1);

% Numero di canali selezionati
n = length(canali);

% Ampiezza massima assoluta di tutti i segnali
ampiezza = max(abs(segnali(:)));

% Distanza verticale tra le tracce
distanza = 2.5 * ampiezza;

% Evita distanza nulla se i segnali sono costanti
if distanza == 0
    distanza = 1;
end

% Posizioni verticali delle tracce
posizioni = (n-1:-1:0) * distanza;

% Creazione della figura
figure;
hold on;

% Colori differenti per ogni canale
colori = lines(n);

% Visualizzazione delle tracce
for i = 1:n

    offset = posizioni(i);

    plot(t, segnali(:,i) + offset, ...
        'Color', colori(i,:), ...
        'LineWidth', 0.8);

end

% Posizioni ordinate dal basso verso l'alto
yticks(fliplr(posizioni));
% Etichette nello stesso ordine
yticklabels(compose('Canale %d',fliplr(canali)));

% Etichette degli assi
xlabel('Tempo (s)');
ylabel('Canali');

% Titolo
title('Confronto temporale dei segnali');

% Limiti temporali
if numel(t) > 1
    xlim([t(1) t(end)]);
end

% Limiti verticali
ylim([-distanza/2, max(posizioni)+distanza/2]);

grid on;
box on;
hold off;

end
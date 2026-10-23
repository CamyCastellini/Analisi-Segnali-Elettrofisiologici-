function h = heatmaps()

h.tempo = @heatmap_tempo;
h.frequenze = @heatmap_frequenze;

end

function heatmap_tempo(timestamps, dati, titolo)

figure;
imagesc(timestamps, 1:size(dati,2), dati.');
axis xy;

xlabel('Tempo (s)');
ylabel('Canale');
title(titolo);

cb = colorbar;
cb.Label.String = 'Ampiezza (uV)';

colormap(turbo);
end

function heatmap_frequenze(f, valori, testoColore, titolo)

figure;
imagesc(f, 1:size(valori,2), valori.');
axis xy;

xlabel('Frequenza (Hz)');
ylabel('Canale');
title(titolo);

cb = colorbar;
cb.Label.String = testoColore;

colormap(turbo);
end
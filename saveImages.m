figures = findall(groot, 'Type', 'figure');

% Ordina le figure per numero
[~, indici] = sort([figures.Number]);
figures = figures(indici);

% Salva nella cartella corrente
nomeFile = fullfile(pwd, 'RasterPlots.pdf');

for i = 1:numel(figures)
    exportgraphics(figures(i), nomeFile, ...
        'ContentType', 'image', ...
        'Resolution',300, ...
        'Append', i > 1);
    fprintf('Salvata figura %d di %d\n', i, numel(figures));
end

function s1 = FFT(dati, n_samples)

% Finestra di Hann
w = hann(n_samples);

% Conversione in double e applicazione della finestra
dati = double(dati);
dati_windowed = dati .* w;

% FFT lungo le righe (una FFT per canale)
Y = fft(dati_windowed, [], 1);

% Spettro bilatero normalizzato
s2 = abs(Y / sum(w));

% Spettro monolatero
s1 = s2(1:floor(n_samples/2)+1, :);

% Raddoppio delle ampiezze
if mod(n_samples, 2) == 0
    s1(2:end-1, :) = 2*s1(2:end-1, :);
else
    s1(2:end, :) = 2*s1(2:end, :);
end

end

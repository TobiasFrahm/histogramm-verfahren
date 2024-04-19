clear all;
close all;
%%

freq = 100;
n_Periods = 10;
Fs = 1000*freq;             % sample frequency
T = 1/Fs;
dT = T;
n_Samples = round(1/freq*n_Periods*Fs);
time = (0:(n_Samples-1))*T;        % Time vector
f = Fs*(0:(n_Samples/2))/n_Samples;

maxVoltage = 1;
minVoltage = -1;

signal_sin1 = 0.9 * sin(2 * pi * freq * time);
signal_sin2 = 1.1 * sin(2 * pi * freq * time);
[signal_sin2_cc, signal_sin2_clean, signal_sin2] = clean_clipped(signal_sin2, minVoltage, maxVoltage);


signal_sat = 10 * sin(2 * pi * freq * time);
[signal_sat_cc, signal_sat_clean, signal_sat] = clean_clipped(signal_sat, minVoltage, maxVoltage);

n_bins = 100;

a = 0.5;
b = a;
x = linspace(0, 1, n_bins);
x_plt = linspace(-0.9, 0.9, n_bins);
bb = betapdf(x, a, b);

%%
n_fig = 0;
n_fig = n_fig + 1;
fig = figure(n_fig);
tiledlayout(3, 2, 'TileSpacing','compact','Padding','compact')

nexttile;
plot(time, signal_sin1)
grid on; 
xlabel("time [s]")
ylabel("voltage [V]")
title(["10k sample | 1Hz | 10 periods",  "no clippping | no noise"])

nexttile;
hold on;
h1 = histogram(signal_sin1, n_bins, 'Normalization', 'probability');
grid on;
xlabel("voltage [V]")
ylabel("hits")
title(["100 bins | 10k samples", "Varianz/Kurtosis: " + round(var(signal_sin1), 2) + " / " + round(kurtosis(signal_sin1), 2)])
plot(x_plt, bb*max(h1.Values)/6, "LineWidth", 2)
legend(["Measured Distrubution", "Beta Distribution"])


nexttile;
plot(time, signal_sin2)
grid on; 
xlabel("time [s]")
ylabel("voltage [V]")
title(["10k sample | 1Hz | 10 periods", "partly clipped | no noise"])

nexttile;
h2 = histogram(signal_sin2, n_bins,  'Normalization', 'probability');
grid on; 
xlabel("voltage [V]")
ylabel("hits")
title(["100 bins | 10k samples", "Varianz/Kurtosis: " + round(var(signal_sin2_clean), 2) + " / " + round(kurtosis(signal_sin2_clean), 2)])

nexttile;
plot(time, signal_sat)
grid on; 
xlabel("time [s]")
ylabel("voltage [V]")
title(["10k sample | 1Hz | 10 periods", "full clipped | no noise"])

nexttile;
h3 = histogram(signal_sat, n_bins,  'Normalization', 'probability');
grid on; 
xlabel("voltage [V]")
ylabel("hits")
title(["100 bins | 10k samples", "Varianz/Kurtosis: " + round(var(signal_sat_clean), 2) + " / " + round(kurtosis(signal_sat_clean), 2)])

% exportgraphics(fig, "img/beta-distribution.pdf")
%% Rauschen

n_fig = n_fig + 1;
fig = figure(n_fig);
tiledlayout(3, 2, 'TileSpacing','compact','Padding','compact')

signal = 0.8 * sin(2 * pi * freq * time);
snr = linspace(40, 10, 3);
for ss = 1:length(snr)
    sig_noise = awgn(signal, snr(ss));
    sig_noise(sig_noise > maxVoltage) = maxVoltage;
    sig_noise(sig_noise < minVoltage) = minVoltage;
    [clipp_count, signal_clean, signal_clipped] = clean_clipped(sig_noise, minVoltage, maxVoltage);
    nexttile;
    plot(time, sig_noise)
    grid on; 
    xlabel("time [s]")
    ylabel("voltage [V]")
    title("10k sample | 1Hz | 10 periods | SNR: " + snr(ss) + " dB")
    
    nexttile;
    h3 = histogram(sig_noise, n_bins,  'Normalization', 'probability');
    grid on; 
    xlabel("voltage [V]")
    ylabel("hits")
    title(["100 bins | 10k samples", "Varianz/Kurtosis: " + round(var(signal_clean), 2) + " / " + round(kurtosis(signal_clean), 2)])
    disp("Schiefe: " + skewness(signal_clipped) + " | " + skewness(sig_noise))
    
end

% exportgraphics(fig, "img/noise-histogramm.pdf")
%% Fehler in der Amplitude durch Sättigung

n_fig = n_fig + 1;
fig = figure(n_fig);
hold on;
signal = 0.9 * sin(2 * pi * freq * time);
freq_idx = 11;
fft_clean = fft(signal/n_Samples);
fft_clean = fft_clean(1:n_Samples/2+1);
fft_clean(2:end-1) = 2*fft_clean(2:end-1);
referenceampl = fft_clean(freq_idx);
snr2 = linspace(40, 5, 1000);

seed = 1234;
rng(seed);

for ss = 1:length(snr2)
    sig_noise = awgn(signal, snr2(ss));
    [cnt, ~, sig_noise] = clean_clipped(sig_noise, minVoltage, maxVoltage);
    fft_sig = fft(sig_noise/n_Samples);
    fft_sig = fft_sig(1:n_Samples/2+1);
    fft_sig(2:end-1) = 2*fft_sig(2:end-1);

    f = Fs*(0:(n_Samples/2))/n_Samples;
    
    noiseampl = fft_sig(freq_idx);
    err(ss) = abs(referenceampl/noiseampl);

   cnts(ss) = cnt;
end

for ss = 1:length(cnts)
    cnts(ss);
    if cnts(ss) > (n_Samples*0.001)
        snr_idx = ss;
        break;
    end
end

semilogy(snr2, err)
xline(snr2(snr_idx), 'r')

grid on; 
xlabel("SNR [dB]")
ylabel("Abweichung der Amplitude")
title("10k sample | 10 periods | \alpha = " + round(abs(referenceampl), 2))
legend("Abweichung")
% Create arrow
annotation(fig,'arrow',[0.732142857142857 0.641071428571428],...
    [0.447619047619048 0.390476190476191]);

% Create textbox
annotation(fig,'textbox',...
    [0.65278571428571 0.44761904761905 0.245428571428576 0.0500000000000004],...
    'String','0.1% Sättigungsgrad',...
    'FitBoxToText','off');


exportgraphics(fig, "img/noise-err.pdf")
% single current pulse to charge/discharge from a starting SoC
% OCV curve is modelled in a very simple manner
%%
clear; clc; close all;
%%
freq_list = logspace(0, 4, 50);
n_Freq = length(freq_list);

% single P45B
R0 = 6E-3; % Ohm
R1 = 4E-3;   % Ohm
C1 = 500E-3;    % capacity, F

%%
nn = 17;

%% noise
snr = [80 ,linspace(40, -5, nn-1)];
n_SNR = length(snr);
snr_fokus = 5;
noise_fit_cutoff = 2.0;

%% gain
gain_V = linspace(120, 180, nn);
n_gain_Voltage = length(gain_V);
gain_Current = 1;

%% quantization

nbitsVoltage = 12;   % 4096 qantization levels
qLevelsVoltage = 2^nbitsVoltage;
signalMinVoltage = 0;  % 0 ... 3 V
signalMaxVoltage = 3.3;
scalingFactorVoltage = (signalMaxVoltage-signalMinVoltage)/qLevelsVoltage;
partition_Voltage = signalMinVoltage:scalingFactorVoltage:signalMaxVoltage; % Length X, to represent X+1 intervals
codebook_Voltage = signalMinVoltage-scalingFactorVoltage:scalingFactorVoltage:signalMaxVoltage; % Length X+1, one entry for each interval

nbitsCurrent = 12;   % 4096 qantization levels
qLevelsCurrent = 2^nbitsCurrent;
signalMinCurrent = 0;  % 0 ... 3 V
signalMaxCurrent = 3.3;
scalingFactorCurrent = (signalMaxCurrent-signalMinCurrent)/qLevelsCurrent;
partition_Current = signalMinCurrent:scalingFactorCurrent:signalMaxCurrent; % Length X, to represent X+1 intervals
codebook_Current = signalMinCurrent-scalingFactorCurrent:scalingFactorCurrent:signalMaxCurrent; % Length X+1, one entry for each interval

%% histogram

n_Bins = 4095;
limit_min_percent = 5;
limit_max_percent = 95;
limit_min = 3.3*limit_min_percent/100;
limit_max = 3.3*limit_max_percent/100;

%% RRC - Modell - Transient
for ss=1:n_SNR
    for ff=1:n_Freq

        freq = freq_list(ff);       % EIS frequency
        n_Periods = 10;
        Fs = 1000*freq;             % sample frequency
        T = 1/Fs;
        dT = T;
        n_Samples = round(1/freq*n_Periods*Fs);
        time{ff} = (0:(n_Samples-1))*T;        % Time vector
        f = Fs*(0:(n_Samples/2))/n_Samples;

        signal_to_sim = 1 * sin(2*pi*freq*time{ff});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % to simulate the measurement: add noise to current input
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        current_clean = signal_to_sim;
        current = awgn(signal_to_sim, snr(ss));

        % status
        disp("   " + ff + " / " + n_Freq + " | f = " + freq + " Hz | " + "SNR: " + snr(ss) + " dB")

        voltage = Transient_RRC(Fs, current, R0, R1, C1);
        voltage_clean = Transient_RRC(Fs, current_clean, R0, R1, C1);
        % voltage_aged = Transient_RRC(Fs, current_clean, R00, R11, C11);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate impedance of the measured signal (no quantization, gain etc.)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f_freq_idx = n_Periods + 1; % select frequency line
        I_cmplx = fft(current);
        U_cmplx = fft(voltage);
        Z2{ss}(ff) = U_cmplx(f_freq_idx)./I_cmplx(f_freq_idx)*1E6;   %Ohm to uOhm

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate reference impedance (current & voltage without any
        % noise, quantization etc.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_cmplx = fft(current_clean);
        U_cmplx = fft(voltage_clean);
        Z2_clean{ss}(ff) = U_cmplx(f_freq_idx)./I_cmplx(f_freq_idx)*1E6;   %Ohm to uOhm

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate impedance of an "aged" cell, i.e. slightly different
        % values for R0, R1 and C1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % I_cmplx = fft(current);
        % U_cmplx = fft(voltage_aged);
        % Z2_aged{ss}(ff) = U_cmplx(f_freq_idx)./I_cmplx(f_freq_idx)*1E6;   %Ohm to uOhm

        if ff==n_Freq || ff==1
            voltage_plot{ss,ff} = voltage;
            current_plot{ss,ff} = current;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % to simulate the measurement: substract DC, scale signal (imitate gain), do a quantization, cut off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for gg=1:n_gain_Voltage
            gain_Voltage = gain_V(gg);

            % substract DC
            DC_Voltage = mean(voltage);
            DC_Current = mean(current);
            DC_Voltage_clean = mean(voltage_clean);
            DC_Current_clean = mean(current_clean);
            % DC_Voltage_aged = mean(voltage_aged);
            % DC_Current_aged = mean(current_clean);

            voltage_Measured = voltage - DC_Voltage;
            current_Measured = current - DC_Current;
            
            voltage_Reference = voltage_clean - DC_Voltage_clean;
            % voltage_Aged = voltage_aged - DC_Voltage_aged;

            voltage_Measured_clean = voltage_clean - DC_Voltage_clean;
            current_Measured_clean = current_clean - DC_Current_clean;


            % Amplify
            voltage_Reference = voltage_Reference * gain_Voltage + 1.65;
            % voltage_Aged = voltage_Aged * gain_Voltage * 1.65;
            
            voltage_Measured = voltage_Measured * gain_Voltage + 1.65; % 1.65 account for ADC measurement from 0 ... 3V3 Volt
            % voltage_Measured = awgn(voltage_Measured, snr(ss)); 
            current_Measured = current_Measured * gain_Current + 1.65;

            voltage_Measured_clean = voltage_Measured_clean * gain_Voltage + 1.65; % 1.65 account for ADC measurement from 0 ... 3V3 Volt
            current_Measured_clean = current_Measured_clean * gain_Current + 1.65;

            % Quantization voltage & current
            [~,voltage_Measured] = quantiz(voltage_Measured, partition_Voltage, codebook_Voltage);
            [~,current_Measured] = quantiz(current_Measured, partition_Current, codebook_Current);
            [~,voltage_Measured_clean] = quantiz(voltage_Measured_clean, partition_Voltage, codebook_Voltage);
            [~,current_Measured_clean] = quantiz(current_Measured_clean, partition_Current, codebook_Current);
  
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Clipping: cut off signal above 3.3V
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            histo_clip_count{ff}(ss, gg) = (length(voltage_Measured(voltage_Measured<signalMinVoltage)) + length(voltage_Measured(voltage_Measured>signalMaxVoltage)));

            [cutoff_values_Voltage_percent{ss,gg}(ff), voltage_Measured] = clip(voltage_Measured, signalMinVoltage, signalMaxVoltage);
            [~, current_Measured] = clip(current_Measured, signalMinVoltage, signalMaxVoltage);

            [cutoff_values_Voltage_percent_clean{ss,gg}(ff), voltage_Measured_clean] = clip(voltage_Measured_clean, signalMinVoltage, signalMaxVoltage);
            [~, current_Measured_clean] = clip(current_Measured_clean, signalMinVoltage, signalMaxVoltage);

            % [cutoff_values_Voltage_percent_aged{ss,gg}(ff), voltage_Aged] = clip(voltage_Aged, signalMinVoltage, signalMaxVoltage);
            % [~, current_Measured_aged] = clip(current_Measured_aged, signalMinVoltage, signalMaxVoltage);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do histogram and calculate values inside interval and outside
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            bin_edges = linspace(signalMinVoltage - signalMaxVoltage/n_Bins, signalMaxVoltage, n_Bins + 1);
            
            histo_Voltage_Counts{ss,gg,ff} = histcounts(voltage_Measured, bin_edges);
            axis_bins = bin_edges(2:end);
            histo_Voltage_outside{ss,gg}(ff) = ...
                ( sum(histo_Voltage_Counts{ss,gg,ff}(axis_bins<limit_min)) + sum(histo_Voltage_Counts{ss,gg,ff}(axis_bins>limit_max)) ) ...
                / sum(histo_Voltage_Counts{ss,gg,ff}) * 100;
            vv_m = voltage_Measured;
            vv_m(vv_m <= signalMinVoltage) = [];
            vv_m(vv_m >= signalMaxVoltage) = [];
            histo_Kurtosis{ff}(ss, gg) = kurtosis(vv_m, 0);
            histo_Varianz{ff}(ss, gg) = var(vv_m);
            

            U_cmplx = fft(voltage_Measured)/gain_Voltage;
            I_cmplx = fft(current_Measured)/gain_Current;
            Z{ss,gg}(ff) = U_cmplx(f_freq_idx)./I_cmplx(f_freq_idx)*1E6;   %Ohm to uOhm
            
            fft_u = fft(voltage_Measured);
            U_fft{ff}(gg, ss) = (fft_u(f_freq_idx));
            fft_ref = fft(voltage_Reference);
            U_fft_ref{ff}(gg, ss) = (fft_ref(f_freq_idx));

            U_cmplx_out{ss, gg}(ff, :) = U_cmplx;
            U_cmplx_abs_out(ss, gg, ff) = abs(U_cmplx(f_freq_idx));
            I_cmplx_out{ss, gg}(ff, :) = I_cmplx;

            voltage_Measured_plot{ss,gg,ff} = voltage_Measured;  % save for plot
            current_Measured_plot{ss,gg,ff} = current_Measured;  % save for plot

            close all

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate factor to correct the voltage amplitude
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sollamp{ff}(ss, gg) = abs(U_fft_ref{ff}(gg, ss))./abs(U_fft{ff}(gg, ss));   % ideal

        end
    end
end


%%
figure(1)
tiledlayout(1,2)
ff = n_Freq;

cc = 1;

nexttile(1)
plot(time{ff},current_plot{cc,ff},'DisplayName',"Current: f= " + num2str(freq_list(ff)) + " Hz")
ylabel('current (A)')
yyaxis right
hold on
plot(time{ff},voltage_plot{cc,ff},'DisplayName',"Voltage: f= " + num2str(freq_list(ff)) + " Hz")
ylabel('voltage (V)')
xlabel('time (s)')
legend
title('Actual current and voltage @ "Battery", i.e. float values')
grid on; box on

nexttile(2)
plot(time{ff},current_Measured_plot{cc,n_gain_Voltage,ff},'DisplayName',"Current: f= " + num2str(freq_list(ff)) + " Hz / SNR= " + num2str(n_gain_Voltage(1)) + " dB")
ylabel('current (A)')
yyaxis right
hold on
plot(time{ff},voltage_Measured_plot{cc,n_gain_Voltage,ff},'DisplayName',"Voltage: f= " + num2str(freq_list(ff)) + " Hz / SNR= " + num2str(gain_V(n_gain_Voltage)) + " dB")
ylabel('voltage (V)')
xlabel('time (s)')
legend
title('Current and voltage values @ "ADC", i.e. quantized and scaled values')
grid on; box on

%%
figure(2)
tiledlayout(2,1,"TileSpacing","tight","Padding","tight")

i_plot1 = 3;
i_plot2 = 4;%ceil(n_gain_Voltage/2);
i_plot3 = n_gain_Voltage;
cc = (n_SNR - snr_fokus);
nexttile
hold on
plot(time{1},voltage_Measured_plot{cc,i_plot1,1}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot1)) + " / f= " + num2str(freq_list(1)) + " Hz");
plot(time{1},voltage_Measured_plot{cc,i_plot2,1}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot2)) + " / f= " + num2str(freq_list(1)) + " Hz");
% plot(time{n_Freq},voltage_Measured_plot{cc,i_plot3,1}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot3)) + " / f= " + num2str(freq_list(1)) + " Hz");
ylim([0 3.3])
ylabel('U (V)')
legend
title('Voltage')
box on; grid on;
hold off


nexttile
hold on
plot(time{n_Freq},voltage_Measured_plot{cc,i_plot1,n_Freq}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot1)) + " / f= " + num2str(freq_list(n_Freq)) + " Hz");
plot(time{n_Freq},voltage_Measured_plot{cc,i_plot2,n_Freq}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot2)) + " / f= " + num2str(freq_list(n_Freq)) + " Hz");
% plot(time{n_Freq},voltage_Measured_plot{cc,i_plot3,n_Freq}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot3)) + " / f= " + num2str(freq_list(n_Freq)) + " Hz");
ylim([0 3.3])
ylabel('U (V)')
legend
title('Voltage')
box on; grid on;
hold off


%%
figure(3)
tiledlayout(3,1,"TileSpacing","tight","Padding","tight")

i_plot1 = 1;
i_plot2 = ceil(n_gain_Voltage/2);
i_plot3 = n_gain_Voltage;
ff = 1; %n_Freq

nexttile
hold on
plot(time{ff},voltage_Measured_plot{ss,i_plot1,ff}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot1)) + " / f= " + num2str(freq_list(ff)) + " Hz");
plot(time{ff},voltage_Measured_plot{ss,i_plot2,ff}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot2)) + " / f= " + num2str(freq_list(ff)) + " Hz");
plot(time{ff},voltage_Measured_plot{ss,i_plot3,ff}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot3)) + " / f= " + num2str(freq_list(ff)) + " Hz");
ylim([0 3.3])
ylabel('U (V)')
legend
title('Voltage')
box on; grid on;
hold off

nexttile
hold on
plot(axis_bins,histo_Voltage_Counts{ss,i_plot1,ff}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot1)) + " / f= " + num2str(freq_list(ff)) + " Hz");
plot(axis_bins,histo_Voltage_Counts{ss,i_plot2,ff}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot2)) + " / f= " + num2str(freq_list(ff)) + " Hz");
plot(axis_bins,histo_Voltage_Counts{ss,i_plot3,ff}, 'DisplayName',"voltage gain= " + num2str(gain_V(i_plot3)) + " / f= " + num2str(freq_list(ff)) + " Hz");
xline([limit_min limit_max],'-',{'lower boundary 5%','upper boundary 95%'},'LineWidth',2, 'HandleVisibility', 'off')
hold off
xlim([0 3.3])
xlabel('U (V)')
ylabel('Counts Histogram')
legend
box on; grid on;

nexttile
for gg=1:n_gain_Voltage
    outside(gg) = histo_Voltage_outside{ss,gg}(ff);
end
plot(gain_V, outside, 'o-', 'DisplayName', "Values outside for f = " + num2str(freq_list(ff)) + " Hz")
xlabel('Gain voltage')
ylabel('Values outside interval (%)')
legend
grid on; box on


%% Erzeugen der LUT
clear lut imlut ls lut_cnt

ls = 100;    % dimension lut = (ls)^3

% clipped samples
min_ccb = min(min(min(cell2mat(histo_clip_count))));
max_ccb = max(max(max(cell2mat(histo_clip_count))));
span_ccb = max_ccb - min_ccb;

% varianz
min_vv = min(min(min(cell2mat(histo_Varianz))));
max_vv = max(max(max(cell2mat(histo_Varianz))));
span_vv = max_vv - min_vv;

% kurtosis
min_kk = min(min(min(cell2mat(histo_Kurtosis))));
max_kk = max(max(max(cell2mat(histo_Kurtosis))));
span_kk = max_kk - min_kk;

% ohne  griddata
for ff = 1:n_Freq
    lut{ff} = zeros(ls,ls,ls);
    for ss = 1:n_SNR
        for gg = 1:n_gain_Voltage
            cr = floor((histo_clip_count{ff}(ss, gg) - min_ccb) / span_ccb*(ls-1)+1);
            vr = floor((histo_Varianz{ff}(ss, gg) - min_vv) / span_vv*(ls-1)+1);
            kr = floor((histo_Kurtosis{ff}(ss, gg) - min_kk) / span_kk*(ls-1)+1);
            lut{ff}(cr,vr,kr) = sollamp{ff}(ss, gg);
        end
    end
end


% separate Schleife notwendig, aufgrund der nicht-Eindeutigkeit der Indices
for ff = 1:n_Freq
    lut_cnt{ff} = zeros(ls,ls,ls);
    for ss = 1:n_SNR
        for gg = 1:n_gain_Voltage
            cr = floor((histo_clip_count{ff}(ss, gg) - min_ccb) / span_ccb*(ls-1)+1);
            vr = floor((histo_Varianz{ff}(ss, gg) - min_vv) / span_vv*(ls-1)+1);
            kr = floor((histo_Kurtosis{ff}(ss, gg) - min_kk) / span_kk*(ls-1)+1);
            imlut{ff}(ss,gg) = lut{ff}(cr,vr,kr);
            lut_cnt{ff}(cr, vr, kr) = lut_cnt{ff}(cr, vr, kr) + 1;

            % Korrektur via lut
            korr = lut{ff}(cr, vr, kr);
            % disp("f = " + freq_list(ff) + " | SNR: " + snr(ss) + " | Gain: " + round(gain_V(gg), 2)+ " | KorFaktor: " + korr + " | " + cr + " | " + vr + " | " + kr + " | ")
            UU = U_cmplx_out{ss, gg}(ff, :) .* korr;
            II = I_cmplx_out{ss, gg}(ff, :);
            Z_corr{ss, gg}(ff) = UU(f_freq_idx)./II(f_freq_idx)*1E6;
        end
    end
end

%%
% plot
ff = 1;
yytikz = [];
xxtikz = [];
for i = 1:length(snr)
    if mod(i, 3)
        yytikz = [yytikz, round(snr(i), 2)];
        xxtikz = [xxtikz, round(gain_V(i), 2)];
    end
end

figure(13)
tiledlayout(2,3,"TileSpacing","tight","Padding","tight")

nexttile;
cm = colormap(parula);
% cm(:,3) = linspace(0,1,256);
% cm(:,2) = linspace(1,0.6,256);
% cm(:,1) = linspace(1,0.3,256);

heatmap((histo_clip_count{ff}/(10000)*100));
ax = gca;
colormap(ax, cm);
ax.YData = round(snr, 2);
ax.XData = round(gain_V, 2);
xlabel("Gain");
ylabel("SNR [dB]");
title("Sättigungsgrad [%]");
colorbar

nexttile;
heatmap((histo_Varianz{ff}));
ax = gca;
colormap(ax, cm);
ax.YData = round(snr, 2);
ax.XData = round(gain_V, 2);
title(["Signalvarianz", "ohne Datenpunkte in Sättigung"]);
xlabel("Gain");
ylabel("SNR [dB]");
colorbar

nexttile;
heatmap((histo_Kurtosis{ff}));
ax = gca;
colormap(ax, cm);
ax.YData = round(snr, 2);
ax.XData = round(gain_V, 2);
title(["Signalkurtosis", "ohne Datenpunkte in Sättigung"]);
xlabel("Gain");
ylabel("SNR [dB]");
colorbar


nexttile;
h = heatmap(sollamp{ff});
h.CellLabelFormat = "%.4f";
h.ColorScaling = 'scaledrows';
% h.ColorLimits = [0 max(max(sollamp{ff}))];
ax = gca;
colormap(ax, cm);
ax.XData = round(gain_V, 2);
ax.YData = round(snr, 2);
xlabel("Gain");
ylabel("SNR [dB]");
title("Referenz - AKF (normiert)");
colorbar

nexttile;
h = heatmap(imlut{ff});
h.CellLabelFormat = "%.2f";
ax = gca;
h.ColorScaling = 'scaledrows';
% h.ColorLimits = [0 max(max(imlut{ff}))];
colormap(ax, cm);
ax.YData = round(snr, 2);
ax.XData = round(gain_V, 2);
title("AKF (normiert)");
xlabel("Gain");
ylabel("SNR [dB]");
colorbar

nexttile;
h = heatmap((imlut{ff}-sollamp{ff}));
h.CellLabelFormat = "%.2f";
% h.ColorScaling = 'scaledrows';
h.Colormap = sky;
ax = gca;
ax.YData = round(snr, 2);
ax.XData = round(gain_V, 2);
title("Differenz");
xlabel("Gain");
ylabel("SNR [dB]");
fontsize(gcf,7,"pixels")
colorbar
exportgraphics(gcf, "img/lut.pdf")
%% EIS Plot
figure(4)

tiledlayout('flow', 'TileSpacing','compact','Padding','compact')

gain_max = 185;
gain_min = 155;
snr_max = 20;
snr_min = 9;
cc = floor(n_SNR - snr_fokus);
for n = 0:0
    cc = cc + n;
    
    nexttile
    hold on
    for gg=1:1:n_gain_Voltage
        if gain_V(gg) < gain_max && gain_V(gg) > gain_min
            plot(real(Z{cc,gg}(:)),-imag(Z{cc,gg}(:)), 'DisplayName',strcat("Voltage Gain " + gain_V(gg)), "Marker", 'x')
        end
    end
    plot(real(Z2_clean{cc}(:)),-imag(Z2_clean{cc}(:)),'o--', 'DisplayName', 'No noise, no quantization etc.')
    title("EIS | SNR: " + snr(cc) + " dB")
    xlabel('Real{Z} (\mu\Omega)')
    ylabel('Imag{Z} (\mu\Omega)')
    lg = legend;
    lg.Location = 'northoutside';
    lg.NumColumns = 2;
    grid on; box on
    
    nexttile
    hold on;
    for gg=1:n_gain_Voltage
        if gain_V(gg) < gain_max && gain_V(gg) > gain_min
            plot(real(Z_corr{cc,gg}(:)),-imag(Z_corr{cc,gg}(:)),'x-', 'DisplayName',strcat("Voltage Gain " + gain_V(gg)))
        end
    end
    
    plot(real(Z2_clean{cc}(:)),-imag(Z2_clean{cc}(:)),'o--', 'DisplayName', 'No noise, no quantization etc.')
    title("EIS (AKF) | SNR: " + snr(cc) + " dB")
    xlabel('Real{Z_{k}} (\mu\Omega)')
    ylabel('Imag{Z_{k}} (\mu\Omega)')
    lg = legend;
    lg.Location = 'northoutside';
    lg.NumColumns = 2;
    grid on; box on
    
    nexttile
    hold on
    title("EIS (with AKF)")
    for ss=1:n_SNR
        for gg=1:n_gain_Voltage
            RMSE_clean(ss, gg) = rmse(Z{ss,gg}(:), Z2_clean{ss}(:));
            RMSE_corr(ss, gg) = rmse(Z_corr{ss,gg}(:), Z2_clean{ss}(:));
        end
        if snr(ss) < snr_max && snr(ss) > snr_min
            semilogy(gain_V, RMSE_corr(ss, :), 'x-', 'DisplayName', "Korrigiert@" + string(snr(ss)) + " dB" )
            semilogy(gain_V, RMSE_clean(ss, :), 'DisplayName', "Referenz@" + string(snr(ss)) + " dB" )
        end
    end
    semilogy(gain_V, RMSE_clean(1, :), 'o--', 'DisplayName', "No noise, no quantization etc." ,'LineWidth',1)
    
    xlabel('Voltage Gain')
    ylabel('RMSE (\mu\Omega)')
    title("RMSE")
    lg = legend;
    set(gca, 'YScale', 'log')
    lg.Location = 'northoutside';
    lg.FontSize = 4;
    lg.NumColumns = 2;
    grid on; box on
    hold off

    nexttile
    hold on
    for ss=1:n_SNR
        for gg=1:n_gain_Voltage
            cutoff_values_Voltage_percent_plot(ss,gg) = cutoff_values_Voltage_percent{ss,gg}(1);
            cutoff_values_Voltage_percent_plot_clean(ss,gg) = cutoff_values_Voltage_percent_clean{ss,gg}(1);
        end
        if snr(ss) < snr_max && snr(ss) > snr_min
            plot(gain_V, cutoff_values_Voltage_percent_plot(ss,:),'x-', 'DisplayName', strcat("SNR = " + snr(ss)+ " dB"))
        end 
    end
    plot(gain_V, cutoff_values_Voltage_percent_plot_clean(1,:),'o--', 'DisplayName', "No noise, no quantization etc." ,'LineWidth', 1)
    
    xlabel('Voltage Gain')
    ylabel('Cutoff Values (<0V or >3.3V) (%)')
    title("Cutoff Values")
    lg = legend;
    lg.Location = 'northoutside';
    lg.NumColumns = 2;
    fontsize(gcf,7,"pixels")
    grid on; box on
end
exportgraphics(gcf, "img/ergebnisse.pdf")
function [clipp_count, signal_clean, signal_clipped] = clean_clipped(signal, min_value, max_value)
    ir = 1;
    clipp_count = 0;
    signal_clean = zeros(size(signal)); % Initialisieren Sie signal_clean vor der Verwendung
    for i = 1:length(signal)
        if signal(i) > max_value
            signal_clipped(i) = max_value;
            clipp_count = clipp_count + 1;
        elseif signal(i) < min_value
            signal_clipped(i) = min_value;
            clipp_count = clipp_count + 1;
        else
            signal_clipped(i) = signal(i); % Behalten Sie das Originalsignal bei
            signal_clean(ir) = signal(i);
            ir = ir + 1;
        end
    end
    signal_clean = signal_clean(1:ir-1); % Kürzen Sie signal_clean auf die tatsächliche Länge
end




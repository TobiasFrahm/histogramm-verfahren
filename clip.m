function [cutoff_values, signal] = clip(signal, signalMinVoltage, signalMaxVoltage)
            below_min = sum(signal < signalMinVoltage);
            above_max = sum(signal > signalMaxVoltage);
            total_length = length(signal);
            cutoff_values = (below_min + above_max) / total_length * 100;
            signal(signal<signalMinVoltage) = signalMinVoltage;   
            signal(signal>signalMaxVoltage) = signalMaxVoltage;    
end


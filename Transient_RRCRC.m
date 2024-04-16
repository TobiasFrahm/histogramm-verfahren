function [voltage_response] = Transient_RRCRC(Fs, current, R0, R1, C1, R2, C2)
% Transiente Simulation R-RC-RC  - Thevenin-Struktur 2. Ordnung
% einfaches Modell Relaxation + Ladung einer Batterie
% gegeben ibat(t) ges. Vbat(i,t)
% % v0.3, Nov. 2022, K.-R. Riemschneider
%
% ---->----R1----R2||C2---R3||C3------->--------
%                                   |
%   ibat(t)          	     Vbat(i,t)
%                                   |
%                                   ---
    vOCV_Ch = 4.2;  % voltage range for OCV curve
    vOCV_DCh = 2.8;
    
    vCutOff_Ch = 4.1;  % cut off voltages
    vCutOff_DCh = 2.9;

    QAh = 4.5;   % battery capacitance; Ah
    QAs = QAh * 3600;
    SOCStart = 60;
    
    tau1 = R1*C1;
    tau2 = R2*C2;
    T = 1/Fs;
    dT = T;

    % model OCV(voltage), very simple
    n_SOC = length(current);
    x = -50:100/n_SOC:50;
    QOCV = 0:100/n_SOC:100;
    vOCV = (x.^3/50^3)*0.7 + 0.7 + 2.8;

    
    % assumption, same for i(Rx)
    t = 1;
    iR1(t) = 0;
    iR2(t) = 0;

    % boundary condition, start at SOC-based charged battery
    vOCV_Initial = vOCV(QOCV==SOCStart);
    Charge(t) = QOCV(QOCV==SOCStart);

    % save OCV per time step
    vOCV_Test = vOCV_Initial;
    vOCV_Step = vOCV_Initial;   % no change in vOCV

    iR1 = zeros(size(current));
    iR2 = zeros(size(current));
    Charge = zeros(size(current));
    vOCV_Test = zeros(size(current));

    % Transiente Simulation
    for t = 2:length(current)
        % R0
        % vR0 = ibat(t) * R0;
        
        % % RC1 Glied R1|C1
        % iC1 = ibat(t) - iR1;
        % Q1 = Q1 + iC1;
        % vR1C1 = 1 /C1 * Q1;
        iR1(t) = exp(-dT/tau1)*iR1(t-1) + (1-exp(-dT/tau1))*current(t);
        
        % RC2 R2|C2
        % iC2 = ibat(t) - iR2;
        % Q2 = Q2 + iC2;
        % vR2C2 = 1 /C2 * Q2;
        % iR2 = vR2C2 /R2;
        iR2(t) = exp(-dT/tau2)*iR2(t-1) + (1-exp(-dT/tau2))*current(t);
        
        % charge; search for Q in QOCV and get vOCV from that finding as vOCV_Step
        Charge(t) = Charge(t-1) + current(t) * dT / QAs * 100;     % normalized to 100

        % get corresponding OCV
        vOCV_Test(t) = vOCV_Step;
        
        % terminal voltage
        voltage_response(t) = vOCV_Step + (R0*current(t) + R1*iR1(t) + R2*iR2(t));
        
        if isnan(voltage_response(t))
            voltage_response(t) = voltage(t-1);
        end
    end
end
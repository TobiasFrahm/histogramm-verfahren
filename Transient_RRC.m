function [voltage_response] = Transient_RRC(Fs, current, R0, R1, C1)
    vOCV_Ch = 4.2;  % voltage range for OCV curve
    vOCV_DCh = 2.8;
    
    vCutOff_Ch = 4.1;  % cut off voltages
    vCutOff_DCh = 2.9;
    
    
    flagChDCh = true;               % charging(true) or discharging(false)?
    QAh = 4.5;   % battery capacitance; Ah
    QAs = QAh * 3600;
    SOCStart = 60;
    
    % single P45B
    tau1 = R1*C1;
    T = 1/Fs;
    dT = T;

    % model OCV(voltage), very simple
    n_SOC = length(current);
    x = -50:100/n_SOC:50;
    QOCV = 0:100/n_SOC:100;
    vOCV = (x.^3/50^3)*0.7 + 0.7 + 2.8;



    % Charge/Discharge cycles

    % first time step
    t = 1;
    
    % assumption, same for i(Rx)
    iR1(t) = 0;

    % boundary condition, start at SOC-based charged battery
    vOCV_Initial = vOCV(QOCV==SOCStart);
    Charge(t) = QOCV(QOCV==SOCStart);

    % save OCV per time step
    vOCV_Test = vOCV_Initial;
    vOCV_Step = vOCV_Initial;   % no change in vOCV

    % voltage
    voltage_response(t) = vOCV_Initial + ( R0*current(t) + R1*(1-exp(-dT/tau1))*current(t) );
    
    iR1 = zeros(size(current));
    Charge = zeros(size(current));
    vOCV_Test = zeros(size(current));
     
    for t = 2:length(current)
        if (flagChDCh && voltage_response(t-1) <  vCutOff_Ch) || (~flagChDCh && voltage_response(t-1) >  vCutOff_DCh)
   
            % current
            iR1(t) = exp(-dT/tau1)*iR1(t-1) + (1-exp(-dT/tau1))*current(t);
    
            % charge; search for Q in QOCV and get vOCV from that finding as vOCV_Step
            Charge(t) = Charge(t-1) + current(t)*dT/QAs*100;     % normalized to 100
    
            % get corresponding OCV
            vOCV_Test(t) = vOCV_Step;
    
            % terminal voltage
            voltage_response(t) = vOCV_Step + (R0*current(t) + R1*iR1(t));

            if isnan(voltage_response(t))
                voltage_response(t) = voltage(t-1);
            end
        end
    end
end


function faceIR = planetIRsimple(param)
    % Constants
    sigma = 5.67e-8; % Stefan-Boltzmann constant (W · m -2 · K -4 )

    %%% Get the satellite half FOV to earth tangent (limb angle: alpha in the schematic)
    alpha = param.orb.prop.limbAngle; % deg

    %%% Get theta angle
    theta = param.orb.prop.albedoAngle.AlbedoAngle_deg_; % deg

    %%% Compute averaged T^4 (see schematic to understand the weighting method)
    dayArcWeight = max((180 - theta - alpha),0);
    nightArcWeight = max((theta - alpha),0);
    Tavg_power4 = (dayArcWeight * (param.Tir_Day)^4 + nightArcWeight * (param.Tir_Night)^4)./(dayArcWeight + nightArcWeight);

    % Compute planetary IR flux
    Pir = sigma * Tavg_power4; % Central body radiation


    % Compute IR flux at each face
    IR = zeros(length(param.pwr.viewFactor(:,1)), length(param.facesMaterial));
    for n=1:length(param.facesMaterial)
        IR(:,n) = param.pwr.viewFactor(:,n) .* Pir * param.emissivity(n) * param.facesArea(n);
    end

    % Transform it to timetable for each face
    time = param.orb.prop.time_Epoch;
    % X+
    faceIR.Xplus = array2table(IR(:,1));
    faceIR.Xplus.time = seconds(time); 
    faceIR.Xplus = table2timetable(faceIR.Xplus);
    % X-
    faceIR.Xminus = array2table(IR(:,2));
    faceIR.Xminus.time = seconds(time); 
    faceIR.Xminus = table2timetable(faceIR.Xminus);
    % Y+
    faceIR.Yplus = array2table(IR(:,3));
    faceIR.Yplus.time = seconds(time); 
    faceIR.Yplus = table2timetable(faceIR.Yplus);
    % Y-
    faceIR.Yminus = array2table(IR(:,4));
    faceIR.Yminus.time = seconds(time); 
    faceIR.Yminus = table2timetable(faceIR.Yminus);
    % Z+
    faceIR.Zplus = array2table(IR(:,5));
    faceIR.Zplus.time = seconds(time); 
    faceIR.Zplus = table2timetable(faceIR.Zplus);
    % Z-
    faceIR.Zminus = array2table(IR(:,6));
    faceIR.Zminus.time = seconds(time); 
    faceIR.Zminus = table2timetable(faceIR.Zminus);
    
end
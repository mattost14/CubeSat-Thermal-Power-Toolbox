function qIR = qIRsubsolar(Latss, albedo, solarFlux, moonEmissivity, Tdark)
% The moon's surface temperature distribution is computed by the equation 
% provided in the Human Landing System Lunar Thermal Analysis Guidebook, 
% volume HLS-UG-001. NASA, 2021. (Pag. 17)
% Latss : degrees latitude from subsolar point (Subsolar Coordinate System)
% albedo: Lunar surface albedo
% solarFlux: solar flux (W/m2)
% moonEmissivity: Lunar Surface Emissivity
% Tdark: unilluminated lunar surface temperature

    sigma = 5.670374419e-8; % W*m2*K^-4
    if(Latss >= 0)
        qIR = sind(Latss)*((1-albedo)*solarFlux - sigma * moonEmissivity * Tdark^4)+  sigma * moonEmissivity * Tdark^4;
    else
        qIR = sigma * moonEmissivity * Tdark^4;
    end
end

function faceSolar = computeSolarRadiation(param)
    solarFlux = param.orb.solarFlux;

    % Compute solar flux at each face
    for n=1:length(param.facesMaterial)

        % Check if user wants to consider shadow losses)
        iluminatedAreaRatio = 1;
        if(param.computePwrWithShadowFlag && param.pwr.shadow.flag) 
            iluminatedAreaRatio = 1 - param.pwr.shadow.shadowFaceAreaRatio(:,n);
        end

        if(param.facesMaterial(n)=="Solar Panel") % If face is solar cell, check whether it is generating electrical energy or it is dissipating outside 
            effectiveAbsoptivity =  min(param.absorptivity(n) - param.pwr.electricalEfficiency * param.solarCellEff(n),0);
            solar(:,n) = iluminatedAreaRatio .* param.pwr.solarLight.(n) .* solarFlux .* effectiveAbsoptivity * param.facesArea(n);
        else
            solar(:,n) = iluminatedAreaRatio .* param.pwr.solarLight.(n) .* solarFlux * param.absorptivity(n) * param.facesArea(n);
        end
    end

    % Transform it to timetable for each face
    time = param.orb.prop.time_Epoch;
    % X+
    faceSolar.Xplus = array2table(solar(:,1));
    faceSolar.Xplus.time = seconds(time); 
    faceSolar.Xplus = table2timetable(faceSolar.Xplus);
    % X-
    faceSolar.Xminus = array2table(solar(:,2));
    faceSolar.Xminus.time = seconds(time); 
    faceSolar.Xminus = table2timetable(faceSolar.Xminus);
    % Y+
    faceSolar.Yplus = array2table(solar(:,3));
    faceSolar.Yplus.time = seconds(time); 
    faceSolar.Yplus = table2timetable(faceSolar.Yplus);
    % Y-
    faceSolar.Yminus = array2table(solar(:,4));
    faceSolar.Yminus.time = seconds(time); 
    faceSolar.Yminus = table2timetable(faceSolar.Yminus);
    % Z+
    faceSolar.Zplus = array2table(solar(:,5));
    faceSolar.Zplus.time = seconds(time); 
    faceSolar.Zplus = table2timetable(faceSolar.Zplus);
    % Z-
    faceSolar.Zminus = array2table(solar(:,6));
    faceSolar.Zminus.time = seconds(time); 
    faceSolar.Zminus = table2timetable(faceSolar.Zminus);
end

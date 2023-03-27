function faceAlbedo = computeAlbedoRadiation(param)
    solarFlux = param.orb.solarFlux;
    albedoFactor = param.Albedo;
    
    % Compute solar flux at each face
    for n=1:length(param.facesMaterial)
        % Albedo view factor (Falb) = Fir * cos(SZA)
        Falb = param.pwr.viewFactor(:,n) .* max(0, cosd(param.orb.prop.albedoAngle.(1)));

        if(param.facesMaterial(n)=="Solar Panel") % If face is solar cell, check whether it is generating electrical energy or it is dissipating outside 
            effectiveAbsoptivity =  param.absorptivity(n) - param.pwr.electricalEfficiency * param.solarCellEff(n);
            albedoRadiation(:,n) = albedoFactor * Falb .* solarFlux .* effectiveAbsoptivity * param.facesArea(n);
        else
            albedoRadiation(:,n) = albedoFactor * Falb .* solarFlux * param.absorptivity(n) * param.facesArea(n);
        end
    end

    % Transform it to timetable for each face
    time = param.orb.prop.time_Epoch;
    % X+
    faceAlbedo.Xplus = array2table(albedoRadiation(:,1));
    faceAlbedo.Xplus.time = seconds(time); 
    faceAlbedo.Xplus = table2timetable(faceAlbedo.Xplus);
    % X-
    faceAlbedo.Xminus = array2table(albedoRadiation(:,2));
    faceAlbedo.Xminus.time = seconds(time); 
    faceAlbedo.Xminus = table2timetable(faceAlbedo.Xminus);
    % Y+
    faceAlbedo.Yplus = array2table(albedoRadiation(:,3));
    faceAlbedo.Yplus.time = seconds(time); 
    faceAlbedo.Yplus = table2timetable(faceAlbedo.Yplus);
    % Y-
    faceAlbedo.Yminus = array2table(albedoRadiation(:,4));
    faceAlbedo.Yminus.time = seconds(time); 
    faceAlbedo.Yminus = table2timetable(faceAlbedo.Yminus);
    % Z+
    faceAlbedo.Zplus = array2table(albedoRadiation(:,5));
    faceAlbedo.Zplus.time = seconds(time); 
    faceAlbedo.Zplus = table2timetable(faceAlbedo.Zplus);
    % Z-
    faceAlbedo.Zminus = array2table(albedoRadiation(:,6));
    faceAlbedo.Zminus.time = seconds(time); 
    faceAlbedo.Zminus = table2timetable(faceAlbedo.Zminus);
end

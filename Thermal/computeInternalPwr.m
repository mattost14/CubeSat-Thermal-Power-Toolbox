function internalPwr = computeInternalPwr(param)
    % This function computes the internal dissipation power profile based
    % on inputs in the power tab

    time = param.orb.prop.time_Epoch;
    pwrDissipationDay = param.pwrDissipationDaySimpleModel;
    pwrDissipationNight = param.pwrDissipationNightSimpleModel;
    sunMagnitude = param.orb.prop.sunMagnitude;
    switch param.pwrDissipationProfile
        case "Constant"
            internalPwr = pwrDissipationDay * ones(length(time),1);
        case "Day/Night"
            internalPwr = pwrDissipationDay .* sunMagnitude;
            internalPwr = internalPwr +  pwrDissipationNight .* (1-sunMagnitude);
    end

    % Internal power as timetable
    internalPwr = array2table(internalPwr);
    internalPwr.time = seconds(time); 
    internalPwr = table2timetable(internalPwr);
end

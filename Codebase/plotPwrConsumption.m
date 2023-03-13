function plotPwrConsumption(appAxis, time, sunMagnitude, pwrDissipationDay, pwrDissipationNight, pwrDissipationProfile)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here
    cla(appAxis)
    switch pwrDissipationProfile
        case "Constant"
            pwrDissipation = pwrDissipationDay * ones(1, length(time));
            plot(appAxis, time/3600, pwrDissipation, 'b', 'DisplayName',"Node 1");
        case "Day/Night"
            pwrDissipation = pwrDissipationDay .* sunMagnitude;
            pwrDissipation = pwrDissipation +  pwrDissipationNight .* (1-sunMagnitude);
            plot(appAxis, time/3600, pwrDissipation, 'b', 'DisplayName',"Node 1");
    end


    % Compute average
    avgPwr = mean(pwrDissipation);

    xL=appAxis.XLim;
    yL=appAxis.YLim;
    text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgPwr,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");

    
    legend(appAxis)
    xlabel(appAxis, "Epoch (hour)")
    ylabel(appAxis, "Power (W)")
    grid(appAxis,"on")
end


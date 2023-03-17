function plotSunlight(appAxis, time, sunMagnitude)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here
    cla(appAxis);
    plot(appAxis, time/3600, sunMagnitude, 'b');
    ylim(appAxis, [0,1])
    xlim(appAxis, "auto")
    xlabel(appAxis, "Epoch (hour)")
    yticks(appAxis, [0 1])
    grid(appAxis,"on")
end


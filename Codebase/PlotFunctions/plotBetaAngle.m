function plotBetaAngle(appAxis, time, betaAngle)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here
    plot(appAxis, time/3600, betaAngle, 'b');
    xlabel(appAxis, "Epoch (hour)")
    ylabel(appAxis, "\beta (°)")
    grid(appAxis,"on")
    ytickformat(appAxis, '%.1f')
end


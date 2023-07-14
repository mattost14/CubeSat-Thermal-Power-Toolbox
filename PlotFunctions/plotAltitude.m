function plotAltitude(appAxis, time, altitude)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here
    plot(appAxis, time/3600, altitude, 'b',"LineWidth",2);
    xlabel(appAxis, "Epoch (hour)")
    ylabel(appAxis, "h (km)")
    grid(appAxis,"on")
    ytickformat(appAxis, '%.1f')
end


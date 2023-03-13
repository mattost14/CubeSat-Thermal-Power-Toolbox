function plotInternalTemperatures(appAxis, thermal)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here

    time = thermal.out.internalNodeTemp.Time;
    temp_C = thermal.out.internalNodeTemp.Data(:)-273;

    cla(appAxis)
    plot(appAxis, time/3600, temp_C, 'b', "DisplayName", "Node 1");

    legend(appAxis);
    ylabel(appAxis, "Temp (oC)")
    xlabel(appAxis, "Epoch (hour)")
    grid(appAxis,"on")
    ytickformat(appAxis, '%.1f')
end
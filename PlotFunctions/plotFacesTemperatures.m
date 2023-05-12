function plotFacesTemperatures(appAxis, thermal)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here

    time = thermal.out.XplusTemp.Time;
    tempXplus_C = thermal.out.XplusTemp.Data(:)-273;
    tempXminus_C = thermal.out.XminusTemp.Data(:)-273;
    tempYplus_C = thermal.out.YplusTemp.Data(:)-273;
    tempYminus_C = thermal.out.YminusTemp.Data(:)-273;
    tempZplus_C = thermal.out.ZplusTemp.Data(:)-273;
    tempZminus_C = thermal.out.ZminusTemp.Data(:)-273;

    cla(appAxis)
    hold(appAxis, "on")
    plot(appAxis, time/3600, tempXplus_C, "DisplayName", "X+", "Color","red");
    plot(appAxis, time/3600, tempXminus_C, "DisplayName", "X-", "Color","green");
    plot(appAxis, time/3600, tempYplus_C,  "DisplayName", "Y+", "Color","blue");
    plot(appAxis, time/3600, tempYminus_C,  "DisplayName", "Y-", "Color","cyan");
    plot(appAxis, time/3600, tempZplus_C,  "DisplayName", "Z+", "Color","magenta");
    plot(appAxis, time/3600, tempZminus_C,  "DisplayName", "Z-", "Color","#EDB120");

    legend(appAxis,'ItemHitFcn',@toggleLegend);
    ylabel(appAxis, "Temp (oC)")
    xlabel(appAxis, "Epoch (hour)")
    grid(appAxis,"on")
    ytickformat(appAxis, '%.1f')
end
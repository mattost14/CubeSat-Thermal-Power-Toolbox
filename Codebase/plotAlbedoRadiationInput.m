function plotAlbedoRadiationInput(appAxis, thermal, plotAllFacesFlag)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here

    % Compute total
    time = seconds(thermal.faceAlbedo.Xplus.time);
    totalAlbedoRadiation =   thermal.faceAlbedo.Xplus.(1) +...
                thermal.faceAlbedo.Xminus.(1) +...
                thermal.faceAlbedo.Yplus.(1) +...
                thermal.faceAlbedo.Yminus.(1) +...
                thermal.faceAlbedo.Zplus.(1) +...
                thermal.faceAlbedo.Zminus.(1);

    % Compute average
    avgRadiation = mean(totalAlbedoRadiation);

    cla(appAxis)
    if(~plotAllFacesFlag)
        plot(appAxis, time/3600, totalAlbedoRadiation, 'b', "DisplayName", "Total");
    else
        hold(appAxis,"on");
        plot(appAxis, time/3600, thermal.faceAlbedo.Xplus.(1), "DisplayName", "X+", "Color","red");
        plot(appAxis, time/3600, thermal.faceAlbedo.Xminus.(1), "DisplayName", "X-", "Color","green");
        plot(appAxis, time/3600, thermal.faceAlbedo.Yplus.(1),  "DisplayName", "Y+", "Color","blue");
        plot(appAxis, time/3600, thermal.faceAlbedo.Yminus.(1),  "DisplayName", "Y-", "Color","cyan");
        plot(appAxis, time/3600, thermal.faceAlbedo.Zplus.(1),  "DisplayName", "Z+", "Color","magenta");
        plot(appAxis, time/3600, thermal.faceAlbedo.Zminus.(1),  "DisplayName", "Z-", "Color","#EDB120");

    end
    legend(appAxis);
    ylabel(appAxis, "Flux (W)")
    xlabel(appAxis, "Epoch (hour)")
    grid(appAxis,"on")
    ytickformat(appAxis, '%.1f')

    xL=appAxis.XLim;
    yL=appAxis.YLim;
    text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");

end
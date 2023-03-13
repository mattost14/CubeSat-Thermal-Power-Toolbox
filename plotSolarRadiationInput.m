function plotSolarRadiationInput(appAxis, thermal, plotAllFacesFlag)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here

    % Compute total
    time = seconds(thermal.faceSolar.Xplus.time);
    totalSolar =   thermal.faceSolar.Xplus.(1) +...
                thermal.faceSolar.Xminus.(1) +...
                thermal.faceSolar.Yplus.(1) +...
                thermal.faceSolar.Yminus.(1) +...
                thermal.faceSolar.Zplus.(1) +...
                thermal.faceSolar.Zminus.(1);

    % Compute average
    avgRadiation = mean(totalSolar);

    cla(appAxis)
    if(~plotAllFacesFlag)
        plot(appAxis, time/3600, totalSolar, 'b', "DisplayName", "Total");
    else
        hold(appAxis,"on");
        plot(appAxis, time/3600, thermal.faceSolar.Xplus.(1), "DisplayName", "X+", "Color","red");
        plot(appAxis, time/3600, thermal.faceSolar.Xminus.(1), "DisplayName", "X-", "Color","green");
        plot(appAxis, time/3600, thermal.faceSolar.Yplus.(1),  "DisplayName", "Y+", "Color","blue");
        plot(appAxis, time/3600, thermal.faceSolar.Yminus.(1),  "DisplayName", "Y-", "Color","cyan");
        plot(appAxis, time/3600, thermal.faceSolar.Zplus.(1),  "DisplayName", "Z+", "Color","magenta");
        plot(appAxis, time/3600, thermal.faceSolar.Zminus.(1),  "DisplayName", "Z-", "Color","#EDB120");
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
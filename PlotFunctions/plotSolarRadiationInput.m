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
        plot(appAxis, time/3600, totalSolar, 'b', "DisplayName", "Total","LineWidth",2);
    else
        hold(appAxis,"on");
        plot(appAxis, time/3600, thermal.faceSolar.Xplus.(1), "DisplayName", "X+", "Color","red","LineWidth",2);
        plot(appAxis, time/3600, thermal.faceSolar.Xminus.(1), "DisplayName", "X-", "Color","green","LineWidth",2);
        plot(appAxis, time/3600, thermal.faceSolar.Yplus.(1),  "DisplayName", "Y+", "Color","blue","LineWidth",2);
        plot(appAxis, time/3600, thermal.faceSolar.Yminus.(1),  "DisplayName", "Y-", "Color","cyan","LineWidth",2);
        plot(appAxis, time/3600, thermal.faceSolar.Zplus.(1),  "DisplayName", "Z+", "Color","magenta","LineWidth",2);
        plot(appAxis, time/3600, thermal.faceSolar.Zminus.(1),  "DisplayName", "Z-", "Color","#EDB120","LineWidth",2);
    end
    legend(appAxis,'ItemHitFcn',@toggleLegend);
    ylabel(appAxis, "Flux (W)")
    xlabel(appAxis, "Epoch (hour)")
    grid(appAxis,"on")
    ytickformat(appAxis, '%.1f')
    
    xL=appAxis.XLim;
    yL=appAxis.YLim;
    text(appAxis, 0.5,1,strcat("\mu = ", num2str(avgRadiation,2), " W"),'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");

%     text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");

end
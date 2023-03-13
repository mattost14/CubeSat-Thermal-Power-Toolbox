function plotBatteryDepletion(appAxis, pwr)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here
    cla(appAxis)
    plot(appAxis, pwr.time_Epoch/3600, -pwr.energyDrawnFromBat_Wh, 'b');
%     legend(appAxis)
    xlabel(appAxis, "Epoch (hour)")
    ylabel(appAxis, "Energy (Wh)")
    grid(appAxis,"on")

    % Compute average
    maxDepletion = max(pwr.energyDrawnFromBat_Wh);

    xL=appAxis.XLim;
    yL=appAxis.YLim;
    text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("max = ", num2str(maxDepletion,2), " Wh"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");


end
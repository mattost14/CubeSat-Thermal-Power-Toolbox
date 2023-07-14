function plotPwrGeneration(appAxis, pwr)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here
    cla(appAxis)
    plot(appAxis, pwr.time_Epoch/3600, pwr.generatedTotalPower, 'b',"LineWidth",2);
%     legend(appAxis)
    xlabel(appAxis, "Epoch (hour)")
    ylabel(appAxis, "Power (W)")
    grid(appAxis,"on")

    % Compute average
    avgPwr = mean(pwr.generatedTotalPower);

    xL=appAxis.XLim;
    yL=appAxis.YLim;
    text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgPwr,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");


end

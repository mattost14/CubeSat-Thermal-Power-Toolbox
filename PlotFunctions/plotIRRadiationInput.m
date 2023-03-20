function plotIRRadiationInput(appAxis, thermal, plotAllFacesFlag, plotAbsorbedIRFlag, plotRadiatedIRFlag)
%PLOTSUNLIGHT Summary of this function goes here
%   Detailed explanation goes here

    % Compute total absorbed/radiated IR
    timeAbsorbedIR = seconds(thermal.faceIR.Xplus.time);
    totalAbsorbedIR =   thermal.faceIR.Xplus.(1) +...
                thermal.faceIR.Xminus.(1) +...
                thermal.faceIR.Yplus.(1) +...
                thermal.faceIR.Yminus.(1) +...
                thermal.faceIR.Zplus.(1) +...
                thermal.faceIR.Zminus.(1);

    timeRadiatedIR = thermal.out.XplusIRfluxOut.Time;
    totalRadiatedIR = thermal.out.XplusIRfluxOut.Data(:) +...
                thermal.out.XminusIRfluxOut.Data(:) +...
                thermal.out.YplusIRfluxOut.Data(:) +...
                thermal.out.YminusIRfluxOut.Data(:) +...
                thermal.out.ZplusIRfluxOut.Data(:) +...
                thermal.out.ZminusIRfluxOut.Data(:);

    % Compute averages
    avgAbsorbedRadiation = mean(totalAbsorbedIR);
    avgRadiatedRadiation = mean(totalRadiatedIR);



    cla(appAxis)
    if(~plotAllFacesFlag)
        if(plotAbsorbedIRFlag && ~plotRadiatedIRFlag) % Only Absorbed IR
            plot(appAxis, timeAbsorbedIR/3600, totalAbsorbedIR, 'b', "DisplayName", "Total absorbed IR");
            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgAbsorbedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        elseif(~plotAbsorbedIRFlag && plotRadiatedIRFlag) % Only Radiate IR
            hold(appAxis,"on");
            plot(appAxis, timeRadiatedIR/3600, totalRadiatedIR, 'r', "DisplayName", "Total radiated IR");
            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgRadiatedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        else  % Both
            plot(appAxis, timeAbsorbedIR/3600, totalAbsorbedIR, 'b', "DisplayName", "Total absorbed IR");
            hold(appAxis,"on");
            plot(appAxis, timeRadiatedIR/3600, totalRadiatedIR, 'r', "DisplayName", "Total radiated IR");
            
            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu_{abs} = ", num2str(avgAbsorbedRadiation,2), " W; ", "\mu_{rad} = ", num2str(avgRadiatedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        end
    else
        if(plotAbsorbedIRFlag)
            hold(appAxis,"on");
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Xplus.(1), "DisplayName", "X+", "Color","red");
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Xminus.(1), "DisplayName", "X-", "Color","green");
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Yplus.(1),  "DisplayName", "Y+", "Color","blue");
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Yminus.(1),  "DisplayName", "Y-", "Color","cyan");
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Zplus.(1),  "DisplayName", "Z+", "Color","magenta");
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Zminus.(1),  "DisplayName", "Z-", "Color","#EDB120");

            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgAbsorbedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        else
            hold(appAxis,"on");
            plot(appAxis, timeRadiatedIR/3600, thermal.out.XplusIRfluxOut.Data(:), "DisplayName", "X+", "Color","red");
            plot(appAxis, timeRadiatedIR/3600, thermal.out.XminusIRfluxOut.Data(:), "DisplayName", "X-", "Color","green");
            plot(appAxis, timeRadiatedIR/3600, thermal.out.YplusIRfluxOut.Data(:),  "DisplayName", "Y+", "Color","blue");
            plot(appAxis, timeRadiatedIR/3600, thermal.out.YminusIRfluxOut.Data(:),  "DisplayName", "Y-", "Color","cyan");
            plot(appAxis, timeRadiatedIR/3600, thermal.out.ZplusIRfluxOut.Data(:),  "DisplayName", "Z+", "Color","magenta");
            plot(appAxis, timeRadiatedIR/3600, thermal.out.ZminusIRfluxOut.Data(:),  "DisplayName", "Z-", "Color","#EDB120");

            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu_{rad} = ", num2str(avgRadiatedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        end
    end
    legend(appAxis);
    ylabel(appAxis, "Flux (W)")
    xlabel(appAxis, "Epoch (hour)")
    grid(appAxis,"on")
    ytickformat(appAxis, '%.1f')

end


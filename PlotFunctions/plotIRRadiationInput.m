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
            plot(appAxis, timeAbsorbedIR/3600, totalAbsorbedIR, 'b', "DisplayName", "Total absorbed IR","LineWidth",2);
            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis, 0.5,1,strcat("\mu = ", num2str(avgAbsorbedRadiation,2), " W"),'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");

%             text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgAbsorbedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        elseif(~plotAbsorbedIRFlag && plotRadiatedIRFlag) % Only Radiate IR
            hold(appAxis,"on");
            plot(appAxis, timeRadiatedIR/3600, totalRadiatedIR, 'r', "DisplayName", "Total radiated IR","LineWidth",2);
            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis,0.5, 1,strcat("\mu = ", num2str(avgRadiatedRadiation,2), " W"),'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");

%             text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgRadiatedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        else  % Both
            plot(appAxis, timeAbsorbedIR/3600, totalAbsorbedIR, 'b', "DisplayName", "Total absorbed IR","LineWidth",2);
            hold(appAxis,"on");
            plot(appAxis, timeRadiatedIR/3600, totalRadiatedIR, 'r', "DisplayName", "Total radiated IR","LineWidth",2);
            
            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis, 0.5, 1,strcat("\mu_{abs} = ", num2str(avgAbsorbedRadiation,2), " W; ", "\mu_{rad} = ", num2str(avgRadiatedRadiation,2), " W"),'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");

%             text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu_{abs} = ", num2str(avgAbsorbedRadiation,2), " W; ", "\mu_{rad} = ", num2str(avgRadiatedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        end
    else
        if(plotAbsorbedIRFlag)
            hold(appAxis,"on");
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Xplus.(1), "DisplayName", "X+", "Color","red","LineWidth",2);
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Xminus.(1), "DisplayName", "X-", "Color","green","LineWidth",2);
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Yplus.(1),  "DisplayName", "Y+", "Color","blue","LineWidth",2);
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Yminus.(1),  "DisplayName", "Y-", "Color","cyan","LineWidth",2);
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Zplus.(1),  "DisplayName", "Z+", "Color","magenta","LineWidth",2);
            plot(appAxis, timeAbsorbedIR/3600, thermal.faceIR.Zminus.(1),  "DisplayName", "Z-", "Color","#EDB120","LineWidth",2);

            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis, 0.5, 1, strcat("\mu = ", num2str(avgAbsorbedRadiation,2), " W"),'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");

%             text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu = ", num2str(avgAbsorbedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        else
            hold(appAxis,"on");
            plot(appAxis, timeRadiatedIR/3600, thermal.out.XplusIRfluxOut.Data(:), "DisplayName", "X+", "Color","red","LineWidth",2);
            plot(appAxis, timeRadiatedIR/3600, thermal.out.XminusIRfluxOut.Data(:), "DisplayName", "X-", "Color","green","LineWidth",2);
            plot(appAxis, timeRadiatedIR/3600, thermal.out.YplusIRfluxOut.Data(:),  "DisplayName", "Y+", "Color","blue","LineWidth",2);
            plot(appAxis, timeRadiatedIR/3600, thermal.out.YminusIRfluxOut.Data(:),  "DisplayName", "Y-", "Color","cyan","LineWidth",2);
            plot(appAxis, timeRadiatedIR/3600, thermal.out.ZplusIRfluxOut.Data(:),  "DisplayName", "Z+", "Color","magenta","LineWidth",2);
            plot(appAxis, timeRadiatedIR/3600, thermal.out.ZminusIRfluxOut.Data(:),  "DisplayName", "Z-", "Color","#EDB120","LineWidth",2);

            xL=appAxis.XLim;
            yL=appAxis.YLim;
            text(appAxis, 0.5, 1,strcat("\mu_{rad} = ", num2str(avgRadiatedRadiation,2), " W"),'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");

%             text(appAxis, 1.2*mean(xL), 0.99*yL(2),strcat("\mu_{rad} = ", num2str(avgRadiatedRadiation,2), " W"),'HorizontalAlignment','right','VerticalAlignment','top','FontWeight','bold','BackgroundColor',"white");
        end
    end
    legend(appAxis,'ItemHitFcn',@toggleLegend);
    ylabel(appAxis, "Flux (W)")
    xlabel(appAxis, "Epoch (hour)")
    grid(appAxis,"on")
    ytickformat(appAxis, '%.1f')

end


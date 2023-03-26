function plotAnimatedPwrCubeSat(param, appAxis, t)
%PLOTANIMATEDCUBESAT Summary of this function goes here
%   Detailed explanation goes here

    [faceVertices_L, fixedDeployableVertices_L, trackingDeployableVertices_L, Sun_L] = getVerticesPoints(param);

    attitude.XplusFaceVector_LVLH = param.pwr.XplusFaceVector_LVLH;
    attitude.YplusFaceVector_LVLH = param.pwr.YplusFaceVector_LVLH;
    attitude.ZplusFaceVector_LVLH = param.pwr.ZplusFaceVector_LVLH;



    shadowBodyColor = [.4 .4 .4];
    iluminatedBodyColor = [.9 .9 .9];

    iluminatedSolarPanelColor = [1 0 0];
    shadowSolarPanelColor = [0.4 0.1 0.1];
    
    deployable = param.deployable;

    % Get closest data point to user selected time
    sampleTime = param.pwr.time_Epoch(2) - param.pwr.time_Epoch(1);
    time = round(t/sampleTime)*sampleTime;


    [i, ~] = find(param.pwr.time_Epoch == time);

%     if(~isempty(i0))
%         if(playFlag)
%             indexArray = i0:length(param.pwr.time_Epoch);
%         else
%             indexArray = i0;
%         end
%     else
%         indexArray = 1;
%     end
    


    cla(appAxis);
    hold(appAxis,"on")
    for n=1:6
        % ## plot cubesat faces ##
        faceX = reshape(faceVertices_L(:,1,n),[2,2]);
        faceY = reshape(faceVertices_L(:,2,n),[2,2]);
        faceZ = reshape(faceVertices_L(:,3,n),[2,2]);
        if(param.facesMaterial(n) == "Solar Panel")
            if(param.pwr.solarLight(i,n).(1)>0)
                surf(appAxis, faceX, faceY, faceZ,"FaceColor",iluminatedSolarPanelColor,"EdgeColor","k")
            else
                surf(appAxis, faceX, faceY, faceZ,"FaceColor",shadowSolarPanelColor,"EdgeColor","k")
            end
        else
            if(param.pwr.solarLight(i,n).(1)>0)
                surf(appAxis, faceX, faceY, faceZ,"FaceColor",max(iluminatedBodyColor*param.pwr.solarLight(i,n).(1), shadowBodyColor*1.3),"EdgeColor","k")
            else
                surf(appAxis, faceX, faceY, faceZ,"FaceColor",shadowBodyColor,"EdgeColor","k")
            end
        end
        % ## plot cubesat fixed panels ##
        if(~all(fixedDeployableVertices_L(:,:,n)==0,"all"))
            faceXVector = fixedDeployableVertices_L(:,1,n);
            faceYVector = fixedDeployableVertices_L(:,2,n);
            faceZVector = fixedDeployableVertices_L(:,3,n);
            [faceXVector_shifted, faceYVector_shifted, faceZVector_shifted] = pushOutFromFixedPanel(faceXVector, faceYVector, faceZVector,deployable.hinge(n), deployable.flipFixedPanel(n), attitude);
            faceX = reshape(faceXVector,[2,2]);
            faceY = reshape(faceYVector,[2,2]);
            faceZ = reshape(faceZVector,[2,2]);
            surf(appAxis, faceX, faceY, faceZ,"FaceColor",iluminatedSolarPanelColor,"EdgeColor","k")
            % Plot panel cover
            faceX = reshape(faceXVector_shifted,[2,2]);
            faceY = reshape(faceYVector_shifted,[2,2]);
            faceZ = reshape(faceZVector_shifted,[2,2]);
            surf(appAxis, faceX, faceY, faceZ,"FaceColor",iluminatedBodyColor,"EdgeColor","k")
        end

        % ## plot cubesat tracking panels ##
        if(~all( trackingDeployableVertices_L(:,:,n,i)==0,"all") )
            faceXVector = trackingDeployableVertices_L(:,1,n,i);
            faceYVector = trackingDeployableVertices_L(:,2,n,i);
            faceZVector = trackingDeployableVertices_L(:,3,n,i);
            [faceXVector_shifted, faceYVector_shifted, faceZVector_shifted] = pushOutFromTrackingPanel(faceXVector, faceYVector, faceZVector, Sun_L(i,:), "cover");
            faceX = reshape(faceXVector,[2,2]);
            faceY = reshape(faceYVector,[2,2]);
            faceZ = reshape(faceZVector,[2,2]);
            if(param.orb.prop.sunMagnitude(i)>0)
                surf(appAxis, faceX, faceY, faceZ,"FaceColor",iluminatedSolarPanelColor,"EdgeColor","k")
            else
                surf(appAxis, faceX, faceY, faceZ,"FaceColor",shadowSolarPanelColor,"EdgeColor","k")
            end
            % Plot panel cover
            faceX = reshape(faceXVector_shifted,[2,2]);
            faceY = reshape(faceYVector_shifted,[2,2]);
            faceZ = reshape(faceZVector_shifted,[2,2]);
            surf(appAxis, faceX, faceY, faceZ,"FaceColor",shadowBodyColor,"EdgeColor","k")
        end
        
        % ### Shadow if computed sucessfully ###
        if(param.computePwrWithShadowFlag && param.pwr.shadow.flag)
            % ## plot shaded area on body faces
            shadeVertices = param.pwr.shadow.shadowFaceVertices_L{i,n};
            if(~isempty(shadeVertices))
                faceXVector = shadeVertices(:,1);
                faceYVector = shadeVertices(:,2);
                faceZVector = shadeVertices(:,3);
                fill3(appAxis, faceXVector + pushOut(faceXVector),faceYVector + pushOut(faceYVector),faceZVector + pushOut(faceZVector),shadowBodyColor)        
            end
            % ## plot shaded area on panels
            shadeVertices = param.pwr.shadow.shadowPanelVertices_L{i,n};
            if(~isempty(shadeVertices))
                faceXVector = shadeVertices(:,1);
                faceYVector = shadeVertices(:,2);
                faceZVector = shadeVertices(:,3);
                [faceXVector, faceYVector, faceZVector] = pushOutFromTrackingPanel(faceXVector, faceYVector, faceZVector, Sun_L(i,:), "shadow");
                fill3(appAxis, faceXVector,faceYVector,faceZVector,shadowSolarPanelColor)        
            end
        end

        % ## Plot Sun vector
        v = 3*Sun_L(i,:); 
        quiver3(appAxis, [0],[0],[0],v(1),v(2),v(3), 'Color',"#EDB120",'LineWidth',5);
        text(appAxis, v(1),v(2),v(3),'Sun','HorizontalAlignment','left','FontSize',12, 'Color',"#EDB120"); 
    end

    axis(appAxis, "equal");
    view(appAxis, [1,1,1])
    axis(appAxis, "off");
    set(appAxis, 'Zdir', 'reverse')
    
    if(param.pwr.flag)
        numH = floor(time/3600);
        numM = round((time - 3600*numH)/60);
        title(appAxis, strcat("t = ", num2str(numH,2), ':',num2str(numM,'%02d'), "  P = ", num2str(param.pwr.generatedTotalPower(i),2), " W"),'FontWeight','normal', 'BackgroundColor','none');
    end

end

function delta = pushOut(v)
    delta = zeros(length(v),1);
    if(sum(diff(v)) == 0)
        % Check direction to push out panel
        delta = sign(v(1)) * 0.001 * ones(length(v),1);
    end
end

function [faceXVector, faceYVector, faceZVector] = pushOutFromFixedPanel(faceXVector, faceYVector, faceZVector, hingeFace, flipFlag, attitude)
    
    hingeFace = extractAfter(hingeFace, " | ");
    hingeFace = char(hingeFace);
    switch hingeFace(1) % X
        case "X"
            if(hingeFace(2) == '+')
                normalVector = attitude.XplusFaceVector_LVLH;
            else
                normalVector = -attitude.XplusFaceVector_LVLH;
            end
        case "Y"
            if(hingeFace(2) == '+')
                normalVector = attitude.YplusFaceVector_LVLH;
            else
                normalVector = -attitude.YplusFaceVector_LVLH;
            end
        case "Z"
            if(hingeFace(2) == '+')
                normalVector = attitude.ZplusFaceVector_LVLH;
            else
                normalVector = -attitude.ZplusFaceVector_LVLH;
            end
    end

    % Check if the panel is flipped
    if(~flipFlag)
        normalVector = -normalVector;
    end

    faceXVector = faceXVector + 0.02 * normalVector(1);
    faceYVector = faceYVector + 0.02 * normalVector(2);
    faceZVector = faceZVector + 0.02 * normalVector(3);

end

function [faceXVector, faceYVector, faceZVector] = pushOutFromTrackingPanel(faceXVector, faceYVector, faceZVector, sun, purpose)
    % Get the normal direction of the panel
    vertices = [faceXVector - mean(faceXVector), faceYVector - mean(faceYVector), faceZVector - mean(faceZVector)];
    panelNormalVector = cross(vertices(2,:), vertices(1,:));
    panelNormalVector = panelNormalVector/norm(panelNormalVector);
    % Check the direction with the Sun vector
    if(dot(panelNormalVector, sun)<0)
        if(purpose == "shadow")
            panelNormalVector = -panelNormalVector;
        end
    else
        if(purpose == "cover")
            panelNormalVector = -panelNormalVector;
        end
    end

    faceXVector = faceXVector + 0.02 * panelNormalVector(1);
    faceYVector = faceYVector + 0.02 * panelNormalVector(2);
    faceZVector = faceZVector + 0.02 * panelNormalVector(3);

end


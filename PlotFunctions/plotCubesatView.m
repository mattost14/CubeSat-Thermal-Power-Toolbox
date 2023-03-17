function plotCubesatView(param, appAxis)
    facesPoints = getCubeSatPlotPoints(param.numberOfUnits);

    cla(appAxis);

    % Plot CubeSat faces
    faceText = ["X+","X-","Y+","Y-","Z+","Z-"];
    for i=1:length(param.facesMaterial)
        h = surf(appAxis, facesPoints(i).x, facesPoints(i).y, facesPoints(i).z, 'CData',getTexture(param.facesMaterial(i)), 'FaceColor','texturemap', 'edgecolor','#656565', 'LineWidth', 3); 
        textPosition = [mean(facesPoints(i).x,"all"), mean(facesPoints(i).y,"all"),  mean(facesPoints(i).z,"all")];
        textPosition(textPosition~=0) = textPosition(textPosition~=0)*1.4;
        text(appAxis, textPosition(1), textPosition(2), textPosition(3),strcat(faceText(i)),'HorizontalAlignment','left',"FontSize",20)
        hold(appAxis, 'on')
    end

%     Plot Deployable solar panels
    for i=1:length(param.facesMaterial)
        if((param.deployable.type(i) == "Fixed") || (param.deployable.type(i) == "Tracking"))
            hinge = param.deployable.hinge(i);
            size = param.deployable.size(i);
            flipFlag = param.deployable.flipFixedPanel(i);
            [deployablePoints, coverPoints] = getDeployablePlotPoints(facesPoints, hinge, size, flipFlag);
            h = surf(appAxis, deployablePoints.x, deployablePoints.y, deployablePoints.z, 'CData',getTexture("Solar Panel"), 'FaceColor','texturemap', 'edgecolor','#656565', 'LineWidth', 3);
            h = surf(appAxis, coverPoints.x, coverPoints.y, coverPoints.z, 'CData',getTexture("White coating"), 'FaceColor','texturemap', 'edgecolor','#656565', 'LineWidth', 3);
        end
    end


    % Plot Nadir vector
    i = find(faceText == param.attitude.nadirFace);
    v = [mean(facesPoints(i).x,"all"),mean(facesPoints(i).y,"all"),mean(facesPoints(i).z,"all")];
    if(v(v~=0)>0)
        v(v~=0) = v(v~=0) + 1;
    else
        v(v~=0) = v(v~=0) - 1;
    end
    q=quiver3(appAxis, [0],[0],[0],v(1),v(2),v(3), 'Color','g','LineWidth',5);
    text(appAxis, v(1),v(2),v(3),'Nadir','HorizontalAlignment','left','FontSize',12, 'Color','g');

    % Plot Ram vector
    i = find(faceText == param.attitude.ramFace);
    v = [mean(facesPoints(i).x,"all"),mean(facesPoints(i).y,"all"),mean(facesPoints(i).z,"all")];
    if(v(v~=0)>0)
        v(v~=0) = v(v~=0) + 1;
    else
        v(v~=0) = v(v~=0) - 1;
    end
    q=quiver3(appAxis, [0],[0],[0],v(1),v(2),v(3), 'Color','r','LineWidth',5);
    text(appAxis, v(1),v(2),v(3),'Ram','HorizontalAlignment','left','FontSize',12, 'Color','r');


    
    
    axis(appAxis, "equal");
    axis(appAxis, "off");
    set(appAxis, 'Xdir', 'reverse')
   
end
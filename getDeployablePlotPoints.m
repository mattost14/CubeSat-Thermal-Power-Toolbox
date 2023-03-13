function [deployablePoints, coverPoints] = getDeployablePlotPoints(facesPoints, hinge, size, flipFlag)
    faceText = ["X+","X-","Y+","Y-","Z+","Z-"];

    numPanels = str2num(extractBefore(size, "P"));

    panelFace = extractBefore(hinge, " | ");
    panelFaceId = find(faceText==panelFace);

    hingeFace = extractAfter(hinge, " | ");
    hingeFaceId = find(faceText==hingeFace);

    % Get panel face points
    panelFacePoints = facesPoints(panelFaceId);
    hingeFacePoints = facesPoints(hingeFaceId);

    % Remove internal points
    panelFacePoints.x = [panelFacePoints.x(1,1), panelFacePoints.x(1,end); panelFacePoints.x(end,1), panelFacePoints.x(end,end)];
    panelFacePoints.y = [panelFacePoints.y(1,1), panelFacePoints.y(1,end); panelFacePoints.y(end,1), panelFacePoints.y(end,end)];
    panelFacePoints.z = [panelFacePoints.z(1,1), panelFacePoints.z(1,end); panelFacePoints.z(end,1), panelFacePoints.z(end,end)];

    
    % Define rotation axis
    listOfFaces = ["X", "Y", "Z"];
    panelFace = char(panelFace); hingeFace = char(hingeFace);
    filterCondition = find((listOfFaces~=panelFace(1)) & (listOfFaces~=hingeFace(1)));
    hingeAxis = listOfFaces(filterCondition);

    % Get panel length along hinge face axis
    switch hingeFace(1)
        case "X"
            panelLength = max(panelFacePoints.x(:)) - min(panelFacePoints.x(:));
        case "Y"
            panelLength = max(panelFacePoints.y(:)) - min(panelFacePoints.y(:));
        case "Z"
            panelLength = max(panelFacePoints.z(:)) - min(panelFacePoints.z(:));
    end 
    panelLength = numPanels * panelLength;


   
    originalPoints = [panelFacePoints.x(:),panelFacePoints.y(:), panelFacePoints.z(:)];
    faceCG =  mean(originalPoints);
    panelCG = (1 + panelLength/(2*norm(faceCG)) ) * faceCG;
    hingeFaceCG = mean([hingeFacePoints.x(:),hingeFacePoints.y(:), hingeFacePoints.z(:)]);

    % Translate points to origin
    originalPointsAtOrigin(:,1) = originalPoints(:,1) - faceCG(1);
    originalPointsAtOrigin(:,2) = originalPoints(:,2) - faceCG(2);
    originalPointsAtOrigin(:,3) = originalPoints(:,3) - faceCG(3);
    

    % Defining the rotation matrix
    switch hingeAxis
        case "X"
            rotM = eul2rotm([deg2rad(90),0,0],"XYZ");
        case "Y"
            rotM = eul2rotm([0,deg2rad(90),0],"XYZ");
        case "Z"
            rotM = eul2rotm([0,0,deg2rad(90)],"XYZ");
    end

    for i=1:length(originalPointsAtOrigin)
        p = originalPointsAtOrigin(i,:);
        % Perform rotation
        rotP = rotM * p';
        % Stretching the pane
        switch panelFace(1)
            case "X"
                rotP(1) = rotP(1) * numPanels;
            case "Y"
                rotP(2) = rotP(2) * numPanels;
            case "Z"
                rotP(3) = rotP(3) * numPanels;
        end      
        % Translate 
        rotP = rotP + panelCG' + hingeFaceCG';
        rotatedPoints(i,:) = rotP;

        % Translate more to create the cover panel
        if(flipFlag)
            rotPcover = rotP + 0.01*(hingeFaceCG'/norm(hingeFaceCG));
        else
            rotPcover = rotP - 0.01*(hingeFaceCG'/norm(hingeFaceCG));
        end

        % Store points
        rotatedPoints(i,:) = rotP;
        rotatedPointsCover(i,:) = rotPcover;
    end

    deployablePoints.x = reshape(rotatedPoints(:,1),2,2);
    deployablePoints.y = reshape(rotatedPoints(:,2),2,2);
    deployablePoints.z = reshape(rotatedPoints(:,3),2,2);

    % Cover points
    coverPoints.x = reshape(rotatedPointsCover(:,1),2,2);
    coverPoints.y = reshape(rotatedPointsCover(:,2),2,2);
    coverPoints.z = reshape(rotatedPointsCover(:,3),2,2);

end
function shadow = computeShadowFactors(param, faceSolarLight)
%COMPUTESHADOWFACTORS computes the the proportion of each face area that is
%covered by the deployable panel's shadow.

    deployable = param.deployable;
    sunMagnitude = param.orb.prop.sunMagnitude;
    
    % Get vertices points in the LVLH frame
    [faceVertices_L, fixedDeployableVertices_L, trackingDeployableVertices_L, Sun_L] = getVerticesPoints(param);
        
    % Initialize arrays
    shadowFaceAreaRatio = zeros(length(Sun_L(:,1)),6);
    shadowPanelAreaRatio = zeros(length(Sun_L(:,1)),6);
    shadowFaceVertices_L = cell(length(Sun_L(:,1)),6);
    shadowPanelVertices_L = cell(length(Sun_L(:,1)),6);
    
    for i=1:length(Sun_L(:,1))
        sun = Sun_L(i,:);
    
        if(sunMagnitude(i))
            % Loop through faces that has either a fixed or tracking panel
            for n=1:6
                facePoints = faceVertices_L(:,:,n);
                if(deployable.type(n) == "Fixed") 
                    panelPoints = fixedDeployableVertices_L(:,:,n);  
                elseif(deployable.type(n) == "Tracking")
                    panelPoints = trackingDeployableVertices_L(:,:,n,i);
                else
                    continue; % skip
                end
    
                %% Shadow from Deployable Panel to Faces
    
                %%% ## Check shadow for the face underneath the panel ##
                if(faceSolarLight(i,n).(1))
                    [shadowFaceAreaRatio(i,n), shadowFaceVertices_L{i,n}] =  computeShadowProjectionPanel2Face(facePoints, panelPoints, sun);
                end
    
                %%% ## Check shadow for the face neighboor to the panel ##
                % Get hinge face vertices
                faceIDText = ["X+","X-","Y+","Y-","Z+","Z-"];
                hingeFace = extractAfter(deployable.hinge(n), " | ");
                hingeFaceIdx = find(faceIDText == hingeFace);
                hingeFacePoints = faceVertices_L(:,:,hingeFaceIdx);
                if(faceSolarLight(i,hingeFaceIdx).(1) && faceSolarLight(i,n).(1))
                    [shadowFaceAreaRatio(i,hingeFaceIdx), shadowFaceVertices_L{i,hingeFaceIdx}] =  computeShadowProjectionPanel2Face(hingeFacePoints, panelPoints, sun);
                end
    
                %% Shadow from Body to Panel
                % If the face that holds the deployable is not illuminated
                % during day light, then check for shadows projected over the
                % panel
                if(~faceSolarLight(i,n).(1))
                    [shadowPanelAreaRatio(i,n), shadowPanelVertices_L{i,n}] = computeShadowProjectionBody2Panel(faceVertices_L, panelPoints, sun, faceSolarLight(i,:));
                end
            end
        end
    end

    % Return shadow analysis results
    shadow.faceVertices_L = faceVertices_L;
    shadow.fixedDeployableVertices_L = fixedDeployableVertices_L;
    shadow.trackingDeployableVertices_L = trackingDeployableVertices_L;
    shadow.Sun_L = Sun_L;
    shadow.shadowFaceVertices_L = shadowFaceVertices_L;
    shadow.shadowPanelVertices_L = shadowPanelVertices_L;
    shadow.shadowFaceAreaRatio = shadowFaceAreaRatio;
    shadow.shadowPanelAreaRatio = shadowPanelAreaRatio;

%     plotAnimatedCubeSat(param, shadowFaceVertices_L, shadowPanelVertices_L)

end

function [shadowPanelAreaRatio, shadowPanelVertices_L] = computeShadowProjectionBody2Panel(bodyPoints, panelPoints, sun, faceSolarLight)
    % Computes the shadow that the body casts over the deployable panel

    % Check for iluminated faces
    iluminatedFaces = 0;
    bodyVertices = [];
    for n=1:6
        if(faceSolarLight(1,n).(1))
            iluminatedFaces = iluminatedFaces + 1;
            bodyVertices = [bodyVertices; bodyPoints(:,:,n)];
        end
    end

    % Remove common vertex
    bodyVertices_L = unique(bodyVertices,'rows');

    % Change the Frame to make x-axis along the sun vector (let's call
    % Sun-View Frame)
    x_S = -sun';
    x_S = x_S/norm(x_S);
    z_S = cross([0;0;1], x_S);
    if(norm(z_S)==0)
        z_S = cross([0;1;0], x_S);
    end
    z_S = z_S/norm(z_S);
    y_S = cross(z_S,x_S);
    y_S = y_S/norm(y_S);
    A_L2S = [x_S'; y_S'; z_S'];

    % Chante the bodyVertices from LVLH to Sun-View Frame
    for vertex=1:length(bodyVertices_L)
        bodyVertices_S = (A_L2S * bodyVertices_L(vertex,:)')'; 
        bodyVertices_2D(vertex,:) = [bodyVertices_S(1,2), bodyVertices_S(1,3)];
    end
    
    % Remove inner most vertices
    distToOrigin = sqrt(bodyVertices_2D(:,1).^2 + bodyVertices_2D(:,2).^2);
    [~, idx] = sort(distToOrigin,"descend");   % sorting by distance 
    bodyVertices_2D = bodyVertices_2D(idx,:); % sorting the given points
    switch iluminatedFaces
        case 1 || 2
            bodyVertices_2D = bodyVertices_2D(1:4,:);
        case 3
            bodyVertices_2D = bodyVertices_2D(1:6,:);
    end

    % Chante the panelVertices from LVLH to Sun-View Frame
    for vertex=1:length(panelPoints)
        panelPoints_S = (A_L2S * panelPoints(vertex,:)')'; 
        panelVertices_2D(vertex,:) = [panelPoints_S(1,2), panelPoints_S(1,3)];
    end

    % Sort the points to correctly close the polygon
    % Face polygon
    c = mean(bodyVertices_2D); % mean/ central point
    d = bodyVertices_2D - [c(1)*ones(length(bodyVertices_2D),1), c(2)*ones(length(bodyVertices_2D),1)]; % vectors connecting the central point and the given points 
    th = atan2d(d(:,2),d(:,1)); % angle above x axis
    [~, idx] = sort(th);   % sorting the angles 
    bodyVertices_2D = bodyVertices_2D(idx,:); % sorting the given points
    % Panel shadow polygon
    c = mean(panelVertices_2D); % mean/ central point
    d = panelVertices_2D-[c(1)*ones(length(panelVertices_2D),1), c(2)*ones(length(panelVertices_2D),1)]; % vectors connecting the central point and the given points 
    th = atan2d(d(:,2),d(:,1)); % angle above x axis
    [~, idx] = sort(th);   % sorting the angles 
    panelVertices_2D = panelVertices_2D(idx,:); % sorting the given points

    bodyProjectionPolygon = polyshape(bodyVertices_2D(:,1), bodyVertices_2D(:,2));
    panelPolygon = polyshape(panelVertices_2D(:,1), panelVertices_2D(:,2));

    shadowPolygon = intersect(bodyProjectionPolygon,panelPolygon);

    % Compute the ratio of the face area covered by the shadow
    shadowPanelAreaRatio = area(shadowPolygon)/area(panelPolygon);

    % Convert from Sun-View Frame to LVLH Frame
    shadowPanelVertices_L = [];
    for vertex=1:length(shadowPolygon.Vertices)
        % Shadow vertices in the YZ plane of the Sun-View Frame
        p_S = [0; shadowPolygon.Vertices(vertex,:)'];
        % Sun-View Frame -> LVLH
        A_S2L = A_L2S'; 
        p_L = (A_S2L * p_S);

        % Line paremeters in the LVLH frame
        u = sun;                    % director vector of the parametric line
        N = p_L';                    % point belonging to the line
        
        % Panel Plane parameters
        CM = mean(panelPoints);                                 % Center of mass of the plane
        v = cross(panelPoints(2,:)-CM, panelPoints(1,:)-CM);    % vector normal to the plane
        v = v/norm(v);
        M = panelPoints(1,:);                                   % point belonging to the plane
        
        % Compute intersection
        [I ,rc] = line_plane_intersection(u, N, v, M);

        if(rc==1)
            shadowPanelVertices_L(vertex,:) = I;
        end

    end

end

function [shadowFaceAreaRatio, shadowFaceVertices_L] = computeShadowProjectionPanel2Face(facePoints, panelPoints, sun)

    % Compute the intersection between the line that passes over each 
    % panel vertex along the sun vector and the plane that contains the
    % face
    for vertex=1:4
        % Line paremeters
        u = sun;                    % director vector of the parametric line
        N = panelPoints(vertex,:);  % point belonging to the line
        % Plane parameters
        CM = mean(facePoints);                              % Center of mass of the plane
        v = cross(facePoints(2,:)-CM, facePoints(1,:)-CM);  % vector normal to the plane
        M = facePoints(1,:);                                % point belonging to the plane
        % Compute intersection
        [I(vertex,:),rc] = line_plane_intersection(u, N, v, M);
    end



    if(rc == 1) % Point intersection (it forms a polygon)
        % Remove the third dimension
        if sum(abs(diff(facePoints(:,1)))) == 0
            excludedCol = 1;
            excludedArray = facePoints(:,1);
            facePoints2D = facePoints(:,2:3);
            I_2D = I(:,2:3);
        elseif sum(abs(diff(facePoints(:,2)))) == 0
            excludedCol = 2;
            excludedArray = facePoints(:,2);
            facePoints2D = [facePoints(:,1), facePoints(:,3)];
            I_2D = [I(:,1), I(:,3)];
        else
            excludedCol = 3;
            excludedArray = facePoints(:,3);
            facePoints2D = [facePoints(:,1:2)];
            I_2D = [I(:,1:2)];
        end
       
        % Sort the points to correctly close the polygon
        % Face polygon
        c = mean(facePoints2D); % mean/ central point
        d = facePoints2D-[c;c;c;c] ; % vectors connecting the central point and the given points 
        th = atan2d(d(:,2),d(:,1)); % angle above x axis
        [~, idx] = sort(th);   % sorting the angles 
        facePoints2D = facePoints2D(idx,:); % sorting the given points
        % Panel shadow polygon
        c = mean(I_2D); % mean/ central point
        d = I_2D-[c;c;c;c] ; % vectors connecting the central point and the given points 
        th = atan2d(d(:,2),d(:,1)); % angle above x axis
        [~, idx] = sort(th);   % sorting the angles 
        I_2D = I_2D(idx,:); % sorting the given points

        facePolygon = polyshape(facePoints2D(:,1), facePoints2D(:,2));
        panelProjectionPolygon = polyshape(I_2D(:,1), I_2D(:,2));

        shadowPolygon = intersect(facePolygon,panelProjectionPolygon);

        % Compute the ratio of the face area covered by the shadow
        shadowFaceAreaRatio = area(shadowPolygon)/area(facePolygon);

        % Save shadow polygon vertices in 3D to be later plotted
        switch excludedCol
            case 1
                shadowFaceVertices_L = [excludedArray(1)*ones(length(shadowPolygon.Vertices),1), shadowPolygon.Vertices];
            case 2
                shadowFaceVertices_L = [shadowPolygon.Vertices(:,1), excludedArray(1)*ones(length(shadowPolygon.Vertices),1), shadowPolygon.Vertices(:,2)];
            case 3
                shadowFaceVertices_L = [shadowPolygon.Vertices, excludedArray(1)*ones(length(shadowPolygon.Vertices),1)];
        end
    end
end




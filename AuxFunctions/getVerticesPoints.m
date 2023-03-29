function [faceVertices_L, fixedDeployableVertices_L, trackingDeployableVertices_L, Sun_L] = getVerticesPoints(param)
    % Get the vertices points in the LVLH Frame

    [XplusFaceVector_LVLH, YplusFaceVector_LVLH, ZplusFaceVector_LVLH] = getCubeSatAttitude(param);

    % Rotation Matrix from Body to LVLH
    A_B2L = [XplusFaceVector_LVLH, YplusFaceVector_LVLH, ZplusFaceVector_LVLH];
  
    %% Faces Vertices
    facesPoints = getCubeSatPlotPoints(param.numberOfUnits);
    % Remove internal points
    faceVertices_L=zeros(4,3,6);
    for n=1:6
        panelFacePoints = facesPoints(n);
        x = [panelFacePoints.x(1,1), panelFacePoints.x(1,end); panelFacePoints.x(end,1), panelFacePoints.x(end,end)];
        y = [panelFacePoints.y(1,1), panelFacePoints.y(1,end); panelFacePoints.y(end,1), panelFacePoints.y(end,end)];
        z = [panelFacePoints.z(1,1), panelFacePoints.z(1,end); panelFacePoints.z(end,1), panelFacePoints.z(end,end)];
        faceVertices_B = [x(:), y(:),z(:)];
        for i=1:4
            faceVertices_L(i,:,n) = (A_B2L * faceVertices_B(i,:)')';
        end
    end



    %% Fixed Deployable Panel Vertices
    fixedDeployableVertices_L=zeros(4,3,6);
    for n=1:6
        if(param.deployable.type(n)=="Fixed")
            hinge = param.deployable.hinge(n);
            size = param.deployable.size(n);
            flipFlag = param.deployable.flipFixedPanel(n);
            [deployablePoints, ~] = getDeployablePlotPoints(facesPoints, hinge, size, flipFlag);
            x = [deployablePoints.x(1,1), deployablePoints.x(1,end); deployablePoints.x(end,1), deployablePoints.x(end,end)];
            y = [deployablePoints.y(1,1), deployablePoints.y(1,end); deployablePoints.y(end,1), deployablePoints.y(end,end)];
            z = [deployablePoints.z(1,1), deployablePoints.z(1,end); deployablePoints.z(end,1), deployablePoints.z(end,end)];
            fixedDeployableVertices_B = [x(:), y(:),z(:)];
            for i=1:4
                fixedDeployableVertices_L(i,:,n) = (A_B2L * fixedDeployableVertices_B(i,:)')';
            end
        end
    end



    %% Tracking Deployable Panel Vertices

    % Get orbit propagation position and velocity in the inertial frame
    r_I = param.orb.prop.r_I;
    v_I = param.orb.prop.v_I;
    time_Epoch = param.orb.prop.time_Epoch;
    
    % Get sun lightining
    sun_I = param.orb.prop.Sun_I;

    trackingDeployableVertices_L=zeros(4,3,6,length(time_Epoch)); 
    % Get tracking panel initial position (t=0)
    atLeastOneTrackingPanel = 0;
    for n=1:6
        if(param.deployable.type(n)=="Tracking")
            atLeastOneTrackingPanel = 1;
            hinge = param.deployable.hinge(n);
            size = param.deployable.size(n);
            flipFlag = param.deployable.flipFixedPanel(n);
            [deployablePoints, ~] = getDeployablePlotPoints(facesPoints, hinge, size, flipFlag);
            x = [deployablePoints.x(1,1), deployablePoints.x(1,end); deployablePoints.x(end,1), deployablePoints.x(end,end)];
            y = [deployablePoints.y(1,1), deployablePoints.y(1,end); deployablePoints.y(end,1), deployablePoints.y(end,end)];
            z = [deployablePoints.z(1,1), deployablePoints.z(1,end); deployablePoints.z(end,1), deployablePoints.z(end,end)];
            P0_B = [x(:), y(:),z(:)];
            for i=1:4
                trackingDeployableVertices_L(i,:,n,1) = (A_B2L * P0_B(i,:)')';
                trackingDeployableVertices_L_default(i,:,n) = trackingDeployableVertices_L(i,:,n,1);
            end
            
        end
    end


    % Compute the Sun vector at the LVLH frame
    Sun_L = zeros(length(r_I),3);
    for i=1:length(r_I)
        % LVLH unit vectors 
        % (z: Nadir, 
        % y: opposite direction of orbit angular momentum, 
        % x: completes right-handed frame ~ velocity direction)
        r = r_I(i,:)'; v = v_I(i,:)'; 
        z_L2I = -r/norm(r);   % Unit vector z
        h_I = cross(r, v);    % Angular momentum 
        y_L2I = -h_I/norm(h_I);             % Unit vector y
        x_L2I = cross(y_L2I, z_L2I);        % Unit vector x
        % Rotation matrix from Inertial Frame to LVLH frame
        A_I2L = [x_L2I'; y_L2I'; z_L2I'];
    
        % Compute sun unit vector at LVLH frame
        sun = sun_I(i,:)';
        sun_L = A_I2L * sun;
        sun_L = sun_L' / norm(sun_L);
        Sun_L(i,:) = sun_L;

         if(atLeastOneTrackingPanel)
            for n=1:6
                if(param.deployable.type(n)=="Tracking")
                    
                    % Initial Normal Vector
                    hingeFace = extractAfter(param.deployable.hinge(n), " | ");
                    hingeFace = char(hingeFace);
                    faceSignal = 1;
                    if(hingeFace(2)=="-")
                        faceSignal = -1;
                    end
                    switch hingeFace(1)
                        case "X"
                            normalVector = faceSignal * XplusFaceVector_LVLH;
                        case "Y"
                            normalVector = faceSignal * YplusFaceVector_LVLH;
                        case "Z"
                            normalVector = faceSignal * ZplusFaceVector_LVLH;
                    end

                    % Panel points in the default position
                    P0_L = trackingDeployableVertices_L_default(:, :, n);

                    % Define the axis of rotaiton in the LVLH frame
                    if(n==1||n==2) % X face
                        faceNormal = XplusFaceVector_LVLH;
                    elseif(n==3||n==4)
                        faceNormal = YplusFaceVector_LVLH;
                    else
                        faceNormal = ZplusFaceVector_LVLH;
                    end

                    if(sum(abs(abs(faceNormal) - [1;0;0])) == 0)
                        lvlh_rotAxis = faceNormal(1); % x-axis
                    elseif(sum(abs(abs(faceNormal) - [0;1;0])) == 0)
                        lvlh_rotAxis = 2*faceNormal(2); % y-axis
                    else
                        lvlh_rotAxis = 3*faceNormal(3); % z-axis
                    end


                    % Define the rotation matrix
                    if(abs(lvlh_rotAxis) == 1) % x-axis rotation in the LVLH frame
%                         rotationAngle = (lvlh_rotAxis/abs(lvlh_rotAxis)) * calculateAngleBetweenVectors([0;sun_L(2);sun_L(3)], normalVector);
                        rotationAngle =  (lvlh_rotAxis/abs(lvlh_rotAxis)) * vecangle360([0;sun_L(2);sun_L(3)], normalVector, [(lvlh_rotAxis/abs(lvlh_rotAxis)); 0; 0]);
                        rotM = [1        0               0;
                                0 cosd(rotationAngle) -sind(rotationAngle);
                                0 sind(rotationAngle) cosd(rotationAngle)];
                    elseif(abs(lvlh_rotAxis) == 2) % y-axis rotation in the LVLH frame
%                         rotationAngle = (lvlh_rotAxis/abs(lvlh_rotAxis)) * calculateAngleBetweenVectors([sun_L(1); 0; sun_L(3)], normalVector);
                        rotationAngle =  (lvlh_rotAxis/abs(lvlh_rotAxis)) * vecangle360([sun_L(1); 0; sun_L(3)], normalVector, [0; (lvlh_rotAxis/abs(lvlh_rotAxis)); 0]);
                        rotM = [cosd(rotationAngle) 0  sind(rotationAngle); 
                                0                   1           0; 
                                -sind(rotationAngle) 0 cosd(rotationAngle)];
                    else % z-axis rotation in the LVLH frame
%                         rotationAngle = (lvlh_rotAxis/abs(lvlh_rotAxis)) * calculateAngleBetweenVectors([sun_L(1); sun_L(2); 0], normalVector);
                        rotationAngle =  (lvlh_rotAxis/abs(lvlh_rotAxis)) * vecangle360([sun_L(1); sun_L(2); 0], normalVector, [0; 0; (lvlh_rotAxis/abs(lvlh_rotAxis))]);
                        rotM = [cosd(rotationAngle) -sind(rotationAngle) 0; 
                                sind(rotationAngle) cosd(rotationAngle) 0; 
                                0                   0                   1];
                    end

                    % Translate to Center of Mass
                    CM = mean(P0_L);
                    P0_L_translated = P0_L - [CM; CM; CM; CM];
                    % Rotate the points using the rotM
                    for k=1:4
                        P0_L_translated_rotated(k,:) = (rotM' * P0_L_translated(k,:)')';
                    end
                    % Translate back from Center of Mass
                    trackingDeployableVertices_L(:, :, n, i) = P0_L_translated_rotated + [CM; CM; CM; CM];                      
                end
            end
        end
    end


end

function a = vecangle360(v1,v2,n)
    x = cross(v1,v2);
    c = sign(dot(x,n)) * norm(x);
    a = atan2d(c,dot(v1,v2));
end
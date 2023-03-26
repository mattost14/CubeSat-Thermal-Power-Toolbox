function solarLight = computeFaceSolarLight(param)
%COMPUTEFACESOLARLIGHT Summary of this function goes here
%   Detailed explanation goes here

% Get attitude 
XplusFaceVector_LVLH = param.attitude.XplusFaceVector_LVLH;
YplusFaceVector_LVLH = param.attitude.YplusFaceVector_LVLH;
ZplusFaceVector_LVLH = param.attitude.ZplusFaceVector_LVLH;

% Get orbit propagation position and velocity in the inertial frame
r_I = param.orb.prop.r_I;
v_I = param.orb.prop.v_I;
time_Epoch = param.orb.prop.time_Epoch;
sampleTime = time_Epoch(2) - time_Epoch(1);

% Get sun lightining
sun_I = param.orb.prop.Sun_I;
sunMagnitude = param.orb.prop.sunMagnitude;

% function that return angle between two vectors in degrees
calculateAngleBetweenVectors = @(u,v) acosd(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1));

% Compute sun lightining at each face
sunMagnitudeXfaces = zeros(length(time_Epoch),1);
sunMagnitudeYfaces = zeros(length(time_Epoch),1);
sunMagnitudeZfaces = zeros(length(time_Epoch),1);
for i=1:length(time_Epoch)
    if(sunMagnitude(i))
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
        sun_L = sun_L / norm(sun_L);
    
        % Compute sun angles at each face
        angle_faceXNormal_Sun = calculateAngleBetweenVectors(sun_L,XplusFaceVector_LVLH);
        angle_faceYNormal_Sun = calculateAngleBetweenVectors(sun_L,YplusFaceVector_LVLH);
        angle_faceZNormal_Sun = calculateAngleBetweenVectors(sun_L,ZplusFaceVector_LVLH);
        
        sunMagnitudeXfaces(i) = cosd(angle_faceXNormal_Sun);
        sunMagnitudeYfaces(i) = cosd(angle_faceYNormal_Sun);
        sunMagnitudeZfaces(i) = cosd(angle_faceZNormal_Sun); 
    end
end

solarLight = [time_Epoch, ...
              max(sunMagnitudeXfaces,0),...
              max(-sunMagnitudeXfaces,0),...
              max(sunMagnitudeYfaces,0),...
              max(-sunMagnitudeYfaces,0),...
              max(sunMagnitudeZfaces,0),...
              max(-sunMagnitudeZfaces,0),...
              sunMagnitude];
          
solarLight = array2table(solarLight);
solarLight.Properties.VariableNames = {'time','x_plus_magnitude','x_minus_magnitude','y_plus_magnitude','y_minus_magnitude', 'z_plus_magnitude', 'z_minus_magnitude','magnitude'};
solarLight.time = seconds(solarLight.time);
solarLight = table2timetable(solarLight);
end


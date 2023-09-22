function [propData, flag] = propagateMoonOrbit(orb)


%% Define Mission Parameters

mission.StartDate = datetime(year(orb.startDate),month(orb.startDate),day(orb.startDate),orb.startHour,orb.startMinute,0, "TimeZone","UTC");
mission.Duration  = hours(orb.durationHours);


%% Configure the Model

mission.mdl = "LunarOrbitPropagator";
% open_system(mission.mdl);

% Use a SimulationInput object to configure the model for our mission.
mission.sim.in = Simulink.SimulationInput(mission.mdl);

% Define the path to the Orbit Propagator block in the model.
sat.blk = mission.mdl + "/Orbit Propagator";

% Load Moon properties into the base workspace.
moon.F = 0.0012;  % Moon ellipticity (flattening) (Ref 1)
moon.R_eq = 1737400; % [m] Lunar radius in meters (Ref 1)
moon.ReferenceEllipsoid = referenceEllipsoid("moon","meter"); % Moon reference ellipsoid
% moon.Data = matfile("lunarGeographicalData.mat"); % Load moon geographical data

calculateAngleBetweenVectors = @(u,v) acosd(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1)); % function that return angle between two vectors in degrees
% calculateSignedAngleBetweenVectors = @(w,v) atan2d(w(2)*v(1)-w(1)*v(2), w(1)*v(1)+w(2)*v(2));

%% Adjust the orbital parameter from the user frame to the ICRF frame used by the Orbit Propagor Block
% The orbit propagor block considers that the inclination angle is relative to the ICRF X-Y plane. 
% The ICRF X-Y axis is normal to Earth's north pole. 
% The axial tilt of Earth relative to the ecliptic is ~23.44 degrees, while the axial tilt of the Moon is ~5.145 degrees. 
% Therefore, the axial tilt of the Moon relative to the ICRF X-Y plane varies between approximately 23.44±5.145 degrees.
% That is very confusing for the user, so we are going to check the
% relative tilt angle of the of the Moon to the Earth to adjust the inclination angle and make it relative to the X-Y plane on the ME frame.
    moonAngles = moonLibration(juliandate(mission.StartDate), '421');
    % Rotation from Fixed frame (Principal Axis) to Inertial
    rotm_PA2I = euler2rotmZXZ(-moonAngles(1), -moonAngles(2), -moonAngles(3));
    % Rotation from  Inertial to Fixed frame (Principal Axis)
    rotm_I2PA = rotm_PA2I';
    % Rotation from PA to mean-Earth/mean-rotation (MER) frame
    rotm_PA2ME = euler2rotmZYX(-deg2rad(dms2degrees([0 0 0.2785])), -deg2rad(dms2degrees([0 1 18.6944])), -deg2rad(dms2degrees([0 1 7.8526])));
    % Rotation from Inertial to MER
    rotm_I2ME = rotm_PA2ME * rotm_I2PA;
    rotm_ME2I = rotm_I2ME';

    % The user define the inclination and RAAN relative to following frame
    % z-axis: same z-axis od the ME frame (rotation axis)
    % x-axis: same direction of the inertial x-axis define by ICRF
    % Let call this frame keplerian frame (K-frame)

    % rotation about ME(z-axis) to align with its x-axis with the ICRF (x-axis)
    xI_ME = rotm_I2ME * [1;0;0];
    % Define the x-axis of the K-frame
    xK = [xI_ME(1); xI_ME(2); 0]; 
    xK = xK/norm(xK);
    zK = [0;0;1];
    yK = cross(zK, xK);
    % Rotation Matrix from ME to K-frame
    A_ME2K = [xK'; yK'; zK'];
    A_K2ME = A_ME2K';
    % Rotation Matrix from K-frame to ICRF
    A_K2I = rotm_ME2I * A_K2ME;
    
    % Get position and velocity in the K-Frame
    [r0_K,v0_K] = orb2rv(orb.semiMajorAxis, orb.eccentricity, deg2rad(orb.inclination), deg2rad(orb.RAAN), deg2rad(orb.argumentOfPeriapsis), deg2rad(orb.trueAnomaly));
    
    % Transform position and velocity from K-Frame to ICRF frame using
    % rotation matrix
    r0_I = A_K2I * r0_K;  % unit: km
    v0_I = A_K2I * v0_K;  % unit: km/s

    % Get orbital parameters in the ICRF frame
    [a,ecc,incl,raan,argPeriapsis, trueAnomaly, ~, ~, ~, ~] = rv2orb(r0_I,v0_I);

    % Enter the adjusted ICRF orbital parameters
    sat.SemiMajorAxis  = a*1e3;         % [m]
    sat.Eccentricity   = ecc;
    sat.Inclination    = rad2deg(incl);    % [deg]
    sat.RAAN           = rad2deg(raan);    % [deg]
    sat.ArgOfPeriapsis = rad2deg(argPeriapsis);    % [deg]
    sat.TrueAnomaly    = rad2deg(trueAnomaly);      % [deg]
%%

% Set CSM initial conditions. To assign the Keplerian orbital element set defined in the previous section, use setBlockParameter.
mission.sim.in = mission.sim.in.setBlockParameter(...
    sat.blk, "startDate", string(juliandate(mission.StartDate)),...
    sat.blk, "stateFormatNum", "Orbital elements",...
    sat.blk, "orbitType", "Keplerian",...
    sat.blk, "semiMajorAxis", string(sat.SemiMajorAxis),...
    sat.blk, "eccentricity", string(sat.Eccentricity),...
    sat.blk, "inclination", string(sat.Inclination ),...
    sat.blk, "raan", string(sat.RAAN),...
    sat.blk, "argPeriapsis", string(sat.ArgOfPeriapsis),...
    sat.blk, "trueAnomaly", string(sat.TrueAnomaly));

% Set the position and velocity output ports of the block to use the Moon-fixed frame.  The fixed-frame for the Moon is the Mean Earth/Pole Axis (ME) reference system.
mission.sim.in = mission.sim.in.setBlockParameter(...
    sat.blk, "centralBody", "Moon",...
    sat.blk, "outportFrame", "ICRF"); % This parameter is ignored if propagator is "Kepler", but it keeps always in the Inertial frame (ICRF)


% Configure the propagator.
mission.sim.in = mission.sim.in.setBlockParameter(...
    sat.blk, "propagator", orb.orbitPropagator,...
    sat.blk, "gravityModel", "Spherical Harmonics",...
    sat.blk, "moonSH", "LP-100K",... % moon spherical harmonic potential model
    sat.blk, "shDegree", "100",... % Spherical harmonic model degree and order
    sat.blk, "useMoonLib", "off");


% Using constant time step 
mission.sim.in = mission.sim.in.setModelParameter(...
FixedStep="60",...
StopTime=string(seconds(mission.Duration)));

% if(orb.orbitPropagator == "Numerical (high precision)")
%     % Apply model-level solver settings using setModelParameter.  
%     % For best performance and accuracy when using a numerical propagator, use a variable-step solver.
%         mission.sim.in = mission.sim.in.setModelParameter(...
%         FixedStep="60",...
%         StopTime=string(seconds(mission.Duration)));
% %     mission.sim.in = mission.sim.in.setModelParameter(...
% %         SolverType="Variable-step",...
% %         SolverName="VariableStepAuto",...
% %         RelTol="1e-6",...
% %         AbsTol="1e-7",...
% %         StopTime=string(seconds(mission.Duration)));
% else
%     mission.sim.in = mission.sim.in.setModelParameter(...
%     "FixedStep","60",...
%     StopTime=string(seconds(mission.Duration)));
% end


% Save model output port data as a dataset of timetable objects.
mission.sim.in = mission.sim.in.setModelParameter(...
    SaveOutput="on",...
    OutputSaveName="yout",...
    SaveFormat="Dataset",...
    DatasetSignalFormat="timetable");

%% RUN the model and collect ephemerides
try
%     options = simset('SrcWorkspace','current', 'DstWorkspace', 'current');
    mission.sim.out = sim(mission.sim.in);
    
    % Extract the position and velocity data from the model output data structure.
    % Expected ref frame: Inertial-Frame always !!!! 
    % ### Inertial-frame ###
    sat.r_I = mission.sim.out.yout{1}.Values;
    sat.v_I = mission.sim.out.yout{2}.Values;
    r_I = sat.r_I.Data; % m
    v_I = sat.v_I.Data; % m/s
    

    % Time vectors (UTC, Epoch, Julian)
    time_Epoch =  seconds(sat.r_I.Time);
    sat.r_I.Properties.StartTime = mission.StartDate;
    sat.v_I.Properties.StartTime = mission.StartDate;
    time_UTC = sat.r_I.Time;
    time_Julian = juliandate(time_UTC);

    % Sun position x,y,z at Moon-Centered (Inertial) frame
    sun_I = planetEphemeris(time_Julian,'Moon','Sun'); % units: km

    % Earth position x,y,z at Moon-Centered (Inertial) frame
    earth_I = planetEphemeris(time_Julian,'Moon','Earth'); % units: km


    % ### Convert Position and Velocity from Inertial to Moon Fixed Frame
    % ###
    moonAngles = moonLibration(time_Julian, '421');
    for incl=1:length(time_Julian)
        % Rotation from Fixed frame (Principal Axis) to Inertial
        rotm_PA2I = euler2rotmZXZ(-moonAngles(incl,1), -moonAngles(incl,2), -moonAngles(incl,3));
        % Rotation from  Inertial to Fixed frame (Principal Axis)
        rotm_I2PA = rotm_PA2I';
        % Rotation from PA to mean-Earth/mean-rotation (MER) frame
        rotm_PA2ME = euler2rotmZYX(-deg2rad(dms2degrees([0 0 0.2785])), -deg2rad(dms2degrees([0 1 18.6944])), -deg2rad(dms2degrees([0 1 7.8526])));
        % Rotation from Inertial to MER
        rotm_I2ME = rotm_PA2ME * rotm_I2PA;
%         rotm_ME2I = rotm_I2ME';
        % Convert position, velocity from inertial to fixed frame
        r_f = rotm_I2ME * r_I(incl,:)';
        v_f = rotm_I2ME * v_I(incl,:)';
        r_ff(incl,:) = r_f';
        v_ff(incl,:) = v_f';
        % Sub-solar coordinates
        sun_vec = sun_I(incl, :) / norm(sun_I(incl, :));
        sun_ff(incl,:) = rotm_I2ME * sun_vec';
    end
    % ### End of Frame Conversion ####


    % Compute latitude, longitude, and altitude using lunar equatorial radius and flattening. Values are displayed in degrees and meters.
    lla = ecef2lla(r_ff, moon.F, moon.R_eq);
    sat.LLA = timetable(time_UTC, ...
        lla(:,1), lla(:,2), lla(:,3), ...
        VariableNames=["Lat", "Lon", "Alt"]);
    altitude = sat.LLA.Alt/1e3;     % Altitude (km)

    % Check if altitude > 0km
    if(any(altitude<0))
        flag = 0;
        propData = [];
        f = errordlg("Orbit intersects the central body",'Orbit Propagation Error');
        return;
    end



    %% Eclipse check calculation
%     r_I = sat.r_I.Data;
    distance_Sat2MoonCenter = sqrt(r_I(:,1).^2+r_I(:,2).^2+r_I(:,3).^2)/1e3; % (km)
    thetaSat2MoonTangent = asind((moon.R_eq/1e3)./distance_Sat2MoonCenter);

    % Earth blocking (Eclipse)
    distance_Sat2EarthCenter = sqrt(earth_I(:,1).^2 + earth_I(:,2).^2 + earth_I(:,3).^2); % km
    thetaSat2EarthTangent = asind((earthRadius/1e3)./distance_Sat2EarthCenter);
    
    
    sunMagnitude = ones(length(time_Epoch),1);
    eclipseFlag = zeros(length(time_Epoch),1);
    for incl=1:length(time_Epoch)
        sat2sun_I = sun_I - r_I/1e3;
        sat2earth_I = earth_I - r_I/1e3;
        angle_Sat2Sun_Sat2Moon = calculateAngleBetweenVectors(sat2sun_I(incl,:),-r_I(incl,:));
        if angle_Sat2Sun_Sat2Moon > thetaSat2MoonTangent(incl)
            % Check if earth is blocking the Sun or not (Eclipse check)
            angle_Sat2Sun_Sat2Earth = calculateAngleBetweenVectors(sat2sun_I(incl,:),sat2earth_I(incl,:));
            if(angle_Sat2Sun_Sat2Earth > thetaSat2EarthTangent(incl))
                sunMagnitude(incl) = 1;
            else
                sunMagnitude(incl) = 0;    % Earth is blocking the Sun (Eclipse)
                eclipseFlag(incl) = 1;
            end
        else 
            sunMagnitude(incl) = 0;
        end
    end

    %% Albedo
    
    % Albedo Angle: angle(Moon->Sat, Moon->Sun)
    AlbedoAngle_deg_ = zeros(length(time_Epoch),1);
    satLat_deg = zeros(length(time_Epoch),1);
    for incl=1:length(time_Epoch)
        AlbedoAngle_deg_(incl) = calculateAngleBetweenVectors(r_I(incl,:), sun_I(incl,:));
        %Calculating satellite underneath latitude
        satLat_deg(incl) = 90-calculateAngleBetweenVectors(r_I(incl,:), [0,0,1]);
    end
    albedoAngle(:,1) = time_Epoch;
    albedoAngle(:,2) = AlbedoAngle_deg_;
    albedoAngle = array2table(albedoAngle);
    albedoAngle.Properties.VariableNames = {'time','AlbedoAngle_deg_'};
    albedoAngle.time = seconds(albedoAngle.time);
    albedoAngle = table2timetable(albedoAngle);

    %% Beta angle

    betaAngle = zeros(length(time_Epoch),1);
    for incl=1:length(time_Epoch)
        h_I = cross(r_I(incl,:)/norm(r_I(incl,:)), v_I(incl,:)/norm(v_I(incl,:))) ;% Angular momentum vector
        beta = calculateAngleBetweenVectors(h_I, sun_I(incl,:));
        betaAngle(incl) = 90-beta;
    end


    %% Return data
    flag = 1;
    % Position/Vel at Inertial Frame 
    propData.r_I = r_I;    % m
    propData.v_I = v_I;    % m/s
    % Position/Vel at Fixed Frame
    propData.r_ff = r_ff;   % m
    propData.v_ff = v_ff;     % m/s
    propData.sat = sat;
    propData.time_Epoch = time_Epoch;
    propData.time_UTC = time_UTC;
    propData.altitude = altitude;
    propData.Sun_I = sun_I;
    propData.Sun_ff = sun_ff; % unit vector
    propData.sunMagnitude = sunMagnitude;
    propData.albedoAngle = albedoAngle;
    propData.limbAngle = thetaSat2MoonTangent;
    propData.centralBodyRadius = moon.R_eq; 
    propData.betaAngle = betaAngle;
    propData.eclipseFlag = eclipseFlag;
    
catch ME
    warning(ME.message);
    if(ME.message)
        f = errordlg(ME.message,'Orbit Propagation Error');
    end
    flag = 0;
    propData = [];
end

end

function rotm = euler2rotmZXZ(phi, tet, psi)

    rotm_z1 = [ cos(phi)    sin(phi)    0;
                -sin(phi)   cos(phi)    0;
                0           0           1];
    
    rotm_x = [  1   0           0;
                0   cos(tet)    sin(tet);
                0   -sin(tet)   cos(tet)];
    
    rotm_z2 = [ cos(psi)    sin(psi)    0;
                -sin(psi)   cos(psi)    0;
                0           0           1];
    
    rotm = rotm_z1 * rotm_x * rotm_z2;
end

function rotm = euler2rotmZYX(phi, tet, psi)

    rotm_x = [  1   0           0;
                0   cos(psi)    sin(psi);
                0   -sin(psi)   cos(psi)];
    
    rotm_y = [  cos(tet)   0  -sin(tet);
                0          1          0;
                sin(tet)   0   cos(tet)];


    rotm_z = [ cos(phi)    sin(phi)    0;
                -sin(phi)   cos(phi)    0;
                0           0           1];
    
    
    rotm = rotm_x * rotm_y * rotm_z;
end

%% Process simulation Data
% Compute latitude, longitude, and altitude using lunar equatorial radius and flattening. Values are displayed in degrees and meters.
% sat.MEPos = [sat.TimetablePos.Data(:,1) ...
%     sat.TimetablePos.Data(:,2) sat.TimetablePos.Data(:,3)];
% lla = ecef2lla(sat.MEPos, moon.F, moon.R_eq);
% sat.LLA = timetable(sat.TimetablePos.Time, ...
%     lla(:,1), lla(:,2), lla(:,3), ...
%     VariableNames=["Lat", "Lon", "Alt"]);
% clear lla;
% disp(csm.LLA);

%% Results

% % Display Trajectory over the 3-D Moon
% figure; axis off; colormap gray; 
% axesm("globe","Geoid", moon.ReferenceEllipsoid);
% geoshow(moon.Data.moonalb20c, moon.Data.moonalb20cR, DisplayType="texturemap");
% plot3(sat.MEPos(:,1), sat.MEPos(:,2), sat.MEPos(:,3),"r");
% view(90,0);
% 
% % Display 2-D Projection
% figure; colormap gray;
% axesm(MapProjection="robinson");
% geoshow(moon.Data.moonalb20c, moon.Data.moonalb20cR, DisplayType="texturemap");
% plotm(sat.LLA.Lat, sat.LLA.Lon, Color="r");





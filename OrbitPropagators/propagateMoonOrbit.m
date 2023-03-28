function [propData, flag] = propagateMoonOrbit(orb)


%% Define Mission Parameters

mission.StartDate = datetime(year(orb.startDate),month(orb.startDate),day(orb.startDate),orb.startHour,orb.startMinute,0, "TimeZone","UTC");
mission.Duration  = hours(orb.durationHours);

sat.SemiMajorAxis  = orb.semiMajorAxis*1e3;     % [m]
sat.Eccentricity   = orb.eccentricity;
sat.Inclination    = orb.inclination;    % [deg]
sat.RAAN           = orb.RAAN;    % [deg]
sat.ArgOfPeriapsis = orb.RAAN;    % [deg]
sat.TrueAnomaly    = orb.trueAnomaly;      % [deg]

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

% Set CSM initial conditions. To assign the Keplerian orbital element set defined in the previous section, use setBlockParameter.
mission.sim.in = mission.sim.in.setBlockParameter(...
    sat.blk, "startDate", string(juliandate(mission.StartDate)),...
    sat.blk, "stateFormatNum", "Orbital elements",...
    sat.blk, "orbitType", "Keplerian",...
    sat.blk, "semiMajorAxis", string(sat.SemiMajorAxis),...
    sat.blk, "eccentricity", string(sat.Eccentricity),...
    sat.blk, "inclination", string(sat.Inclination),...
    sat.blk, "raan", string(sat.RAAN),...
    sat.blk, "argPeriapsis", string(sat.ArgOfPeriapsis),...
    sat.blk, "trueAnomaly", string(sat.TrueAnomaly));

% Set the position and velocity output ports of the block to use the Moon-fixed frame.  The fixed-frame for the Moon is the Mean Earth/Pole Axis (ME) reference system.
mission.sim.in = mission.sim.in.setBlockParameter(...
    sat.blk, "centralBody", "Moon",...
    sat.blk, "outportFrame", "ICRF");


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
%     assignin('base', 'moon', moon);
%     assignin('base', 'mission', mission);
    mission.sim.out = sim(mission.sim.in);
    
    % Extract the position and velocity data from the model output data structure.
    % Expected ref frame: Fixed-Frame always !!!! 
    % ### Fixed-frame ###
    sat.r_ff = mission.sim.out.yout{1}.Values;
    sat.r_ff.Data = sat.r_ff.Data./1e3; % m->km
    sat.v_ff = mission.sim.out.yout{2}.Values;


    % Time vectors (UTC, Epoch, Julian)
    time_Epoch =  seconds(sat.r_ff.Time);
    sat.r_ff.Properties.StartTime = mission.StartDate;
    sat.v_ff.Properties.StartTime = mission.StartDate;
    time_UTC = sat.r_ff.Time;
    time_Julian = juliandate(time_UTC);

    % Sun position x,y,z at Moon-Centered (Inertial) frame
    sun_I = planetEphemeris(time_Julian,'Moon','Sun'); % units: km


    % Convert to Moon Fixed Frame
    moonAngles = moonLibration(time_Julian, '421');


    for i=1:length(time_Julian)
        % Rotation from Fixed frame (Principal Axis) to Inertial
        rotm_PA2I = euler2rotmZXZ(-moonAngles(i,1), -moonAngles(i,2), -moonAngles(i,3));
        % Rotation from  Inertial to Fixed frame (Principal Axis)
        rotm_I2PA = rotm_PA2I';
        % Rotation from PA to mean-Earth/mean-rotation (MER) frame
        rotm_PA2ME = euler2rotmZYX(-deg2rad(dms2degrees([0 0 0.2785])), -deg2rad(dms2degrees([0 1 18.6944])), -deg2rad(dms2degrees([0 1 7.8526])));
        % Rotation from Inertial to MER
        rotm_I2ME = rotm_PA2ME * rotm_I2PA;
        rotm_ME2I = rotm_I2ME';
        % Convert position, velocity from inertial to fixed frame
        r_f = sat.r_ff.Data(i,:)'; v_f = sat.v_ff.Data(i,:)';
        r_i = rotm_ME2I * r_f;
        v_i = rotm_ME2I * v_f;
        r_I(i,:) = r_i';
        v_I(i,:) = v_i';
        % Sub-solar coordinates
        sun_vec = sun_I(i, :) / norm(sun_I(i, :));
        sun_ff(i,:) = rotm_I2ME * sun_vec';
    end


    % Compute latitude, longitude, and altitude using lunar equatorial radius and flattening. Values are displayed in degrees and meters.
    r_ff = sat.r_ff.Data;
    lla = ecef2lla(r_ff * 1e3, moon.F, moon.R_eq);
    sat.LLA = timetable(sat.r_ff.Time, ...
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
    distance_Sat2MoonCenter = sqrt(r_I(:,1).^2+r_I(:,2).^2+r_I(:,3).^2); % (km)
    thetaSat2MoonTangent = asind((moon.R_eq/1e3)./distance_Sat2MoonCenter);
    
    calculateAngleBetweenVectors = @(u,v) acosd(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1)); % function that return angle between two vectors in degrees
    
    sunMagnitude = ones(length(time_Epoch),1);
    for i=1:length(time_Epoch)
        sat2sun_I = sun_I - r_I;
        angle_Sat2Sun_Sat2Moon = calculateAngleBetweenVectors(sat2sun_I(i,:),-r_I(i,:));
        if angle_Sat2Sun_Sat2Moon > thetaSat2MoonTangent(i)
            sunMagnitude(i) = 1;
        else 
            sunMagnitude(i) = 0;
        end
    end

    %% Albedo
    
    % Albedo Angle: angle(Moon->Sat, Moon->Sun)
    AlbedoAngle_deg_ = zeros(length(time_Epoch),1);
    satLat_deg = zeros(length(time_Epoch),1);
    for i=1:length(time_Epoch)
        AlbedoAngle_deg_(i) = calculateAngleBetweenVectors(r_I(i,:), sun_I(i,:));
        %Calculating satellite underneath latitude
        satLat_deg(i) = 90-calculateAngleBetweenVectors(r_I(i,:), [0,0,1]);
    end
    albedoAngle(:,1) = time_Epoch;
    albedoAngle(:,2) = AlbedoAngle_deg_;
    albedoAngle = array2table(albedoAngle);
    albedoAngle.Properties.VariableNames = {'time','AlbedoAngle_deg_'};
    albedoAngle.time = seconds(albedoAngle.time);
    albedoAngle = table2timetable(albedoAngle);

    %% Beta angle

    betaAngle = zeros(length(time_Epoch),1);
    for i=1:length(time_Epoch)
        h_I = cross(r_I(i,:)/norm(r_I(i,:)), v_I(i,:)/norm(v_I(i,:))) ;% Angular momentum vector
        beta = calculateAngleBetweenVectors(h_I, sun_I(i,:));
        betaAngle(i) = 90-beta;
    end


    %% Return data
    flag = 1;
    % Position/Vel at Inertial Frame 
    propData.r_I = r_I*1e3;         % m
    propData.v_I = v_I;    % m/s
    % Position/Vel at Fixed Frame
    propData.r_ff = r_ff*1e3;   % m
    propData.v_ff = sat.v_ff.Data;     % m/s

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





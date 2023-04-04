function [propData, flag] = propagateEarthOrbit(orb)

    startTime = datetime(year(orb.startDate),month(orb.startDate),day(orb.startDate),orb.startHour,orb.startMinute,0, "TimeZone","UTC");
    stopTime = startTime + hours(orb.durationHours);
    
    sc = satelliteScenario(startTime,stopTime,orb.sampleTime);
    satName = "SAT";
    
    try
        % Prop satellite
        sat = satellite(sc, ...
            orb.semiMajorAxis*1e3,...
            orb.eccentricity,...
            orb.inclination, ...
            orb.RAAN,...
            orb.argumentOfPeriapsis,...
            orb.trueAnomaly,...
            'Name',satName,...
            "OrbitPropagator", orb.orbitPropagator);

        % Time vectors (UTC, Epoch, Julian)
        startTimeUTC = datetime(sc.StartTime,"TimeZone","UTC");
        stopTimeUTC = datetime(sc.StopTime,"TimeZone","UTC");
        time_Epoch = [0:sc.SampleTime:seconds(sc.StopTime - sc.StartTime)]';
        time_UTC = [startTimeUTC:seconds(sc.SampleTime):stopTimeUTC]';
        time_Julian = juliandate(time_UTC);
    
        % Satellite position x,y,z and velocity xdot, ydot, zdot at ECI (Inertial) frame
        [r_I,v_I] = states(sat,'CoordinateFrame', 'inertial'); %m
        r_I = (r_I')./1e3; % m -> km
        v_I = (v_I')./1e3; % m/s -> km/s
    
        % Sun position x,y,z at ECI (Inertial) frame
        [sun_I, ~] = planetEphemeris(time_Julian,'Earth','Sun'); % units: km

        % Convert sun vector from ECI to ECEF (Much faster than using the eci2ecef function)
        sun_ff = zeros(size(sun_I));
        A_ECI2ECEF = dcmeci2ecef('IAU-2000/2006',time_UTC); 
        for i=1:length(time_Julian)
            sun_f = A_ECI2ECEF(:,:,i) * sun_I(i,:)';
            sun_ff(i,:) = sun_f / norm(sun_f);
        end
        
        % Moon position x,y,z at ECI (Inertial) frame
        moon_I = planetEphemeris(time_Julian,'Earth','Moon'); % units: km
    
    
        %% Eclipse check calculation
        distance_Sat2EarthCenter = sqrt(r_I(:,1).^2+r_I(:,2).^2+r_I(:,3).^2); % (km)
        altitude = distance_Sat2EarthCenter - earthRadius/1e3;
        thetaSat2EarthTangent = asind((earthRadius/1e3)./distance_Sat2EarthCenter);
        
        calculateAngleBetweenVectors = @(u,v) acosd(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1)); % function that return angle between two vectors in degrees
        
        sunMagnitude = ones(length(time_Epoch),1);
        for i=1:length(time_Epoch)
            sat2sun_I = sun_I - r_I;
            angle_Sat2Sun_Sat2Earth = calculateAngleBetweenVectors(sat2sun_I(i,:),-r_I(i,:));
            if angle_Sat2Sun_Sat2Earth > thetaSat2EarthTangent(i)
                sunMagnitude(i) = 1;
            else 
                sunMagnitude(i) = 0;
            end
        end
    
        %% Albedo
    
        % Albedo Angle: angle(Earth->Sat, Earth->Sun)
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
        [r_I, v_I] = states(sat,'CoordinateFrame', 'inertial'); %m
        r_I = r_I'; v_I = v_I';
        [r_ecef, v_ecef] = states(sat,'CoordinateFrame', 'ecef'); %m

        betaAngle = zeros(length(time_Epoch),1);
        for i=1:length(time_Epoch)
            h_I = cross(r_I(i,:)/norm(r_I(i,:)), v_I(i,:)/norm(v_I(i,:))) ;% Angular momentum vector
            beta = calculateAngleBetweenVectors(h_I, sun_I(i,:));
            betaAngle(i) = 90-beta;
        end

        
        %% Return data


        flag = 1;
        % Position/Vel at Inertial Frame
        propData.r_I = r_I;     % m
        propData.v_I = v_I;     % m/s
        % Position/Vel at Fixed Frame
        propData.r_ff = r_ecef'; % m
        propData.v_ff = v_ecef'; % m/s
        %
        propData.sat = sat;
        propData.time_Epoch = time_Epoch;
        propData.time_UTC = time_UTC;
        propData.altitude = altitude;
        propData.Sun_I = sun_I;
        propData.Sun_ff = sun_ff; % unit vector
        propData.sunMagnitude = sunMagnitude;
        propData.Moon_I = moon_I;
        propData.albedoAngle = albedoAngle;
        propData.limbAngle = thetaSat2EarthTangent;
        propData.centralBodyRadius = earthRadius; 
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
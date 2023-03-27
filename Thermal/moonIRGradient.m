function faceIR = moonIRGradient(param)
% This function computs each face IR input integratrin the Moon's surface
% temperature distribution
% The moon's surface temperature distribution is computed by the equation 
% provided in the Human Landing System Lunar Thermal Analysis Guidebook, 
% volume HLS-UG-001. NASA, 2021. (Pag. 17)

    % Moon's surface emissivity
    moonEmissivity = 0.96;
    Rmoon = param.orb.prop.centralBodyRadius;
    
    % Satellite position, velocity at fixed frame
    r_ff = param.orb.prop.r_ff;
    v_ff = param.orb.prop.v_ff;
    
    % Sun position at fixed frame
    sun_ff = param.orb.prop.Sun_ff;
    
    % function that return angle between two vectors in degrees
    calculateAngleBetweenVectors = @(u,v) acosd(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1));
    
    % Moon discretization steps
    delta_lat = 1;                      % [degree]
    delta_long = 1;                     % [degree]
    LONG = -180:delta_long:180;
    LAT = -90:delta_lat:90;

    % Loop through time
    time = param.orb.prop.time_Epoch;
    IR = zeros(length(time), 6);
    for i=1:length(time)
         r_sat = r_ff(i,:)'; % sat position at time i
         v_sat = v_ff(i,:)'; % sat velocity at time i
         sun = sun_ff(i,:)'; % sun position at time i
         % Loop through longitude angle (360 degrees)
        for long=LONG
            % Loop through latitude angle (180 degrees)
            for lat = LAT  
                
                % Compute unit vector of the discrete area centered on the Moon
                i_patch = [cosd(lat)*cosd(long); cosd(lat)*sind(long); sind(lat)];

                % Surface patch area
                dA = (Rmoon)^2 * cosd(lat) * deg2rad(delta_lat) * deg2rad(delta_long);    % units: m^2
                
                % Check if satellite is visible to patch
                patch2sat = r_sat - i_patch * Rmoon;
                satDeclinationAngle = calculateAngleBetweenVectors(i_patch, patch2sat);
                if(satDeclinationAngle < 90 && dA~=0)
    
                    % Compute the subsolar latitude angle
                    Lat_ss = 90-calculateAngleBetweenVectors(sun, i_patch);
    
                    % Compute QIR from patch 
                    QIR = qIRsubsolar(Lat_ss, param.Albedo, param.orb.solarFlux, moonEmissivity, param.Tir_Night); % units: W/m2
    
                    % Compute view factors for each face
                    F = computeViewFactorMoonPatch2Face(r_sat, v_sat, i_patch, Rmoon, param);
    
                    % Compute accumulated IR flux at each face
                    for n=1:6
                        IR(i,n) = IR(i,n) + F(n) .* QIR * param.emissivity(n) * dA;
                    end
                end
            end
        end
    end

    % Transform it to timetable for each face
    time = param.orb.prop.time_Epoch;
    % X+
    faceIR.Xplus = array2table(IR(:,1));
    faceIR.Xplus.time = seconds(time); 
    faceIR.Xplus = table2timetable(faceIR.Xplus);
    % X-
    faceIR.Xminus = array2table(IR(:,2));
    faceIR.Xminus.time = seconds(time); 
    faceIR.Xminus = table2timetable(faceIR.Xminus);
    % Y+
    faceIR.Yplus = array2table(IR(:,3));
    faceIR.Yplus.time = seconds(time); 
    faceIR.Yplus = table2timetable(faceIR.Yplus);
    % Y-
    faceIR.Yminus = array2table(IR(:,4));
    faceIR.Yminus.time = seconds(time); 
    faceIR.Yminus = table2timetable(faceIR.Yminus);
    % Z+
    faceIR.Zplus = array2table(IR(:,5));
    faceIR.Zplus.time = seconds(time); 
    faceIR.Zplus = table2timetable(faceIR.Zplus);
    % Z-
    faceIR.Zminus = array2table(IR(:,6));
    faceIR.Zminus.time = seconds(time); 
    faceIR.Zminus = table2timetable(faceIR.Zminus);

end
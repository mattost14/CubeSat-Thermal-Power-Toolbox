function [thermal,flag] = computeThermal(param)
%COMPUTETHERMAL Summary of this function goes here
%   Detailed explanation goes here

    sigma = 5.67e-8; % Stefan-Boltzmann constant (W · m -2 · K -4 )
    param.krad = param.emissivity*sigma;

    % Mass paramters
    param.facesMass = param.facesMassDistribution * (1-param.internalNode2TotalMassRatio) * param.satelliteMass;
    param.internalNodeMass = param.internalNode2TotalMassRatio * param.satelliteMass;

    % ### Internal loads ###
    param.thermal.internalPwr = computeInternalPwr(param);

    % ### Compute external thermal loads ###
    % IR Radiation
    switch param.IRmodel
        case "Day/Night Temp"
            param.thermal.faceIR = planetIRsimple(param);
        case "Gradient"
            param.thermal.faceIR = moonIRGradient(param);
    end
    % Solar Radiation (Direct)
    param.thermal.faceSolar = computeSolarRadiation(param);
    % Solar Radiation (Albedo)
    param.thermal.faceAlbedo = computeAlbedoRadiation(param);

    % ######## RUN SIMSCAPE MODEL ########
    thermalModel = prepareSimpleThermalModelInputs(param);
    assignin('base', 'param', param)
    out = sim(thermalModel.sim.in);


    % Return results
    flag = 1;
    thermal.faceIR = param.thermal.faceIR;
    thermal.faceSolar = param.thermal.faceSolar;
    thermal.faceAlbedo = param.thermal.faceAlbedo;
    thermal.out = out;
end

function thermalModel = prepareSimpleThermalModelInputs(param)
    % ######## PREPARE SIMSCAPE simpleThermalModel MODEL ########
    thermalModel.mdl = "simpleThermalModel_v2";

    % Use a SimulationInput object to configure the model for our simulation.
    thermalModel.sim.in = Simulink.SimulationInput(thermalModel.mdl);

    % Solver setttings
    thermalModel.sim.in = thermalModel.sim.in.setModelParameter(...
        SolverType="Variable-step",...
        SolverName="VariableStepAuto",...
        RelTol="1e-6",...
        AbsTol="1e-7",...
        StopTime=string(param.orb.durationHours*3600));
    
    
    % Define thermal resistance inputs
    blk_RXplus = thermalModel.mdl + "/RX+";
    blk_RXminus = thermalModel.mdl + "/RX-";
    blk_RYplus = thermalModel.mdl + "/RY+";
    blk_RYminus = thermalModel.mdl + "/RY-";
    blk_RZplus = thermalModel.mdl + "/RZ+";
    blk_RZminus = thermalModel.mdl + "/RZ-";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_RXplus, "resistance", string(param.resistance(1)),...
    blk_RXminus, "resistance", string(param.resistance(2)),...
    blk_RYplus, "resistance", string(param.resistance(3)),...
    blk_RYminus, "resistance", string(param.resistance(4)),...
    blk_RZplus, "resistance", string(param.resistance(5)),...
    blk_RZminus, "resistance", string(param.resistance(6)));

    % Define inputs for the Internal Node 1
    blk_thermalMass = thermalModel.mdl + "/Internal Node 1";
    blk_internalNodePwr = thermalModel.mdl + "/Internal Node Pwr";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_thermalMass, "mass", string(param.internalNodeMass(1)),...
    blk_thermalMass, "sp_heat", string(param.internalNodeCp(1)),...
    blk_internalNodePwr, "VariableName", 'param.thermal.internalPwr');


    % Define inputs for Face X+ block
    faceId = 1;
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face X+/Radiative Heat Transfer";
    blk_thermalMass = thermalModel.mdl + "/Face X+/Thermal Mass";
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X+/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X+/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X+/Albedo Flux";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
    blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
    blk_thermalMass, "mass", string(param.facesMass(faceId)),...
    blk_thermalMass, "sp_heat", string(param.facesCp(faceId)),...
    blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Xplus',...
    blk_IRFlux, "VariableName", 'param.thermal.faceIR.Xplus',...
    blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Xplus');

    % Define inputs for Face X- block
    faceId = 2;
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face X-/Radiative Heat Transfer";
    blk_thermalMass = thermalModel.mdl + "/Face X-/Thermal Mass";
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X-/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X-/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment X-/Albedo Flux";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
    blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
    blk_thermalMass, "mass", string(param.facesMass(faceId)),...
    blk_thermalMass, "sp_heat", string(param.facesCp(faceId)),...
    blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Xminus',...
    blk_IRFlux, "VariableName", 'param.thermal.faceIR.Xminus',...
    blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Xminus');

    % Define inputs for Face Y+ block
    faceId = 3;
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face Y+/Radiative Heat Transfer";
    blk_thermalMass = thermalModel.mdl + "/Face Y+/Thermal Mass";
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y+/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y+/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y+/Albedo Flux";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
    blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
    blk_thermalMass, "mass", string(param.facesMass(faceId)),...
    blk_thermalMass, "sp_heat", string(param.facesCp(faceId)),...
    blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Yplus',...
    blk_IRFlux, "VariableName", 'param.thermal.faceIR.Yplus',...
    blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Yplus');

    % Define inputs for Face Y- block
    faceId = 4;
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face Y-/Radiative Heat Transfer";
    blk_thermalMass = thermalModel.mdl + "/Face Y-/Thermal Mass";
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y-/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y-/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Y-/Albedo Flux";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
    blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
    blk_thermalMass, "mass", string(param.facesMass(faceId)),...
    blk_thermalMass, "sp_heat", string(param.facesCp(faceId)),...
    blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Yminus',...
    blk_IRFlux, "VariableName", 'param.thermal.faceIR.Yminus',...
    blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Yminus');

    % Define inputs for Face Z+ block
    faceId = 5;
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face Z+/Radiative Heat Transfer";
    blk_thermalMass = thermalModel.mdl + "/Face Z+/Thermal Mass";
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z+/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z+/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z+/Albedo Flux";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
    blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
    blk_thermalMass, "mass", string(param.facesMass(faceId)),...
    blk_thermalMass, "sp_heat", string(param.facesCp(faceId)),...
    blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Zplus',...
    blk_IRFlux, "VariableName", 'param.thermal.faceIR.Zplus',...
    blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Zplus');

    % Define inputs for Face Z- block
    faceId = 6;
    blk_radiativeHeatTransfer = thermalModel.mdl + "/Face Z-/Radiative Heat Transfer";
    blk_thermalMass = thermalModel.mdl + "/Face Z-/Thermal Mass";
    blk_solarFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z-/Solar Flux";
    blk_IRFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z-/IR Flux";
    blk_AlbedoFlux =   thermalModel.mdl + "/Thermal Loads from Enviroment Z-/Albedo Flux";

    thermalModel.sim.in = thermalModel.sim.in.setBlockParameter(...
    blk_radiativeHeatTransfer, "area", string(param.facesArea(faceId)),...
    blk_radiativeHeatTransfer, "rad_tr_coeff", string(param.krad(faceId)),...
    blk_thermalMass, "mass", string(param.facesMass(faceId)),...
    blk_thermalMass, "sp_heat", string(param.facesCp(faceId)),...
    blk_solarFlux, "VariableName", 'param.thermal.faceSolar.Zminus',...
    blk_IRFlux, "VariableName", 'param.thermal.faceIR.Zminus',...
    blk_AlbedoFlux, "VariableName",'param.thermal.faceAlbedo.Zminus');
end

function internalPwr = computeInternalPwr(param)
    % This function computes the internal dissipation power profile based
    % on inputs in the power tab

    time = param.orb.prop.time_Epoch;
    pwrDissipationDay = param.pwrDissipationDaySimpleModel;
    pwrDissipationNight = param.pwrDissipationNightSimpleModel;
    sunMagnitude = param.orb.prop.sunMagnitude;
    switch param.pwrDissipationProfile
        case "Constant"
            internalPwr = pwrDissipationDay * ones(length(time),1);
        case "Day/Night"
            internalPwr = pwrDissipationDay .* sunMagnitude;
            internalPwr = internalPwr +  pwrDissipationNight .* (1-sunMagnitude);
    end

    % Internal power as timetable
    internalPwr = array2table(internalPwr);
    internalPwr.time = seconds(time); 
    internalPwr = table2timetable(internalPwr);
end


function faceIR = planetIRsimple(param)
    % Constants
    sigma = 5.67e-8; % Stefan-Boltzmann constant (W · m -2 · K -4 )

    %%% Planetary IR Power
    daytime = (param.orb.prop.albedoAngle.AlbedoAngle_deg_ - 90 < 0);
    Pir = daytime*sigma*(param.Tir_Day)^4 + (1-daytime)*sigma*(param.Tir_Night)^4; % Central body radiation

    % Compute IR flux at each face
    for n=1:length(param.facesMaterial)
        IR(:,n) = param.pwr.viewFactor(:,n) .* Pir * param.emissivity(n) * param.facesArea(n);
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
                    F = computeViewFactorPatch2Face(r_sat, v_sat, i_patch, Rmoon, param);
    
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

function F = computeViewFactorPatch2Face(r_sat, v_sat, i_patch, R, param)

    % function that return angle between two vectors in degrees
    calculateAngleBetweenVectors = @(u,v) acosd(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1));

    % Get face vector at LVLH frame
    XplusFaceVector_LVLH = param.pwr.XplusFaceVector_LVLH;
    YplusFaceVector_LVLH = param.pwr.YplusFaceVector_LVLH;
    ZplusFaceVector_LVLH = param.pwr.ZplusFaceVector_LVLH;

    % Get rotation matrix from LVLH to Fixed Frame
    z_L2F = -r_sat/norm(r_sat);   % Unit vector z
    h_F = cross(r_sat, v_sat);    % Angular momentum 
    y_L2F = -h_F/norm(h_F);             % Unit vector y
    x_L2F = cross(y_L2F, z_L2F);        % Unit vector x
    A_F2L = [x_L2F'; y_L2F'; z_L2F'];   % Rotation matrix from Fixed Frame to LVLH frame
    A_L2F = A_F2L';                     % Rotation matrix from LVLH frame to Fixed Frame
    
    % Get face vector at Fixed Frame
    XplusFaceVector_F = A_L2F * XplusFaceVector_LVLH;
    YplusFaceVector_F = A_L2F * YplusFaceVector_LVLH;
    ZplusFaceVector_F = A_L2F * ZplusFaceVector_LVLH;

    % View factor function (ref:
    % http://imartinez.etsiae.upm.es/~isidoro/tc3/Radiation%20View%20factors.pdf,
    % Eq. 1, Pag. 3)
    computeViewFactor = @(beta1, beta2, A2, r12) cosd(beta1)*cosd(beta2)*A2/(pi*r12^2);

    % ### Compute View Factor ###
    patch2sat = r_sat - i_patch * R;
    i_patch2sat = patch2sat/norm(patch2sat);
    beta1 = calculateAngleBetweenVectors(i_patch, patch2sat);

    beta2(1) =  calculateAngleBetweenVectors(XplusFaceVector_F, -i_patch2sat);
    beta2(2) =  calculateAngleBetweenVectors(-XplusFaceVector_F, -i_patch2sat);
    beta2(3) =  calculateAngleBetweenVectors(YplusFaceVector_F, -i_patch2sat);
    beta2(4) =  calculateAngleBetweenVectors(-YplusFaceVector_F, -i_patch2sat);
    beta2(5) =  calculateAngleBetweenVectors(ZplusFaceVector_F, -i_patch2sat);
    beta2(6) =  calculateAngleBetweenVectors(-ZplusFaceVector_F, -i_patch2sat);

    F = zeros(6,1);
    for n=1:6
        if(beta2(n)<90)
            F(n) = computeViewFactor(beta1, beta2(n), param.facesArea(n), norm(patch2sat));
        end
    end
end

function qIR = qIRsubsolar(Latss, albedo, solarFlux, moonEmissivity, Tdark)
% The moon's surface temperature distribution is computed by the equation 
% provided in the Human Landing System Lunar Thermal Analysis Guidebook, 
% volume HLS-UG-001. NASA, 2021. (Pag. 17)
% Latss : degrees latitude from subsolar point (Subsolar Coordinate System)
% albedo: Lunar surface albedo
% solarFlux: solar flux (W/m2)
% moonEmissivity: Lunar Surface Emissivity
% Tdark: unilluminated lunar surface temperature

    sigma = 5.670374419e-8; % W*m2*K^-4
    if(Latss >= 0)
        qIR = sind(Latss)*((1-albedo)*solarFlux - sigma * moonEmissivity * Tdark^4)+  sigma * moonEmissivity * Tdark^4;
    else
        qIR = sigma * moonEmissivity * Tdark^4;
    end
end


function faceSolar = computeSolarRadiation(param)
    solarFlux = param.orb.solarFlux;

    % Compute solar flux at each face
    for n=1:length(param.facesMaterial)
        if(param.facesMaterial(n)=="Solar Panel") % If face is solar cell, check whether it is generating electrical energy or it is dissipating outside 
            effectiveAbsoptivity =  param.absorptivity(n) - param.pwr.electricalEfficiency * param.solarCellEff(n);
            solar(:,n) = param.pwr.solarLight.(n) .* solarFlux .* effectiveAbsoptivity * param.facesArea(n);
        else
            solar(:,n) = param.pwr.solarLight.(n) .* solarFlux * param.absorptivity(n) * param.facesArea(n);
        end
    end

    % Transform it to timetable for each face
    time = param.orb.prop.time_Epoch;
    % X+
    faceSolar.Xplus = array2table(solar(:,1));
    faceSolar.Xplus.time = seconds(time); 
    faceSolar.Xplus = table2timetable(faceSolar.Xplus);
    % X-
    faceSolar.Xminus = array2table(solar(:,2));
    faceSolar.Xminus.time = seconds(time); 
    faceSolar.Xminus = table2timetable(faceSolar.Xminus);
    % Y+
    faceSolar.Yplus = array2table(solar(:,3));
    faceSolar.Yplus.time = seconds(time); 
    faceSolar.Yplus = table2timetable(faceSolar.Yplus);
    % Y-
    faceSolar.Yminus = array2table(solar(:,4));
    faceSolar.Yminus.time = seconds(time); 
    faceSolar.Yminus = table2timetable(faceSolar.Yminus);
    % Z+
    faceSolar.Zplus = array2table(solar(:,5));
    faceSolar.Zplus.time = seconds(time); 
    faceSolar.Zplus = table2timetable(faceSolar.Zplus);
    % Z-
    faceSolar.Zminus = array2table(solar(:,6));
    faceSolar.Zminus.time = seconds(time); 
    faceSolar.Zminus = table2timetable(faceSolar.Zminus);
end

function faceAlbedo = computeAlbedoRadiation(param)
    solarFlux = param.orb.solarFlux;
    albedoFactor = param.Albedo;
    
    % Compute solar flux at each face
    for n=1:length(param.facesMaterial)
        % Albedo view factor (Falb) = Fir * cos(SZA)
        Falb = param.pwr.viewFactor(:,n) .* max(0, cosd(param.orb.prop.albedoAngle.(1)));

        if(param.facesMaterial(n)=="Solar Panel") % If face is solar cell, check whether it is generating electrical energy or it is dissipating outside 
            effectiveAbsoptivity =  param.absorptivity(n) - param.pwr.electricalEfficiency * param.solarCellEff(n);
            albedoRadiation(:,n) = albedoFactor * Falb .* solarFlux .* effectiveAbsoptivity * param.facesArea(n);
        else
            albedoRadiation(:,n) = albedoFactor * Falb .* solarFlux * param.absorptivity(n) * param.facesArea(n);
        end
    end

    % Transform it to timetable for each face
    time = param.orb.prop.time_Epoch;
    % X+
    faceAlbedo.Xplus = array2table(albedoRadiation(:,1));
    faceAlbedo.Xplus.time = seconds(time); 
    faceAlbedo.Xplus = table2timetable(faceAlbedo.Xplus);
    % X-
    faceAlbedo.Xminus = array2table(albedoRadiation(:,2));
    faceAlbedo.Xminus.time = seconds(time); 
    faceAlbedo.Xminus = table2timetable(faceAlbedo.Xminus);
    % Y+
    faceAlbedo.Yplus = array2table(albedoRadiation(:,3));
    faceAlbedo.Yplus.time = seconds(time); 
    faceAlbedo.Yplus = table2timetable(faceAlbedo.Yplus);
    % Y-
    faceAlbedo.Yminus = array2table(albedoRadiation(:,4));
    faceAlbedo.Yminus.time = seconds(time); 
    faceAlbedo.Yminus = table2timetable(faceAlbedo.Yminus);
    % Z+
    faceAlbedo.Zplus = array2table(albedoRadiation(:,5));
    faceAlbedo.Zplus.time = seconds(time); 
    faceAlbedo.Zplus = table2timetable(faceAlbedo.Zplus);
    % Z-
    faceAlbedo.Zminus = array2table(albedoRadiation(:,6));
    faceAlbedo.Zminus.time = seconds(time); 
    faceAlbedo.Zminus = table2timetable(faceAlbedo.Zminus);
end
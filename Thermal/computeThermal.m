function [thermal,flag] = computeThermal(param)
%COMPUTETHERMAL Summary of this function goes here
%   Detailed explanation goes here

%     sigma = 5.67e-8; % Stefan-Boltzmann constant (W · m -2 · K -4 )
%     param.krad = param.emissivity*sigma;
% 
%     % Mass paramters
%     param.facesMass = param.facesMassDistribution * (1-param.internalNode2TotalMassRatio) * param.satelliteMass;
%     param.internalNodeMass = param.internalNode2TotalMassRatio * param.satelliteMass;

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
%     thermalModel = prepareSimscapeModelInputs(param);
    thermalModel = prepareSimscapeModelInputs_v2(param);
    assignin('base', 'param', param)
    out = sim(thermalModel.sim.in);


    % Return results
    flag = 1;
    thermal.faceIR = param.thermal.faceIR;
    thermal.faceSolar = param.thermal.faceSolar;
    thermal.faceAlbedo = param.thermal.faceAlbedo;
    thermal.out = out;
end













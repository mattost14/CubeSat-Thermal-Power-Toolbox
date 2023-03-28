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

    % ### Blanket Setup ###
    % If neighboor faces has the blanket material, then set a low thermal
    % resistance between faces
    blanketMaterialName = "Gold coating";
    faceText = ["X+","X-","Y+","Y-","Z+","Z-"];
    for n=1:6
        if(param.facesMaterial(n) == blanketMaterialName) % if the face material is "MLI - Blanket"
            face = char(faceText(n));
            param.resistanceBlanket(n) = 5;    % K/W (High thermal resistance crossing the blanket)
            for otherFacesID=1:6
                otherFace = char(faceText(otherFacesID));
                % If the neighboor face is also covered with MLI - Blanket
                if(otherFace(1)~=face(1) && param.facesMaterial(otherFacesID) == blanketMaterialName)
                    param.resistanceInterFacesBlanket(n,otherFacesID) = .1;        % K/W (Low thermal resistance along the blanket surface)
                elseif(otherFace(1)==face(1))
                    param.resistanceInterFacesBlanket(n,otherFacesID) = nan;        % Same or oposite face (Not applicable)
                else
                    param.resistanceInterFacesBlanket(n,otherFacesID) = 1e3;       % K/W (High thermal resistance)
                end
            end
        else
            param.resistanceBlanket(n) = 0;    % K/W (No blanket -> no resistance)
        end
    end

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













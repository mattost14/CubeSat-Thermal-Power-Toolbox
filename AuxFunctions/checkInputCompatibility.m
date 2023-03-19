function inputCompatible = checkInputCompatibility(sim)
%CHECKINPUTCOMPATABILITY Summary of this function goes here
%   Detailed explanation goes here
    requiredFields = [
        "orb"
        "Tir_Day"
        "Tir_Night"
        "Albedo"
        "earthIRmodelList"
        "moonIRmodelList"
        "IRmodel"
        "moonOrbPropagatorList"
        "earthOrbPropagatorList"
        "plot"
        "numberOfUnits"
        "facesMaterial"
        "satelliteMass"
        "solarCellEff"
        "effectiveAreaRatio"
        "resistance"
        "internalNodeCp"
        "internalNode2TotalMassRatio"
        "facesMassDistribution"
        "facesArea"
        "attitude"
        "deployable"
        "absorptivity"
        "emissivity"
        "facesCp"
        "pwrDissipationProfileList"
        "pwrDissipationProfile"
        "pwrDissipationDaySimpleModel"
        "pwrDissipationNightSimpleModel"
        "pwr"
        "thermal"];

        inputCompatible = 1; 
        for i=1:length(requiredFields)
            if(~isfield(sim,requiredFields(i)))
                inputCompatible = 0; 
                errordlg('Invalid input file');
                break
            end
        end
end


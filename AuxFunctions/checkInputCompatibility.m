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
    "IRmodelItemList"
    "IRmodel"
    "moonOrbPropagatorList"
    "earthOrbPropagatorList"
    "propagatorItemList"
    "plot"
    "numberOfUnits"
    "facesMaterial"
    "satelliteMass"
    "solarCellEff"
    "effectiveAreaRatio"
    "conductionResistance1U"
    "resistance"
    "internalNodeCp"
    "internalNode2TotalMassRatio"
    "facesMassDistribution"
    "facesArea"
    "attitude"
    "deployable"
    "resistanceBlanket"
    "facesBlanketMass"
    "facesBlanketCp"
    "kradFace2Blanket"
    "resistanceInterFacesBlanket"
    "absorptivity"
    "emissivity"
    "facesCp"
    "resistanceInterFaces"
    "pwrDissipationProfileList"
    "pwrDissipationProfile"
    "pwrDissipationDaySimpleModel"
    "pwrDissipationNightSimpleModel"
    "pwr"
    "computePwrWithShadowFlag"
    "thermal"
    "appVersion"];

        inputCompatible = 1; 
        for i=1:length(requiredFields)
            if(~isfield(sim,requiredFields(i)))
                inputCompatible = 0; 
                errordlg('Invalid input file');
                break
            end
        end
end


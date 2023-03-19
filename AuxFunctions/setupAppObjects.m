function setupAppObjects(app)
%SETUPAPPOBJECTS Summary of this function goes here
%   Detailed explanation goes here
    
    %% ### Orb Tab ###
    app.CentralBodyDropDown.Value = app.param.orb.centralBody;                  % Central Body
    app.SolarfluxWm2EditField.Value = app.param.orb.solarFlux;                  % Solar Flux
    app.IRmodelDropDown.Items = app.param.IRmodelItemList;                      % IR model items
    app.IRmodelDropDown.Value = app.param.IRmodel;                              % IR model
    app.DayTempKEditField.Value = app.param.Tir_Day;                            % T Day (K)
    app.NightTempKEditField.Value = app.param.Tir_Night;                        % T Day (K)
    app.EpochDatePicker.Value = app.param.orb.startDate;                        % Epoch day
    app.HSpinner.Value = app.param.orb.startHour;                               % Epoch hour
    app.MSpinner.Value = app.param.orb.startMinute;                             % Epoch minute
    app.DurationhoursSpinner.Value = app.param.orb.durationHours;               % Duration (hours)
    app.SemimajoraxiskmEditField.Value = app.param.orb.semiMajorAxis;           % Semi-major axis (km)
    app.EccentricityEditField.Value = app.param.orb.eccentricity;               % Eccentricity
    app.InclinationEditField.Value = app.param.orb.inclination;                 % Inclination
    app.RAANEditField.Value = app.param.orb.RAAN;                               % RAAN
    app.ArgumentofPeriapsisEditField.Value = app.param.orb.argumentOfPeriapsis; % Argument of Periapsis
    app.TrueAnomalyEditField.Value = app.param.orb.trueAnomaly;                 % True anomaly
    app.CheckBoxShowSunVector.Value = app.param.plot.orbitPlot.showSunVectorFlag;   % Sun vector flag
    app.PropagatorDropDown.Items = app.param.propagatorItemList;                % Propagator Items
    app.PropagatorDropDown.Value = app.param.orb.orbitPropagator;               % Propagator
    switch app.IRmodelDropDown.Value
        case "Day/Night Temp"
            app.DayTempKEditField.Visible = "on";
            app.DayTempKEditFieldLabel.Visible = "on";
        case "Gradient"
            app.DayTempKEditField.Visible = "off";
            app.DayTempKEditFieldLabel.Visible = "off";
    end

    % Flags goes to zero (require the user to re-run orbit and power 
    app.param.orb.prop.flag = 0;
    app.param.pwr.flag = 0; 
    app.UpdateOrbitStatusLabel.Text = '';
    %% ### Model Tab ###
    app.UnitsDropDown.Value = strcat(num2str(app.param.numberOfUnits),'U');     % #Units
    app.MasskgSpinner.Value = app.param.satelliteMass;                          % Mass (kg)
    % Internal node
    app.IntermalnodetototalmassratioSpinner.Value = app.param.internalNode2TotalMassRatio*100;  % Internal Node Mass Ratio
    app.IntermalnodeSpecificHeatJkgKSpinner.Value = app.param.internalNodeCp;                   % Internal Node Cp
    % Attitude
    app.NadirFaceDropDown.Value = app.param.attitude.nadirFace;                 % Nadir face
    app.RamFaceDropDown.Items = app.param.attitude.ramFaceItems;                % Ram face items
    app.RamFaceDropDown.Value = app.param.attitude.ramFace;                     % Ram face 
    % Resistance
    app.SpinnerRXplus.Value =   app.param.resistance(1); 
    app.SpinnerRXminus.Value =  app.param.resistance(2); 
    app.SpinnerRYplus.Value =   app.param.resistance(3); 
    app.SpinnerRYminus.Value =  app.param.resistance(4); 
    app.SpinnerRZplus.Value =   app.param.resistance(5); 
    app.SpinnerRZminus.Value =  app.param.resistance(6); 
    %% ### X+ Face
    faceID = 1;
    app.SurfaceTypeDropDownXplus.Value = app.param.facesMaterial(faceID);                           % Surface type
    app.CellEfficiencySpinnerXplus.Value = app.param.solarCellEff(faceID)*100;                      % Cell efficiency
    app.EffecAreaRatioSpinnerXplus.Value = app.param.effectiveAreaRatio(faceID)*100;                % Packing factor
    app.AbsorptivitySpinnerXplus.Value = app.param.absorptivity(faceID);                            % Absorptivity
    app.EmissivitySpinnerXplus.Value = app.param.emissivity(faceID);                                % Emissivity
    app.RatioEditFieldXplus.Value = app.param.absorptivity(faceID)/app.param.emissivity(faceID);    % Ratio abs/emiss
    app.CpSpinnerXplus.Value = app.param.facesCp(faceID);                                           % Face Cp  
    app.HingeDropDownXplus.Value = app.param.deployable.hinge(faceID);                              % Hinge 
    app.SizeDropDownXplus.Value = app.param.deployable.size(faceID);                                % Size
    app.FlipPanelCheckBoxXplus.Value = app.param.deployable.flipFixedPanel(faceID);                 % Flip panel
    app.CellEfficiencySpinnerXplusDeployable.Value = app.param.deployable.solarCellEff(faceID)*100; % Cell Efficiency (Deployable)
    app.EffecAreaRatioSpinnerXplusDeployable.Value = app.param.deployable.effectiveAreaRatio(faceID)*100;% Packing factor (Deployable)
    switch app.param.deployable.type(faceID)                                                        % Deployable type
        case "None"
            app.PanelTypeButtonGroupXplus.SelectedObject = app.NoneButtonXplus;
            app.NoneButtonXplus.Value = 1;
        case "Fixed"
            app.PanelTypeButtonGroupXplus.SelectedObject = app.FixedButtonXplus;
            app.FixedButtonXplus.Value = 1;
        case "Tracking"
            app.PanelTypeButtonGroupXplus.SelectedObject = app.TrackingButtonXplus;
            app.TrackingButtonXplus.Value = 1;
    end
    if(app.SurfaceTypeDropDownXplus.Value == "Solar Panel")
        app.CellEfficiencySpinnerXplus.Visible = "on";
        app.CellEfficiencySpinnerLabelXplus.Visible = "on";
        app.EffecAreaRatioSpinnerXplus.Visible = "on";
        app.EffecAreaRatioSpinnerLabelXplus.Visible = "on";
    else
        app.CellEfficiencySpinnerXplus.Visible = "off";
        app.CellEfficiencySpinnerLabelXplus.Visible = "off";
        app.EffecAreaRatioSpinnerXplus.Visible = "off";
        app.EffecAreaRatioSpinnerLabelXplus.Visible = "off";
    end
    if(app.SurfaceTypeDropDownXplus.Value == "Custom")
        app.EmissivitySpinnerXplus.Enable = "on";
        app.AbsorptivitySpinnerXplus.Enable = "on";
        app.CpSpinnerXplus.Enable = "on";
    else
        app.EmissivitySpinnerXplus.Enable = "off";
        app.AbsorptivitySpinnerXplus.Enable = "off";
        app.CpSpinnerXplus.Enable = "off";
    end
    if(app.PanelTypeButtonGroupXplus.SelectedObject == "Fixed")
        app.FlipPanelCheckBoxXplus.Enable = "on";
    else
        app.FlipPanelCheckBoxXplus.Enable = "off";
    end
    if((app.PanelTypeButtonGroupXplus.SelectedObject == "Fixed") || (app.PanelTypeButtonGroupXplus.SelectedObject == "Tracking"))
        app.HingeDropDownXplus.Enable = "on";
        app.SizeDropDownXplus.Enable = "on";
        app.CellEfficiencySpinnerXplusDeployable.Enable = "on";
        app.EffecAreaRatioSpinnerXplusDeployable.Enable = "on";
    else
        app.HingeDropDownXplus.Enable = "off";
        app.SizeDropDownXplus.Enable = "off";
        app.CellEfficiencySpinnerXplusDeployable.Enable = "off";
        app.EffecAreaRatioSpinnerXplusDeployable.Enable = "off";
        app.FlipPanelCheckBoxXplus.Enable = "off";
    end

    %% ### X- Face
    faceID = 2;
    app.SurfaceTypeDropDownXminus.Value = app.param.facesMaterial(faceID);                           % Surface type
    app.CellEfficiencySpinnerXminus.Value = app.param.solarCellEff(faceID)*100;                      % Cell efficiency
    app.EffecAreaRatioSpinnerXminus.Value = app.param.effectiveAreaRatio(faceID)*100;                % Packing factor
    app.AbsorptivitySpinnerXminus.Value = app.param.absorptivity(faceID);                            % Absorptivity
    app.EmissivitySpinnerXminus.Value = app.param.emissivity(faceID);                                % Emissivity
    app.RatioEditFieldXminus.Value = app.param.absorptivity(faceID)/app.param.emissivity(faceID);    % Ratio abs/emiss
    app.CpSpinnerXminus.Value = app.param.facesCp(faceID);                                           % Face Cp  
    app.HingeDropDownXminus.Value = app.param.deployable.hinge(faceID);                              % Hinge 
    app.SizeDropDownXminus.Value = app.param.deployable.size(faceID);                                % Size
    app.FlipPanelCheckBoxXminus.Value = app.param.deployable.flipFixedPanel(faceID);                 % Flip panel
    app.CellEfficiencySpinnerXminusDeployable.Value = app.param.deployable.solarCellEff(faceID)*100; % Cell Efficiency (Deployable)
    app.EffecAreaRatioSpinnerXminusDeployable.Value = app.param.deployable.effectiveAreaRatio(faceID)*100;% Packing factor (Deployable)
    switch app.param.deployable.type(faceID)                                                        % Deployable type
        case "None"
            app.PanelTypeButtonGroupXminus.SelectedObject = app.NoneButtonXminus;
            app.NoneButtonXminus.Value = 1;
        case "Fixed"
            app.PanelTypeButtonGroupXminus.SelectedObject = app.FixedButtonXminus;
            app.FixedButtonXminus.Value = 1;
        case "Tracking"
            app.PanelTypeButtonGroupXminus.SelectedObject = app.TrackingButtonXminus;
            app.TrackingButtonXminus.Value = 1;
    end
    if(app.SurfaceTypeDropDownXminus.Value == "Solar Panel")
        app.CellEfficiencySpinnerXminus.Visible = "on";
        app.CellEfficiencySpinnerLabelXminus.Visible = "on";
        app.EffecAreaRatioSpinnerXminus.Visible = "on";
        app.EffecAreaRatioSpinnerLabelXminus.Visible = "on";
    else
        app.CellEfficiencySpinnerXminus.Visible = "off";
        app.CellEfficiencySpinnerLabelXminus.Visible = "off";
        app.EffecAreaRatioSpinnerXminus.Visible = "off";
        app.EffecAreaRatioSpinnerLabelXminus.Visible = "off";
    end
    if(app.SurfaceTypeDropDownXminus.Value == "Custom")
        app.EmissivitySpinnerXminus.Enable = "on";
        app.AbsorptivitySpinnerXminus.Enable = "on";
        app.CpSpinnerXminus.Enable = "on";
    else
        app.EmissivitySpinnerXminus.Enable = "off";
        app.AbsorptivitySpinnerXminus.Enable = "off";
        app.CpSpinnerXminus.Enable = "off";
    end
    if(app.PanelTypeButtonGroupXminus.SelectedObject == "Fixed")
        app.FlipPanelCheckBoxXminus.Enable = "on";
    else
        app.FlipPanelCheckBoxXminus.Enable = "off";
    end
    if((app.PanelTypeButtonGroupXminus.SelectedObject == "Fixed") || (app.PanelTypeButtonGroupXminus.SelectedObject == "Tracking"))
        app.HingeDropDownXminus.Enable = "on";
        app.SizeDropDownXminus.Enable = "on";
        app.CellEfficiencySpinnerXminusDeployable.Enable = "on";
        app.EffecAreaRatioSpinnerXminusDeployable.Enable = "on";
    else
        app.HingeDropDownXminus.Enable = "off";
        app.SizeDropDownXminus.Enable = "off";
        app.CellEfficiencySpinnerXminusDeployable.Enable = "off";
        app.EffecAreaRatioSpinnerXminusDeployable.Enable = "off";
        app.FlipPanelCheckBoxXminus.Enable = "off";
    end
    %% ### Y+ Face
    faceID = 3;
    app.SurfaceTypeDropDownYplus.Value = app.param.facesMaterial(faceID);                           % Surface type
    app.CellEfficiencySpinnerYplus.Value = app.param.solarCellEff(faceID)*100;                      % Cell efficiency
    app.EffecAreaRatioSpinnerYplus.Value = app.param.effectiveAreaRatio(faceID)*100;                % Packing factor
    app.AbsorptivitySpinnerYplus.Value = app.param.absorptivity(faceID);                            % Absorptivity
    app.EmissivitySpinnerYplus.Value = app.param.emissivity(faceID);                                % Emissivity
    app.RatioEditFieldYplus.Value = app.param.absorptivity(faceID)/app.param.emissivity(faceID);    % Ratio abs/emiss
    app.CpSpinnerYplus.Value = app.param.facesCp(faceID);                                           % Face Cp  
    app.HingeDropDownYplus.Value = app.param.deployable.hinge(faceID);                              % Hinge 
    app.SizeDropDownYplus.Value = app.param.deployable.size(faceID);                                % Size
    app.FlipPanelCheckBoxYplus.Value = app.param.deployable.flipFixedPanel(faceID);                 % Flip panel
    app.CellEfficiencySpinnerYplusDeployable.Value = app.param.deployable.solarCellEff(faceID)*100; % Cell Efficiency (Deployable)
    app.EffecAreaRatioSpinnerYplusDeployable.Value = app.param.deployable.effectiveAreaRatio(faceID)*100;% Packing factor (Deployable)
    switch app.param.deployable.type(faceID)                                                        % Deployable type
        case "None"
            app.PanelTypeButtonGroupYplus.SelectedObject = app.NoneButtonYplus;
            app.NoneButtonYplus.Value = 1;
        case "FiYed"
            app.PanelTypeButtonGroupYplus.SelectedObject = app.FiYedButtonYplus;
            app.FiYedButtonYplus.Value = 1;
        case "Tracking"
            app.PanelTypeButtonGroupYplus.SelectedObject = app.TrackingButtonYplus;
            app.TrackingButtonYplus.Value = 1;
    end
    if(app.SurfaceTypeDropDownYplus.Value == "Solar Panel")
        app.CellEfficiencySpinnerYplus.Visible = "on";
        app.CellEfficiencySpinnerLabelYplus.Visible = "on";
        app.EffecAreaRatioSpinnerYplus.Visible = "on";
        app.EffecAreaRatioSpinnerLabelYplus.Visible = "on";
    else
        app.CellEfficiencySpinnerYplus.Visible = "off";
        app.CellEfficiencySpinnerLabelYplus.Visible = "off";
        app.EffecAreaRatioSpinnerYplus.Visible = "off";
        app.EffecAreaRatioSpinnerLabelYplus.Visible = "off";
    end
    if(app.SurfaceTypeDropDownYplus.Value == "Custom")
        app.EmissivitySpinnerYplus.Enable = "on";
        app.AbsorptivitySpinnerYplus.Enable = "on";
        app.CpSpinnerYplus.Enable = "on";
    else
        app.EmissivitySpinnerYplus.Enable = "off";
        app.AbsorptivitySpinnerYplus.Enable = "off";
        app.CpSpinnerYplus.Enable = "off";
    end
    if(app.PanelTypeButtonGroupYplus.SelectedObject == "FiYed")
        app.FlipPanelCheckBoxYplus.Enable = "on";
    else
        app.FlipPanelCheckBoxYplus.Enable = "off";
    end
    if((app.PanelTypeButtonGroupYplus.SelectedObject == "FiYed") || (app.PanelTypeButtonGroupYplus.SelectedObject == "Tracking"))
        app.HingeDropDownYplus.Enable = "on";
        app.SizeDropDownYplus.Enable = "on";
        app.CellEfficiencySpinnerYplusDeployable.Enable = "on";
        app.EffecAreaRatioSpinnerYplusDeployable.Enable = "on";
    else
        app.HingeDropDownYplus.Enable = "off";
        app.SizeDropDownYplus.Enable = "off";
        app.CellEfficiencySpinnerYplusDeployable.Enable = "off";
        app.EffecAreaRatioSpinnerYplusDeployable.Enable = "off";
        app.FlipPanelCheckBoxYplus.Enable = "off";
    end

    %% ### Y- Face
    faceID = 4;
    app.SurfaceTypeDropDownYminus.Value = app.param.facesMaterial(faceID);                           % Surface type
    app.CellEfficiencySpinnerYminus.Value = app.param.solarCellEff(faceID)*100;                      % Cell efficiency
    app.EffecAreaRatioSpinnerYminus.Value = app.param.effectiveAreaRatio(faceID)*100;                % Packing factor
    app.AbsorptivitySpinnerYminus.Value = app.param.absorptivity(faceID);                            % Absorptivity
    app.EmissivitySpinnerYminus.Value = app.param.emissivity(faceID);                                % Emissivity
    app.RatioEditFieldYminus.Value = app.param.absorptivity(faceID)/app.param.emissivity(faceID);    % Ratio abs/emiss
    app.CpSpinnerYminus.Value = app.param.facesCp(faceID);                                           % Face Cp  
    app.HingeDropDownYminus.Value = app.param.deployable.hinge(faceID);                              % Hinge 
    app.SizeDropDownYminus.Value = app.param.deployable.size(faceID);                                % Size
    app.FlipPanelCheckBoxYminus.Value = app.param.deployable.flipFixedPanel(faceID);                 % Flip panel
    app.CellEfficiencySpinnerYminusDeployable.Value = app.param.deployable.solarCellEff(faceID)*100; % Cell Efficiency (Deployable)
    app.EffecAreaRatioSpinnerYminusDeployable.Value = app.param.deployable.effectiveAreaRatio(faceID)*100;% Packing factor (Deployable)
    switch app.param.deployable.type(faceID)                                                        % Deployable type
        case "None"
            app.PanelTypeButtonGroupYminus.SelectedObject = app.NoneButtonYminus;
            app.NoneButtonYminus.Value = 1;
        case "FiYed"
            app.PanelTypeButtonGroupYminus.SelectedObject = app.FiYedButtonYminus;
            app.FiYedButtonYminus.Value = 1;
        case "Tracking"
            app.PanelTypeButtonGroupYminus.SelectedObject = app.TrackingButtonYminus;
            app.TrackingButtonYminus.Value = 1;
    end
    if(app.SurfaceTypeDropDownYminus.Value == "Solar Panel")
        app.CellEfficiencySpinnerYminus.Visible = "on";
        app.CellEfficiencySpinnerLabelYminus.Visible = "on";
        app.EffecAreaRatioSpinnerYminus.Visible = "on";
        app.EffecAreaRatioSpinnerLabelYminus.Visible = "on";
    else
        app.CellEfficiencySpinnerYminus.Visible = "off";
        app.CellEfficiencySpinnerLabelYminus.Visible = "off";
        app.EffecAreaRatioSpinnerYminus.Visible = "off";
        app.EffecAreaRatioSpinnerLabelYminus.Visible = "off";
    end
    if(app.SurfaceTypeDropDownYminus.Value == "Custom")
        app.EmissivitySpinnerYminus.Enable = "on";
        app.AbsorptivitySpinnerYminus.Enable = "on";
        app.CpSpinnerYminus.Enable = "on";
    else
        app.EmissivitySpinnerYminus.Enable = "off";
        app.AbsorptivitySpinnerYminus.Enable = "off";
        app.CpSpinnerYminus.Enable = "off";
    end
    if(app.PanelTypeButtonGroupYminus.SelectedObject == "FiYed")
        app.FlipPanelCheckBoxYminus.Enable = "on";
    else
        app.FlipPanelCheckBoxYminus.Enable = "off";
    end
    if((app.PanelTypeButtonGroupYminus.SelectedObject == "FiYed") || (app.PanelTypeButtonGroupYminus.SelectedObject == "Tracking"))
        app.HingeDropDownYminus.Enable = "on";
        app.SizeDropDownYminus.Enable = "on";
        app.CellEfficiencySpinnerYminusDeployable.Enable = "on";
        app.EffecAreaRatioSpinnerYminusDeployable.Enable = "on";
    else
        app.HingeDropDownYminus.Enable = "off";
        app.SizeDropDownYminus.Enable = "off";
        app.CellEfficiencySpinnerYminusDeployable.Enable = "off";
        app.EffecAreaRatioSpinnerYminusDeployable.Enable = "off";
        app.FlipPanelCheckBoxYminus.Enable = "off";
    end

    %% ### Z+ Face
    faceID = 5;
    app.SurfaceTypeDropDownZplus.Value = app.param.facesMaterial(faceID);                           % Surface type
    app.CellEfficiencySpinnerZplus.Value = app.param.solarCellEff(faceID)*100;                      % Cell efficiency
    app.EffecAreaRatioSpinnerZplus.Value = app.param.effectiveAreaRatio(faceID)*100;                % Packing factor
    app.AbsorptivitySpinnerZplus.Value = app.param.absorptivity(faceID);                            % Absorptivity
    app.EmissivitySpinnerZplus.Value = app.param.emissivity(faceID);                                % Emissivity
    app.RatioEditFieldZplus.Value = app.param.absorptivity(faceID)/app.param.emissivity(faceID);    % Ratio abs/emiss
    app.CpSpinnerZplus.Value = app.param.facesCp(faceID);                                           % Face Cp  
    app.HingeDropDownZplus.Value = app.param.deployable.hinge(faceID);                              % Hinge 
    app.SizeDropDownZplus.Value = app.param.deployable.size(faceID);                                % Size
    app.FlipPanelCheckBoxZplus.Value = app.param.deployable.flipFixedPanel(faceID);                 % Flip panel
    app.CellEfficiencySpinnerZplusDeployable.Value = app.param.deployable.solarCellEff(faceID)*100; % Cell Efficiency (Deployable)
    app.EffecAreaRatioSpinnerZplusDeployable.Value = app.param.deployable.effectiveAreaRatio(faceID)*100;% Packing factor (Deployable)
    switch app.param.deployable.type(faceID)                                                        % Deployable type
        case "None"
            app.PanelTypeButtonGroupZplus.SelectedObject = app.NoneButtonZplus;
            app.NoneButtonZplus.Value = 1;
        case "FiYed"
            app.PanelTypeButtonGroupZplus.SelectedObject = app.FiYedButtonZplus;
            app.FiYedButtonZplus.Value = 1;
        case "Tracking"
            app.PanelTypeButtonGroupZplus.SelectedObject = app.TrackingButtonZplus;
            app.TrackingButtonZplus.Value = 1;
    end
    if(app.SurfaceTypeDropDownZplus.Value == "Solar Panel")
        app.CellEfficiencySpinnerZplus.Visible = "on";
        app.CellEfficiencySpinnerLabelZplus.Visible = "on";
        app.EffecAreaRatioSpinnerZplus.Visible = "on";
        app.EffecAreaRatioSpinnerLabelZplus.Visible = "on";
    else
        app.CellEfficiencySpinnerZplus.Visible = "off";
        app.CellEfficiencySpinnerLabelZplus.Visible = "off";
        app.EffecAreaRatioSpinnerZplus.Visible = "off";
        app.EffecAreaRatioSpinnerLabelZplus.Visible = "off";
    end
    if(app.SurfaceTypeDropDownZplus.Value == "Custom")
        app.EmissivitySpinnerZplus.Enable = "on";
        app.AbsorptivitySpinnerZplus.Enable = "on";
        app.CpSpinnerZplus.Enable = "on";
    else
        app.EmissivitySpinnerZplus.Enable = "off";
        app.AbsorptivitySpinnerZplus.Enable = "off";
        app.CpSpinnerZplus.Enable = "off";
    end
    if(app.PanelTypeButtonGroupZplus.SelectedObject == "FiYed")
        app.FlipPanelCheckBoxZplus.Enable = "on";
    else
        app.FlipPanelCheckBoxZplus.Enable = "off";
    end
    if((app.PanelTypeButtonGroupZplus.SelectedObject == "FiYed") || (app.PanelTypeButtonGroupZplus.SelectedObject == "Tracking"))
        app.HingeDropDownZplus.Enable = "on";
        app.SizeDropDownZplus.Enable = "on";
        app.CellEfficiencySpinnerZplusDeployable.Enable = "on";
        app.EffecAreaRatioSpinnerZplusDeployable.Enable = "on";
    else
        app.HingeDropDownZplus.Enable = "off";
        app.SizeDropDownZplus.Enable = "off";
        app.CellEfficiencySpinnerZplusDeployable.Enable = "off";
        app.EffecAreaRatioSpinnerZplusDeployable.Enable = "off";
        app.FlipPanelCheckBoxZplus.Enable = "off";
    end

    %% ### Z- Face
    faceID = 6;
    app.SurfaceTypeDropDownZminus.Value = app.param.facesMaterial(faceID);                           % Surface type
    app.CellEfficiencySpinnerZminus.Value = app.param.solarCellEff(faceID)*100;                      % Cell efficiency
    app.EffecAreaRatioSpinnerZminus.Value = app.param.effectiveAreaRatio(faceID)*100;                % Packing factor
    app.AbsorptivitySpinnerZminus.Value = app.param.absorptivity(faceID);                            % Absorptivity
    app.EmissivitySpinnerZminus.Value = app.param.emissivity(faceID);                                % Emissivity
    app.RatioEditFieldZminus.Value = app.param.absorptivity(faceID)/app.param.emissivity(faceID);    % Ratio abs/emiss
    app.CpSpinnerZminus.Value = app.param.facesCp(faceID);                                           % Face Cp  
    app.HingeDropDownZminus.Value = app.param.deployable.hinge(faceID);                              % Hinge 
    app.SizeDropDownZminus.Value = app.param.deployable.size(faceID);                                % Size
    app.FlipPanelCheckBoxZminus.Value = app.param.deployable.flipFixedPanel(faceID);                 % Flip panel
    app.CellEfficiencySpinnerZminusDeployable.Value = app.param.deployable.solarCellEff(faceID)*100; % Cell Efficiency (Deployable)
    app.EffecAreaRatioSpinnerZminusDeployable.Value = app.param.deployable.effectiveAreaRatio(faceID)*100;% Packing factor (Deployable)
    switch app.param.deployable.type(faceID)                                                        % Deployable type
        case "None"
            app.PanelTypeButtonGroupZminus.SelectedObject = app.NoneButtonZminus;
            app.NoneButtonZminus.Value = 1;
        case "FiYed"
            app.PanelTypeButtonGroupZminus.SelectedObject = app.FiYedButtonZminus;
            app.FiYedButtonZminus.Value = 1;
        case "Tracking"
            app.PanelTypeButtonGroupZminus.SelectedObject = app.TrackingButtonZminus;
            app.TrackingButtonZminus.Value = 1;
    end
    if(app.SurfaceTypeDropDownZminus.Value == "Solar Panel")
        app.CellEfficiencySpinnerZminus.Visible = "on";
        app.CellEfficiencySpinnerLabelZminus.Visible = "on";
        app.EffecAreaRatioSpinnerZminus.Visible = "on";
        app.EffecAreaRatioSpinnerLabelZminus.Visible = "on";
    else
        app.CellEfficiencySpinnerZminus.Visible = "off";
        app.CellEfficiencySpinnerLabelZminus.Visible = "off";
        app.EffecAreaRatioSpinnerZminus.Visible = "off";
        app.EffecAreaRatioSpinnerLabelZminus.Visible = "off";
    end
    if(app.SurfaceTypeDropDownZminus.Value == "Custom")
        app.EmissivitySpinnerZminus.Enable = "on";
        app.AbsorptivitySpinnerZminus.Enable = "on";
        app.CpSpinnerZminus.Enable = "on";
    else
        app.EmissivitySpinnerZminus.Enable = "off";
        app.AbsorptivitySpinnerZminus.Enable = "off";
        app.CpSpinnerZminus.Enable = "off";
    end
    if(app.PanelTypeButtonGroupZminus.SelectedObject == "FiYed")
        app.FlipPanelCheckBoxZminus.Enable = "on";
    else
        app.FlipPanelCheckBoxZminus.Enable = "off";
    end
    if((app.PanelTypeButtonGroupZminus.SelectedObject == "FiYed") || (app.PanelTypeButtonGroupZminus.SelectedObject == "Tracking"))
        app.HingeDropDownZminus.Enable = "on";
        app.SizeDropDownZminus.Enable = "on";
        app.CellEfficiencySpinnerZminusDeployable.Enable = "on";
        app.EffecAreaRatioSpinnerZminusDeployable.Enable = "on";
    else
        app.HingeDropDownZminus.Enable = "off";
        app.SizeDropDownZminus.Enable = "off";
        app.CellEfficiencySpinnerZminusDeployable.Enable = "off";
        app.EffecAreaRatioSpinnerZminusDeployable.Enable = "off";
        app.FlipPanelCheckBoxZminus.Enable = "off";
    end

    % ### Update Satellite View
    plotCubesatView(app.param, app.UIAxesCubeSat);

    %% ### Power Tab ###
    app.ProfileDropDown.Value = app.param.pwrDissipationProfile;                % Dissipation profile type
    app.PwrDayWSpinner.Value = app.param.pwrDissipationDaySimpleModel;                 % Pwr Day
    app.PwrNightWSpinner.Value = app.param.pwrDissipationNightSimpleModel;      % Pwr Night
    switch app.param.pwrDissipationProfile 
        case "Constant"
            app.PwrNightWSpinner.Enable = "off";
        case "Day/Night"
            app.PwrNightWSpinner.Enable = "on";
    end
    % ### Update Dissipation Profile Plot
    plotPwrConsumption(app.UIAxesPwrConsumption, app.param.orb.prop.time_Epoch, app.param.orb.prop.sunMagnitude, app.param.pwrDissipationDaySimpleModel, app.param.pwrDissipationNightSimpleModel, app.param.pwrDissipationProfile)
    app.ComputePowerStatusLabel.Text = '';

    %% ### Thermal Tab ###
    app.ComputeThermalStatusLabel.Text = '';
end


function [thermal,flag] = computeThermal(param)
%COMPUTETHERMAL performs the thermal computation preparing and calling the
%Simscape Thermal Model


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
    tic
    out = sim(thermalModel.sim.in);
    timeElapsedProcessing = toc;
%     % Try fixed-step
%     % Solver setttings
% %     1. ode23t
% % 2. ode15s
% % 3. daessc
% % 4. ode14x
% % 5. ode1be
%     thermalModel.sim.in = thermalModel.sim.in.setModelParameter(...
%         SolverType="Variable-step",...
%         SolverName="ode23t",...
%         RelTol="1e-4",...
%         AbsTol="1e-5",...
%         StopTime=string(param.orb.durationHours*3600));
%     tic;
%     out2 = sim(thermalModel.sim.in);
%     disp("Fixed-step time:")
%     disp(toc)
% 
%     close all;
%     figure;
%     subplot(6,1,1)
%     plot(out.XplusTemp); hold on; plot(out2.XplusTemp);
%     subplot(6,1,2)
%     plot(out.XminusTemp); hold on; plot(out2.XminusTemp); 
%         subplot(6,1,3)
%     plot(out.YplusTemp); hold on; plot(out2.YplusTemp);
%     subplot(6,1,4)
%     plot(out.YminusTemp); hold on; plot(out2.YminusTemp); 
%         subplot(6,1,5)
%     plot(out.ZplusTemp); hold on; plot(out2.ZplusTemp);
%     subplot(6,1,6)
%     plot(out.ZminusTemp); hold on; plot(out2.ZminusTemp); 


    % Return results
    flag = 1;
    thermal.faceIR = param.thermal.faceIR;
    thermal.faceSolar = param.thermal.faceSolar;
    thermal.faceAlbedo = param.thermal.faceAlbedo;
    thermal.out = out;
    thermal.timeElapsedProcessing = timeElapsedProcessing;
end













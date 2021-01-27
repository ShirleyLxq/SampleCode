function [E_population_V_matrix, E_populationSpikesMatrix]= updateVoltageAndSpike_SingleInput(timeStepInd, SInputForEachE_neuron,...
    E_population_V_matrix, E_populationSpikesMatrix, J_EX, K, tau, delta_t, V_threshold)

sumOf_X_SInputForEachE_neuron=sum(SInputForEachE_neuron, 2); % a column vector
E_population_V_vector =E_population_V_matrix(:, timeStepInd);
w=J_EX/sqrt(K);
gradientForEachE_neuron=-E_population_V_vector/tau+...
    w*sumOf_X_SInputForEachE_neuron; % a column vector
newVForEachE_neuron=E_population_V_vector+delta_t*gradientForEachE_neuron; % a column vector
% if voltage exceeds threshold reset
newVHasExceededV_thresholdForEachE_neuron= (newVForEachE_neuron>V_threshold);
newVForEachE_neuron(newVHasExceededV_thresholdForEachE_neuron)= 0; % a column vector
E_population_V_matrix(:,timeStepInd+1)=newVForEachE_neuron;
E_populationSpikesMatrix(newVHasExceededV_thresholdForEachE_neuron, timeStepInd+1)=1/delta_t;


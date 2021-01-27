function [X_population_V_matrix, E_population_V_matrix, I_population_V_matrix,...
    X_populationSpikesMatrix, E_populationSpikesMatrix, I_populationSpikesMatrix]= ...
    updateVoltageAndSpike_MultipleInputs(N,K,timeStepInd,rate_x,...
    X_population_V_matrix, E_population_V_matrix, I_population_V_matrix,...
    X_populationSpikesMatrix, E_populationSpikesMatrix, I_populationSpikesMatrix,...
    indOfTheInputFromXToEachE, indOfTheInputFromEToEachE, indOfTheInputFromIToEachE,...
    indOfTheInputFromXToEachI, indOfTheInputFromEToEachI, indOfTheInputFromIToEachI)
    %indOfTheInputFromXToEachX, indOfTheInputFromEToEachX, indOfTheInputFromIToEachX,...
J_EE=1; J_IE=1; J_EI=-2; J_II=-1.8; J_EX=1; J_IX=0.8;
delta_t=0.1*10^(-3); % 0.1ms
tau=20*10^(-3); % 20ms
V_threshold=1;

spikeInputFromX=X_populationSpikesMatrix(:, timeStepInd);
spikeInputFromI=I_populationSpikesMatrix(:, timeStepInd);
spikeInputFromE=E_populationSpikesMatrix(:, timeStepInd);
%% 1. Update population E at time step timeStepInd+1
spikeInputFromXToEachE_neuron = spikeInputFromX(indOfTheInputFromXToEachE);
sumOfSpikeInputFromXToEachE_neuron=sum(spikeInputFromXToEachE_neuron, 2); % a column vector

spikeInputFromIToEachE_neuron = spikeInputFromI(indOfTheInputFromIToEachE);
sumOfSpikeInputFromIToEachE_neuron=sum(spikeInputFromIToEachE_neuron, 2); % a column vector

spikeInputFromEToEachE_neuron = spikeInputFromE(indOfTheInputFromEToEachE);
sumOfSpikeInputFromEToEachE_neuron=sum(spikeInputFromEToEachE_neuron, 2); % a column vector

E_population_V_vector =E_population_V_matrix(:, timeStepInd);

tempForEachE=J_EX*sumOfSpikeInputFromXToEachE_neuron+...
    J_EI*sumOfSpikeInputFromIToEachE_neuron+...
    J_EE*sumOfSpikeInputFromEToEachE_neuron;
gradientForEachE_neuron=-E_population_V_vector/tau+1/sqrt(K)*tempForEachE; % a column vector
newVForEachE_neuron=E_population_V_vector+delta_t*gradientForEachE_neuron; % a column vector
% if voltage exceeds threshold reset
newVHasExceededV_thresholdForEachE_neuron= (newVForEachE_neuron>V_threshold);
newVForEachE_neuron(newVHasExceededV_thresholdForEachE_neuron)= 0; % a column vector
E_population_V_matrix(:,timeStepInd+1)=newVForEachE_neuron;
E_populationSpikesMatrix(newVHasExceededV_thresholdForEachE_neuron, timeStepInd+1)=1/delta_t;

%% 2. Update population I at time step timeStepInd+1
spikeInputFromXToEachI_neuron = spikeInputFromX(indOfTheInputFromXToEachI);
sumOfSpikeInputFromXToEachI_neuron=sum(spikeInputFromXToEachI_neuron, 2); % a column vector

spikeInputFromIToEachI_neuron = spikeInputFromI(indOfTheInputFromIToEachI);
sumOfSpikeInputFromIToEachI_neuron=sum(spikeInputFromIToEachI_neuron, 2); % a column vector

spikeInputFromEToEachI_neuron = spikeInputFromE(indOfTheInputFromEToEachI);
sumOfSpikeInputFromEToEachI_neuron=sum(spikeInputFromEToEachI_neuron, 2); % a column vector

I_population_V_vector =I_population_V_matrix(:, timeStepInd);

tempForEachI=J_IX*sumOfSpikeInputFromXToEachI_neuron+...
    J_II*sumOfSpikeInputFromIToEachI_neuron+...
    J_IE*sumOfSpikeInputFromEToEachI_neuron;
gradientForEachI_neuron=-I_population_V_vector/tau+1/sqrt(K)*tempForEachI; % a column vector
newVForEachI_neuron=I_population_V_vector+delta_t*gradientForEachI_neuron; % a column vector
% if voltage exceeds threshold reset
newVHasExceededV_thresholdForEachI_neuron= (newVForEachI_neuron>V_threshold);
newVForEachI_neuron(newVHasExceededV_thresholdForEachI_neuron)= 0; % a column vector
I_population_V_matrix(:, timeStepInd+1)=newVForEachI_neuron;
I_populationSpikesMatrix(newVHasExceededV_thresholdForEachI_neuron, timeStepInd+1)=1/delta_t;

%% 3. Update population X at time step timeStepInd+1
% a new series of poisson events
X_population_SInputVector=fullNetwork_generateSInputVector(N, rate_x, delta_t); % a column vector
X_populationSpikesMatrix(:, timeStepInd+1)=X_population_SInputVector;
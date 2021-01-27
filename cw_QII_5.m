%% Question II.5 full network
function [rate_E,rate_I]=cw_QII_5(rate_x)
% [rate_E,rate_I]=cw_QII_5(isRunningRemotely, rate_x)
delta_t=0.1*10^(-3); % 0.1ms
t=2; %2s
%rate_x=10; % 10Hz
numTimeSteps=(t-0)/delta_t+1;
N=100; % num of neurons 1000
K=99; % 99
times=0:delta_t:t;
timeIndex=0.1/delta_t+1; % 100ms
assert(length(times)==numTimeSteps);

tau=20*10^(-3); % 20ms
V_threshold=1;
J_EE=1; J_IE=1; J_EI=-2; J_II=-1.8; J_EX=1; J_IX=0.8;

X_spikeTimesCellArray = cell(N,1);
E_spikeTimesCellArray = cell(N,1);
I_spikeTimesCellArray = cell(N,1);

X_population_V_matrix=zeros(N, numTimeSteps);
E_population_V_matrix=zeros(N, numTimeSteps);
I_population_V_matrix=zeros(N, numTimeSteps);

X_populationSpikesMatrix=zeros(N, numTimeSteps);
E_populationSpikesMatrix=zeros(N, numTimeSteps);
I_populationSpikesMatrix=zeros(N, numTimeSteps);

%[indOfTheInputFromXToEachX, indOfTheInputFromEToEachX, indOfTheInputFromIToEachX,...
    [indOfTheInputFromXToEachE, indOfTheInputFromEToEachE, indOfTheInputFromIToEachE,...
    indOfTheInputFromXToEachI, indOfTheInputFromEToEachI, indOfTheInputFromIToEachI]= fullNetwork_generateConnectivityMatrices(N,K);

%% A. INITIALISATION
% at the initial time step only the X pop elicits spikes while the E and I populations remain silent
timeStepInd=1;
X_population_SInputVector=fullNetwork_generateSInputVector(N, rate_x, delta_t); % a column vector
X_populationSpikesMatrix(:, timeStepInd)=X_population_SInputVector;
%% Update time step 2 based on 1
% 1. Update population E at time step 2
SInputForEachE_neuron = X_population_SInputVector(indOfTheInputFromXToEachE);
[E_population_V_matrix, E_populationSpikesMatrix]= updateVoltageAndSpike_SingleInput(timeStepInd, SInputForEachE_neuron,...
    E_population_V_matrix, E_populationSpikesMatrix, J_EX, K, tau, delta_t, V_threshold);

% sumOf_X_SInputForEachE_neuron=sum(SInputForEachE_neuron, 2); % a column vector
% E_population_V_vector =E_population_V_matrix(:, timeStepInd);
% w=J_EX/sqrt(K);
% gradientForEachE_neuron=-E_population_V_vector/tau+...
%     w*sumOf_X_SInputForEachE_neuron; % a column vector
% newVForEachE_neuron=E_population_V_vector+delta_t*gradientForEachE_neuron; % a column vector
% % if voltage exceeds threshold reset
% newVHasExceededV_thresholdForEachE_neuron= (newVForEachE_neuron>V_threshold);
% newVForEachE_neuron(newVHasExceededV_thresholdForEachE_neuron)= 0; % a column vector
% E_population_V_matrix(:,timeStepInd+1)=newVForEachE_neuron;
% E_populationSpikesMatrix(newVHasExceededV_thresholdForEachE_neuron, timeStepInd+1)=1/delta_t;

% 2. Update population I at time step 2
% presyn_X_neuronIndicesForEachI_neuron =fullNetwork_generateConnectivityMatrix(N,K);
SInputForEachI_neuron = X_population_SInputVector(indOfTheInputFromXToEachI);
[I_population_V_matrix, I_populationSpikesMatrix]= updateVoltageAndSpike_SingleInput(timeStepInd, SInputForEachI_neuron,...
    I_population_V_matrix, I_populationSpikesMatrix, J_IX, K, tau, delta_t, V_threshold);

% 3. Update population X at time step 2
% a new series of poisson events
X_population_SInputVector=fullNetwork_generateSInputVector(N, rate_x, delta_t); % a column vector
X_populationSpikesMatrix(:, timeStepInd+1)=X_population_SInputVector;

%% B. UPDATE
% now three sources of input rather than one
for timeStepInd=2:numTimeSteps-1
    if (timeStepInd==2 || mod(timeStepInd,1000)==0)
        fprintf('timeStep=%d, time=%f...\n',timeStepInd, times(timeStepInd));
    end
   % 1. Update population E at time step 2+1
 [X_population_V_matrix, E_population_V_matrix, I_population_V_matrix,...
    X_populationSpikesMatrix, E_populationSpikesMatrix, I_populationSpikesMatrix]= ...
    updateVoltageAndSpike_MultipleInputs(N,K,timeStepInd,rate_x,...
    X_population_V_matrix, E_population_V_matrix, I_population_V_matrix,...
    X_populationSpikesMatrix, E_populationSpikesMatrix, I_populationSpikesMatrix,...
    indOfTheInputFromXToEachE, indOfTheInputFromEToEachE, indOfTheInputFromIToEachE,...
    indOfTheInputFromXToEachI, indOfTheInputFromEToEachI, indOfTheInputFromIToEachI); 
%     if(mod(timeStepInd, 20)==0)
%         figure(1);plot(E_population_V_matrix(:, timeStepInd));
%         figure(2);plot(I_population_V_matrix(:, timeStepInd));
%         figure(3);plot(E_populationSpikesMatrix(:, timeStepInd));
%         figure(4);plot(I_populationSpikesMatrix(:, timeStepInd));
%     end
end

%% compute the firing rates
averageNumSpikesPerE_neuron=sum(sum(E_populationSpikesMatrix~=0))/N;
averageNumSpikesPerI_neuron=sum(sum(I_populationSpikesMatrix~=0))/N;
% averageNumSpikesPerX_neuron=sum(sum(X_populationSpikesMatrix~=0))/N;
% expectedNumSpikes=rate_x*t;
rate_E=averageNumSpikesPerE_neuron/t;
rate_I=averageNumSpikesPerI_neuron/t;

% exclude the first column of the spikes matrix which corresponds to time
% zero
for i=1:N
    X_spikeTimesCellArray{i}=times(X_populationSpikesMatrix(i, :)~=0);
    E_spikeTimesCellArray{i}=times(E_populationSpikesMatrix(i, :)~=0);
    I_spikeTimesCellArray{i}=times(I_populationSpikesMatrix(i, :)~=0);
end
%% SPIKE PLOTS
XSpikesFigHandle=figure;
plotSpikeRaster(X_spikeTimesCellArray,'PlotType','vertline');
title({'The raster plot for the Poisson neurons in population X',...
    sprintf('Average firing rate $r_X=$%.3f',rate_x)},'Interpreter','latex','FontSize',14);
ylabel('Neuron index $1\leq i\leq N$','Interpreter','latex','FontSize',14);
xlabel('Time $0\leq k\delta_t\leq 2$ (s)', 'Interpreter','latex','FontSize',14);
XSpikesFigName=sprintf('XSpikesRaster,rx=%d-%s', rate_x, datestr(now));
saveas(XSpikesFigHandle, XSpikesFigName, 'jpg');
saveas(XSpikesFigHandle, XSpikesFigName, 'fig');

ESpikesFigHandle=figure;
plotSpikeRaster(E_spikeTimesCellArray,'PlotType','vertline');
title({'The raster plot for the excitatory neurons in population E',...
    sprintf('Average firing rate $r_E=$%.3f',rate_E)},'Interpreter','latex','FontSize',14);
ylabel('Neuron index $1\leq i\leq N$','Interpreter','latex','FontSize',14);
xlabel('Time $0\leq k\delta_t\leq 2$ (s)', 'Interpreter','latex','FontSize',14);
ESpikesFigName=sprintf('ESpikesRaster,rx=%d-%s', rate_x, datestr(now));
saveas(ESpikesFigHandle, ESpikesFigName, 'jpg');
saveas(ESpikesFigHandle, ESpikesFigName, 'fig');

ISpikesFigHandle=figure;
plotSpikeRaster(I_spikeTimesCellArray,'PlotType','vertline');
title({'The raster plot for the inhibitory neurons in population I',...
    sprintf('Average firing rate $r_I=$%.3f',rate_I)},'Interpreter','latex','FontSize',14);
ylabel('Neuron index $1\leq i\leq N$','Interpreter','latex','FontSize',14);
xlabel('Time $0\leq k\delta_t\leq 2$ (s)', 'Interpreter','latex','FontSize',14);
ISpikesFigName=sprintf('ISpikesRaster,rx=%d-%s', rate_x, datestr(now));
saveas(ISpikesFigHandle, ISpikesFigName, 'jpg');
saveas(ISpikesFigHandle, ISpikesFigName, 'fig');

%% MEMBRANE VOLTAGE PLOTS
selectedPlotIndices=[1:1000:20001];
numSelectedPlots=length(selectedPlotIndices);
colourList=colormap(jet(numSelectedPlots));

EVoltageFigureHandle=figure;
for ind=1:numSelectedPlots
    plotInd=selectedPlotIndices(ind);
%     if(ind==numSelectedPlots)
%         plot(1:N, E_population_V_matrix(:, plotInd), 'Color' ,colourList(ind, :), 'LineWidth',2); hold on;
%     else
        plot(1:N, E_population_V_matrix(:, plotInd), 'Color' ,colourList(ind, :)); hold on;
%     end
end
colormap(jet);
cb=colorbar; %('Location', 'SouthOutside', 'Limits', [1 numSelectedPlots/(numSelectedPlots+10)]);
caxis([0 2]);
set(cb, 'YTick',[0,2]); % 0s and 2s
title('The membrane potential of the excitatory neurons in population E','Interpreter','latex','FontSize',14);
ylabel('$\hat{V}(k)$','Interpreter','latex','FontSize',14);
xlabel('Neuron index $1\leq i\leq N$', 'Interpreter','latex','FontSize',14);
EVoltageFigName=sprintf('EVoltage-%s', datestr(now));
saveas(EVoltageFigureHandle, EVoltageFigName, 'jpg');
saveas(EVoltageFigureHandle, EVoltageFigName, 'fig');

IVoltageFigureHandle=figure;
for ind=1:numSelectedPlots
    plotInd=selectedPlotIndices(ind);
%     if(ind==numSelectedPlots)
%         plot(1:N, E_population_V_matrix(:, plotInd), 'Color' ,colourList(ind, :), 'LineWidth',2); hold on;
%     else
        plot(1:N, I_population_V_matrix(:, plotInd), 'Color' ,colourList(ind, :)); hold on;
%     end
end
colormap(jet);
cb=colorbar; %('Location', 'SouthOutside', 'Limits', [1 numSelectedPlots/(numSelectedPlots+10)]);
caxis([0 2]);
set(cb, 'YTick',[0,2]); % 0s and 2s
title('The membrane potential of the inhibitory neurons in population I','Interpreter','latex','FontSize',14);
ylabel('$\hat{V}(k)$','Interpreter','latex','FontSize',14);
xlabel('Neuron index $1\leq i\leq N$', 'Interpreter','latex','FontSize',14);
IVoltageFigName=sprintf('IVoltage-%s', datestr(now));
saveas(IVoltageFigureHandle, IVoltageFigName, 'jpg');
saveas(IVoltageFigureHandle, IVoltageFigName, 'fig');

% if(isRunningRemotely)
%     exit;
% end
%% Question II.2 and Question II.3 (a)
delta_t=0.1*10^(-3); % 0.1ms
t=2; %2s
rate_x=10; % 10Hz
numTimeSteps=(t-0)/delta_t+1;
N=1000; % num of neurons
times=0:delta_t:t;
assert(length(times)==numTimeSteps);

% w=0.9;
w=1;
tau=20*10^(-3); % 20ms
V_threshold=1;

%% input spike train
SInputArray=1/delta_t*(rand(1,numTimeSteps)<rate_x*delta_t); % uniform distribution between 0 and 1
inputSpikeIsPresentArray=(SInputArray~=0);
inputSpikeTimes=times(inputSpikeIsPresentArray);
inputSpikeTimesCellArray{1}=inputSpikeTimes;

%% post synaptic voltage
SOutputArray=zeros(1, numTimeSteps); % Si array
V_array=zeros(1, numTimeSteps); % output voltage array
V=0;
V_array(1)=V;
for k=2:numTimeSteps
gradient=-V_array(k-1)/tau+w*SInputArray(k-1);
newV=V_array(k-1)+delta_t*gradient;
V_array(k)=newV;
if(newV>V_threshold)
    SOutputArray(k)=1/delta_t;
    % if exceeding the threshold reset V(k)
    %newV=0;
    %V_array(k)=newV;
end
end
outputSpikeIsPresentArray = (SOutputArray~=0);
outputSpikeTimes=times(outputSpikeIsPresentArray);
outputSpikeTimesCellArray{1}=outputSpikeTimes;

%%
figure;
subplot(4,1,1);
plotSpikeRaster(inputSpikeTimesCellArray,'PlotType','vertline');
title('The single input Poisson spike train to the LIF neuron $\tilde{S}_j(k)$','Interpreter','latex','FontSize',14);
ylabel('$\tilde{S}_j(k)$','Interpreter','latex','FontSize',12);
xlim([0,t]);

subplot(4,1,[2,3]);
plot(times, V_array);
title('The membrane potential of the LIF neuron $\tilde{V}(k)$','Interpreter','latex','FontSize',14);
ylabel('$\tilde{V}(k)$','Interpreter','latex','FontSize',12);
xlim([0,t]);

subplot(4,1,4);
plotSpikeRaster(outputSpikeTimesCellArray,'PlotType','vertline');
title('The output spike train of the LIF neuron $\tilde{S}_i(k)$','Interpreter','latex','FontSize',14);
ylabel('$\tilde{S}_i(k)$','Interpreter','latex','FontSize',12);
xlim([0,t]);
xlabel('Time $0\leq k\delta_t\leq 2$ (s)', 'Interpreter','latex','FontSize',14);
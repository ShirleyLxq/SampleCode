%% Question II.4
delta_t=0.1*10^(-3); % 0.1ms
t=10; %2s
rate_x=10; % 10Hz
numTimeSteps=(t-0)/delta_t+1;
N=1000; % num of neurons
K=100;
times=0:delta_t:t;
timeIndex=0.1/delta_t+1; % 100ms
assert(length(times)==numTimeSteps);

w=1.55;%1;
tau=20*10^(-3); % 20ms
V_threshold=1;

%% excitatory input spike trains
excitatorySInputMatrix=1/delta_t*(rand(K,numTimeSteps)<rate_x*delta_t); % uniform distribution between 0 and 1
sumOfExcitatorySInputOverCells=sum(excitatorySInputMatrix, 1); % sum over K presyn partners

%% inhibitory input spike trains
inhibitorySInputMatrix=1/delta_t*(rand(K,numTimeSteps)<rate_x*delta_t); % uniform distribution between 0 and 1
sumOfInhibitorySInputOverCells=sum(inhibitorySInputMatrix, 1); % sum over K presyn partners

%% post synaptic voltage
SOutputArray=zeros(1, numTimeSteps); % Si array
V_array=zeros(1, numTimeSteps); % output voltage array
V=0;
V_array(1)=V;
for k=2:numTimeSteps
    gradient=-V_array(k-1)/tau+w/sqrt(K)*sumOfExcitatorySInputOverCells(k-1)-w/sqrt(K)*sumOfInhibitorySInputOverCells(k-1);
    newV=V_array(k-1)+delta_t*gradient;
    V_array(k)=newV;
    if(newV>V_threshold)
        SOutputArray(k)=1/delta_t;
        % if exceeding the threshold reset V(k)
        newV=0;
        V_array(k)=newV;
    end
end
outputSpikeIsPresentArray = (SOutputArray~=0);
outputSpikeTimes=times(outputSpikeIsPresentArray);
outputSpikeTimesCellArray{1}=outputSpikeTimes;
%% empirical average firing rate
averageNumSpikesPerSecond=sum(outputSpikeIsPresentArray)/t % this has the unit of Hz
%%
figure;
subplot(3,1,[1,2]);
plot(times, V_array); hold on;
title('The membrane potential of the LIF neuron $\tilde{V}(k)$','Interpreter','latex','FontSize',14);
ylabel('$\tilde{V}(k)$','Interpreter','latex','FontSize',14);
%xlabel('Time $k$ (s)', 'Interpreter','latex','FontSize',14);
xlim([0,t]);

subplot(3,1,3);
plotSpikeRaster(outputSpikeTimesCellArray,'PlotType','vertline');
title('The output spike train of the LIF neuron $\tilde{S}_i(k)$','Interpreter','latex','FontSize',14);
ylabel('$\tilde{S}_i(k)$','Interpreter','latex','FontSize',12);
xlabel({'Time $0\leq k\delta_t\leq 10$ (s)', sprintf('w=%.3f, Average firing rate=%.3f Hz', w, averageNumSpikesPerSecond)}, 'Interpreter','latex','FontSize',14);
xlim([0,t]);

%% compute the fano factor
window=100*10^(-3); % 100ms
windowWidth=window/delta_t;
timeStepIndMatrix=[];
for startOfWindowInd = 1:ceil(windowWidth/2):(numTimeSteps-windowWidth+1)
    timeStepIndMatrix=[timeStepIndMatrix; startOfWindowInd:startOfWindowInd+windowWidth-1];
end
outputSpikeIsPresentMatrix=outputSpikeIsPresentArray(timeStepIndMatrix);
numSpikesInEachWindow=sum(outputSpikeIsPresentMatrix, 2);
mean_NumSpikesInEachWindow=mean(numSpikesInEachWindow);
variance_NumSpikesInEachWindow=std(numSpikesInEachWindow).^2;
fanoFactor=variance_NumSpikesInEachWindow/mean_NumSpikesInEachWindow
%% Question II.1 generating poisson spike trains
delta_t=0.1*10^(-3); % 0.1ms
t=2; %2s
rate_x=10; % 10Hz
numTimeSteps=(t-0)/delta_t+1;
N=1000; % num of neurons
times=0:delta_t:t;
assert(length(times)==numTimeSteps);
numSpikesArray=zeros(N,1);
spikeTimesCellArray = cell(N,1);

for i=1:N % for cell i
valueInEachBin=1/delta_t*(rand(1,numTimeSteps)<rate_x*delta_t); % uniform distribution between 0 and 1
numSpikesArray(i)=sum(valueInEachBin~=0);
spikeTimes=times(valueInEachBin~=0);
spikeTimesCellArray{i}=spikeTimes;
end
%%
figure;
plotSpikeRaster(spikeTimesCellArray,'PlotType','vertline');
title('The raster plot for the Poisson neurons in population X','Interpreter','latex','FontSize',14);
ylabel('Neuron index $1\leq i\leq N$','Interpreter','latex','FontSize',14);
xlabel('Time $0\leq k\delta_t\leq 2$ (s)', 'Interpreter','latex','FontSize',14);
%% average num of spikes
averageNumSpikes=sum(numSpikesArray)/N
expectedNumSpikes=rate_x*t
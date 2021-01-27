%% Question II.3
%% (a)
delta_t=0.1*10^(-3); % 0.1ms
t=10; %2s
rate_x=10; % 10Hz
numTimeSteps=(t-0)/delta_t+1;
N=1000; % num of neurons
K=100;
times=0:delta_t:t;
timeIndex=0.1/delta_t+1; % 100ms
assert(length(times)==numTimeSteps);

w=5;%1;
tau=20*10^(-3); % 20ms
V_threshold=1;

%% input spike train
SInputMatrix=1/delta_t*(rand(K,numTimeSteps)<rate_x*delta_t); % uniform distribution between 0 and 1
sumOfSInputOverCells=sum(SInputMatrix, 1); % sum over K presyn partners

%% post synaptic voltage
V_array=zeros(1, numTimeSteps); % output voltage array
V=0;
V_array(1)=V;
for k=2:numTimeSteps
    gradient=-V_array(k-1)/tau+w/K*sumOfSInputOverCells(k-1);
    newV=V_array(k-1)+delta_t*gradient;
    V_array(k)=newV;
end

%%
figure;
plot(times, V_array); hold on;
title('The membrane potential of the LIF neuron $\tilde{V}(t)$','Interpreter','latex','FontSize',14);
ylabel('$\tilde{V}(t)$','Interpreter','latex','FontSize',14);
xlabel('Time t (s)', 'Interpreter','latex','FontSize',14);
xlim([0,t]);

%% (d)
voltageMean=mean(V_array(timeIndex+1:end));
%voltageVariance=std(V_array(timeIndex+1:end))^2;
plot(times, voltageMean*ones(size(times)), 'LineWidth',2);
legend({'Membrane voltage', sprintf('Mean=%f', voltageMean)}, 'Interpreter','latex','FontSize',14);

%% (c)
%% theoretical mu and sigma2
K_theoreticalArray=1:N;
mu_h=w*rate_x;
% sigma_h2=w^2*rate_x./K_theoreticalArray*(-rate_x+1/delta_t);
sigma_h2=w^2*rate_x./K_theoreticalArray;
mu=tau*mu_h;
sigma2=tau*sigma_h2/2;

%% empirical mu and sigma2
K_expArray=[1,10,100,1000];
numK =length(K_expArray);
voltageMeanArray=zeros(1, numK); % mean of the membrane voltage
voltageVarianceArray=zeros(1, numK);
for KInd=1:numK
    K=K_expArray(KInd);
    SInputMatrix=1/delta_t*(rand(K,numTimeSteps)<rate_x*delta_t); % uniform distribution between 0 and 1
    sumOfSInputOverCells=sum(SInputMatrix, 1); % sum over K presyn partners
    
    V_array=zeros(1, numTimeSteps); % output voltage array
    V=0;
    V_array(1)=V;
    for k=2:numTimeSteps
        gradient=-V_array(k-1)/tau+w/K*sumOfSInputOverCells(k-1);
        newV=V_array(k-1)+delta_t*gradient;
        V_array(k)=newV;
    end
    voltageMeanArray(KInd)=mean(V_array(timeIndex+1:end));
    voltageVarianceArray(KInd)=std(V_array(timeIndex+1:end))^2;
end
%%
figure;
plot(K_theoreticalArray, mu*ones(size(K_theoreticalArray)), 'LineWidth',2); hold on;
plot(K_expArray, voltageMeanArray, 'o--');
title('The mean of the membrane potential $\mu$','Interpreter','latex','FontSize',14);
ylabel('$\mu$','Interpreter','latex','FontSize',14);
xlabel('The number of presynaptic neurons $K$','Interpreter','latex','FontSize',14);
xlim([0,N]);
legend({'Theoretical', 'Simulated'}, 'Interpreter','latex','FontSize',14);

figure;
plot(K_theoreticalArray, sigma2, 'LineWidth',2); hold on;
plot(K_expArray, voltageVarianceArray,'o--');
title('The variance of the membrane potential $\sigma^2$','Interpreter','latex','FontSize',14);
ylabel('$\sigma^2$','Interpreter','latex','FontSize',14);
xlabel('The number of presynaptic neurons $K$','Interpreter','latex','FontSize',14);
xlim([0,N]);
legend({'Theoretical', 'Simulated'}, 'Interpreter','latex','FontSize',14);


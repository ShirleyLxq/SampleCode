%% (c)
% boxcar
gamma=1; epsilon=0.01; % 0.01;
N=201; % ntrials
delta_t=0.5;
timeStepsWithinTrial=0:delta_t:25; % time within a trial
numTimeStepsWithinTrial=length(timeStepsWithinTrial);
Tmem=12;
numFeatureDetectors=(12-0)/0.5+1;
firstPartOfS=zeros(20,25);
middlePartOfS=ones(25,25);
middlePartOfS=triu(middlePartOfS);
lastPartOfS=zeros(6,25);
S_compositeMatrix=[firstPartOfS; middlePartOfS; lastPartOfS];
assert(isequal(size(S_compositeMatrix), [numTimeStepsWithinTrial, numFeatureDetectors]));
% each column transposed is the 1-by-25 S
S_compositeMatrix_oneStepBackInTime=[zeros(1, numFeatureDetectors); S_compositeMatrix(1:end-1, :)];
reward_oneStepBackInTime=reward(timeStepsWithinTrial- delta_t);

weights=zeros(1, numFeatureDetectors); % row vector
valueFnMatrix=zeros(N, numTimeStepsWithinTrial); % store a row vector after each trial
TDMatrix=zeros(N, numTimeStepsWithinTrial);
TD_learningErrorMatrix=zeros(N, numTimeStepsWithinTrial);

for n=1:N % N trials
    %valueFnVector=weights*S_compositeMatrix'; % producing a row vector   
    for rowInd=1:numTimeStepsWithinTrial
        S_current=S_compositeMatrix(rowInd,:); % a row vector S
        valueFn=weights*S_current'; % for the current time
        S_future=S_compositeMatrix(min(rowInd+1, numTimeStepsWithinTrial),:); % a row vector S
        futureValueFn=weights*S_future'; % for one time step into the future
        t= timeStepsWithinTrial(rowInd);
        % update weights
        r_t=reward(t);
        TD=gamma*futureValueFn-valueFn;
        %TD_learningError=reward(t-delta_t)+TD;
        TD_update=epsilon*(r_t+TD);
        weights=weights+TD_update*S_current;
    end
    valueFnVector=weights*S_compositeMatrix';
    valueFnMatrix(n,:)=valueFnVector;
    TDVector=weights*S_compositeMatrix'-weights*S_compositeMatrix_oneStepBackInTime';
    TDMatrix(n,:)=TDVector;
    TD_learningErrorVector=reward_oneStepBackInTime+TDVector;
    TD_learningErrorMatrix(n,:)=TD_learningErrorVector;
    if(n==1)
        weights
    end
end
%% plot three graphs
figure(4);
subplot(3,1,1);
selectedPlotIndices=1:10:N;
numSelectedPlots=length(selectedPlotIndices);
colourList=colormap(jet(numSelectedPlots));
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, valueFnMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The value function','Interpreter','latex','FontSize',14);
%xlabel('Time t','Interpreter','latex','FontSize',12);
ylabel('$\hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,3]);
    
subplot(3,1,2);
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, TDMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The TD of the value','Interpreter','latex','FontSize',14);
%xlabel('Time t','Interpreter','latex','FontSize',12);
ylabel('$\Delta \hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,5]);

subplot(3,1,3);
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, TD_learningErrorMatrix(plotInd, :), 'Color' ,colourList(ind, :));  hold on;
end
title('The TD learning error','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
ylabel('$\delta(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,3]);
figure(4);
cb=colorbar; %('Location', 'SouthOutside', 'Limits', [1 numSelectedPlots/(numSelectedPlots+10)]);
caxis([1 numSelectedPlots]);
set(cb, 'YTick',[1,numSelectedPlots])

%% Plots for only the first and the last trials
figure(5);
subplot(3,1,1);
plot(timeStepsWithinTrial, valueFnMatrix(1, :), 'b'); hold on;
plot(timeStepsWithinTrial, valueFnMatrix(end, :), 'r');
title('The value function','Interpreter','latex','FontSize',14);
ylabel('$\hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);

subplot(3,1,2);
plot(timeStepsWithinTrial, TDMatrix(1, :), 'b'); hold on;
plot(timeStepsWithinTrial, TDMatrix(end, :), 'r');
title('The TD of the value','Interpreter','latex','FontSize',14);
ylabel('$\Delta \hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);

subplot(3,1,3);
plot(timeStepsWithinTrial, TD_learningErrorMatrix(1, :), 'b'); hold on;
plot(timeStepsWithinTrial, TD_learningErrorMatrix(end, :), 'r');
legend({'Example trial 1', 'Example trial 21'},'Interpreter','latex','FontSize',10);
title('The TD learning error','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
ylabel('$\delta(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);

%% (c) plot three graphs
figure(6);
subplot(3,1,1);
selectedPlotIndices=1:10:31;
numSelectedPlots=length(selectedPlotIndices);
colourList=colormap(jet(numSelectedPlots));
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, valueFnMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The value function','Interpreter','latex','FontSize',14);
%xlabel('Time t','Interpreter','latex','FontSize',12);
ylabel('$\hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,3]);
    
subplot(3,1,2);
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, TDMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The TD of the value','Interpreter','latex','FontSize',14);
%xlabel('Time t','Interpreter','latex','FontSize',12);
ylabel('$\Delta \hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,5]);

subplot(3,1,3);
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, TD_learningErrorMatrix(plotInd, :), 'Color' ,colourList(ind, :));  hold on;
end
title('The TD learning error','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
ylabel('$\delta(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,3]);
figure(6);
cb=colorbar; %('Location', 'SouthOutside', 'Limits', [1 numSelectedPlots/(numSelectedPlots+10)]);
caxis([1 numSelectedPlots]);
set(cb, 'YTick',[1,numSelectedPlots])


%%
figure(7);
subplot(3,1,1);
selectedPlotIndices=1:10;
numSelectedPlots=length(selectedPlotIndices);
colourList=colormap(jet(numSelectedPlots));
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, valueFnMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The value function','Interpreter','latex','FontSize',14);
%xlabel('Time t','Interpreter','latex','FontSize',12);
ylabel('$\hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,3]);
    
subplot(3,1,2);
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, TDMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The TD of the value','Interpreter','latex','FontSize',14);
%xlabel('Time t','Interpreter','latex','FontSize',12);
ylabel('$\Delta \hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,5]);

subplot(3,1,3);
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, TD_learningErrorMatrix(plotInd, :), 'Color' ,colourList(ind, :));  hold on;
end
title('The TD learning error','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
ylabel('$\delta(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,3]);
figure(7);
cb=colorbar; %('Location', 'SouthOutside', 'Limits', [1 numSelectedPlots/(numSelectedPlots+10)]);
caxis([1 numSelectedPlots]);
set(cb, 'YTick',[1,numSelectedPlots])

%% compare different learnig rates
figure(100);
if(epsilon==0.2)
subplot(2,1,1);
selectedPlotIndices=1:20;
numSelectedPlots=length(selectedPlotIndices);
colourList=colormap(jet(numSelectedPlots));
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, valueFnMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The value function ($\epsilon=0.2$)','Interpreter','latex','FontSize',14);
%xlabel('Time t','Interpreter','latex','FontSize',12);
ylabel('$\hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
cb=colorbar; %('Location', 'SouthOutside', 'Limits', [1 numSelectedPlots/(numSelectedPlots+10)]);
caxis([1 numSelectedPlots]);
set(cb, 'YTick',[1,numSelectedPlots]);

else
    assert(epsilon==0.01);
    subplot(2,1,2);
selectedPlotIndices=1:100;
numSelectedPlots=length(selectedPlotIndices);
colourList=colormap(jet(numSelectedPlots));
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, valueFnMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The value function ($\epsilon=0.01$)','Interpreter','latex','FontSize',14);
%xlabel('Time t','Interpreter','latex','FontSize',12);
ylabel('$\hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
cb=colorbar; %('Location', 'SouthOutside', 'Limits', [1 numSelectedPlots/(numSelectedPlots+10)]);
caxis([1 numSelectedPlots]);
set(cb, 'YTick',[1,numSelectedPlots]);
end

xlabel('Time t','Interpreter','latex','FontSize',14);
%%
function r_t=reward(t)
sigma=1;
r_t=1/2*exp( -(t-20).^2/(2*sigma^2) );
end
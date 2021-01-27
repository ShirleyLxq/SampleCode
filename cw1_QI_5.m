%% Partial reinforcement learning
% boxcar representation is used
%%
figure;
subplot(2,1,1);
plot([10, 10], [0, 1], 'k', 'linewidth',2); % stimulus
title('The stimulus','Interpreter','latex','FontSize',14);
%xlabel('Time t','Interpreter','latex','FontSize',12);
ylabel('$s(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
ylim([0,2]);

subplot(2,1,2);
rewardTimeSpan=20+[-5:0.01:5];
rewardMagnitude=1/2*exp(-(rewardTimeSpan-20).^2/(2));
plot(rewardTimeSpan, rewardMagnitude, 'r', 'LineWidth',2); hold on;
plot([0,25], [0,0], 'b', 'LineWidth',2);
legend({'Reward present', 'Reward absent'}, 'Interpreter','latex','FontSize',10);
title('The reward','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
ylabel('$r(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
ylim([0,2]);
%%
gamma=1; epsilon=0.01; % 0.01;
N=1000; % ntrials
p=0.25; % chance of p that reward is present; 1-p it is not
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

weights=zeros(1, numFeatureDetectors); % row vector
valueFnMatrix=zeros(N, numTimeStepsWithinTrial); % store a row vector after each trial
TDMatrix=zeros(N, numTimeStepsWithinTrial);
TD_learningErrorMatrix=zeros(N, numTimeStepsWithinTrial);
rewardIsPresentInEachTrial=zeros(N, 1);
 
for n=1:N % N trials
    % determine whether reward is presence or absence
    rewardIsPresent=(rand(1)<p); % uniform distribution between 0 and 1
    rewardIsPresentInEachTrial(n)=rewardIsPresent;
    if(rewardIsPresent)
        reward_oneStepBackInTime=reward(timeStepsWithinTrial- delta_t); 
    else
        reward_oneStepBackInTime=zeros(size(timeStepsWithinTrial));
    end
    for rowInd=1:numTimeStepsWithinTrial
        S_current=S_compositeMatrix(rowInd,:); % a row vector S
        valueFn=weights*S_current'; % for the current time
        S_future=S_compositeMatrix(min(rowInd+1, numTimeStepsWithinTrial), :); % a row vector S
        futureValueFn=weights*S_future'; % for one time step into the future
        t= timeStepsWithinTrial(rowInd);
        % update weights
        if(rewardIsPresent)
            r_t=reward(t);
        else
            r_t=0;
        end
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
end

%%
selectedPlotIndices=1:20;
selectedPlotIndices=(N-20+1):N;
plotEvolutionTripleGraphs(timeStepsWithinTrial, selectedPlotIndices, valueFnMatrix, TDMatrix, TD_learningErrorMatrix)

%% plot three graphs
totalNumTrials=100;
selectedTrials= [N-totalNumTrials+1: N]; % retain only the data for the last 100 trials
%selectedTrials= [1:totalNumTrials]; % retain only the data for the last 100 trials

rewardIsPresentInEachTrial=logical(rewardIsPresentInEachTrial(selectedTrials));
valueFnMatrix=valueFnMatrix(selectedTrials, :);
TDMatrix=TDMatrix(selectedTrials, :);
TD_learningErrorMatrix=TD_learningErrorMatrix(selectedTrials, :);
numRewardedTrials=sum(rewardIsPresentInEachTrial);
numUnrewardedTrials=totalNumTrials-numRewardedTrials;

%%
figure;
% the 3 variables averaged for the rewarded trials
try
averageValueFnVector_ForRewardedTrials=sum(valueFnMatrix(rewardIsPresentInEachTrial, :))/numRewardedTrials;
catch
end
averageTDVector_ForRewardedTrials=sum(TDMatrix(rewardIsPresentInEachTrial, :))/numRewardedTrials;
averageTD_learningErrorVector_ForRewardedTrials=sum(TD_learningErrorMatrix(rewardIsPresentInEachTrial, :))/numRewardedTrials;

% for the unrewarded trials
averageValueFnVector_ForUnrewardedTrials=sum(valueFnMatrix(~rewardIsPresentInEachTrial, :))/numUnrewardedTrials;
averageTDVector_ForUnrewardedTrials=sum(TDMatrix(~rewardIsPresentInEachTrial, :))/numUnrewardedTrials;
averageTD_learningErrorVector_ForUnrewardedTrials=sum(TD_learningErrorMatrix(~rewardIsPresentInEachTrial, :))/numUnrewardedTrials;

% for all trials
averageValueFnVector_ForAllTrials=sum(valueFnMatrix)/totalNumTrials;
averageTDVector_ForAllTrials=sum(TDMatrix)/totalNumTrials;
averageTD_learningErrorVector_ForAllTrials=sum(TD_learningErrorMatrix)/totalNumTrials;

%%
subplot(3,1,1);
plot(timeStepsWithinTrial, averageValueFnVector_ForRewardedTrials, 'r'); hold on;
plot(timeStepsWithinTrial, averageValueFnVector_ForUnrewardedTrials, 'b');
plot(timeStepsWithinTrial, averageValueFnVector_ForAllTrials,'k--', 'LineWidth',2);
title('The value function','Interpreter','latex','FontSize',14);
ylabel('$\hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);

subplot(3,1,2);
plot(timeStepsWithinTrial, averageTDVector_ForRewardedTrials, 'r'); hold on;
plot(timeStepsWithinTrial, averageTDVector_ForUnrewardedTrials, 'b');
plot(timeStepsWithinTrial, averageTDVector_ForAllTrials, 'k--', 'LineWidth',2);
title('The TD of the value','Interpreter','latex','FontSize',14);
ylabel('$\Delta \hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);

subplot(3,1,3);
plot(timeStepsWithinTrial, averageTD_learningErrorVector_ForRewardedTrials, 'r'); hold on;
plot(timeStepsWithinTrial, averageTD_learningErrorVector_ForUnrewardedTrials, 'b');
plot(timeStepsWithinTrial, averageTD_learningErrorVector_ForAllTrials,'k--', 'LineWidth',2);
title('The TD learning error','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
ylabel('$\delta(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
legend({'Rewarded trials average', 'Unrewarded trials average', 'All trials average'}, 'Interpreter','latex','FontSize',10);

%% 5(c)
figure;
TD_learningErrorSamples=-3:0.01:3;
dopamineActivitySamples=DA(TD_learningErrorSamples);
plot(TD_learningErrorSamples, dopamineActivitySamples,'k','LineWidth',2);
hold on; plot([0.27,0.27], [0,0.27], 'k:','LineWidth',2);
axis equal;

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title('The dopamine activity as a function of the TD learning error','Interpreter','latex','FontSize',14);
xlabel('TD learning error $x$ $(\delta (t))$','Interpreter','latex','FontSize',14);
ylabel('The dopamine activity DA(x)','Interpreter','latex','FontSize',14);

%%
figure;
% subplot(2,1,1);
% plot(timeStepsWithinTrial, averageTD_learningErrorVector_ForRewardedTrials, 'r'); hold on;
% plot(timeStepsWithinTrial, averageTD_learningErrorVector_ForUnrewardedTrials, 'b');
% plot(timeStepsWithinTrial, averageTD_learningErrorVector_ForAllTrials,'k--', 'LineWidth',2);
% title('The TD learning error','Interpreter','latex','FontSize',14);
% xlabel('Time t','Interpreter','latex','FontSize',14);
% ylabel('$\delta(t)$','Interpreter','latex','FontSize',14);
% xlim([0,25]);
% legend({'Rewarded trials average', 'Unrewarded trials average', 'All trials average'}, 'Interpreter','latex','FontSize',10);

%subplot(2,1,2);
dopamineActivityMatrix=DA(TD_learningErrorMatrix); % this has already gone down to 100 entries rather than 1000

averageDopamineActivityVector_ForRewardedTrials=sum(dopamineActivityMatrix(rewardIsPresentInEachTrial, :))/numRewardedTrials;
averageDopamineActivityVector_ForUnrewardedTrials=sum(dopamineActivityMatrix(~rewardIsPresentInEachTrial, :))/numUnrewardedTrials;
averageDopamineActivityVector_ForAllTrials=sum(dopamineActivityMatrix)/totalNumTrials;

plot(timeStepsWithinTrial, averageTD_learningErrorVector_ForRewardedTrials, 'r:'); hold on;
plot(timeStepsWithinTrial, averageTD_learningErrorVector_ForUnrewardedTrials, 'b:');
plot(timeStepsWithinTrial, averageTD_learningErrorVector_ForAllTrials,'k:', 'LineWidth',2);

plot(timeStepsWithinTrial, averageDopamineActivityVector_ForRewardedTrials, 'r'); hold on;
plot(timeStepsWithinTrial, averageDopamineActivityVector_ForUnrewardedTrials, 'b');
plot(timeStepsWithinTrial, averageDopamineActivityVector_ForAllTrials, 'k--', 'LineWidth',2);

title('The learning error and the dopamine activity','Interpreter','latex','FontSize',14);
ylabel('The learning error and the dopamine activity','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
xlim([0,25]);
legend({'$\delta(t)$: Rewarded trials average','$\delta(t)$: Unrewarded trials average','$\delta(t)$: All trials average',...
    'DA(t): Rewarded trials average', 'DA(t): Unrewarded trials average', 'DA(t): All trials average'},...
    'Interpreter','latex','FontSize',10);

%%
figure;
numSelectedPlots=length(selectedPlotIndices);
for ind=1:length(selectedPlotIndices)
    plot(timeStepsWithinTrial, TD_learningErrorMatrix(ind, :), 'k'); hold on;
end
title('The learning error','Interpreter','latex','FontSize',14);
ylabel('$\delta (t)$','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
xlim([0,25]);
ylim([-0.5,2]);

%%
figure;
numSelectedPlots=length(selectedPlotIndices);
for ind=1:length(selectedPlotIndices)
    plot(timeStepsWithinTrial, dopamineActivityMatrix(ind, :), 'k'); hold on;
end
title('The dopamine activity','Interpreter','latex','FontSize',14);
ylabel('$DA(t)$','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
xlim([0,25]);
ylim([-0.5,2]);

% for ind=1:length(selectedPlotIndices)
%     plotInd=ind;
%     plot(timeStepsWithinTrial, dopamineActivityMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
% end
% title('The dopamine activity','Interpreter','latex','FontSize',14);
% ylabel('$DA(x_t)$','Interpreter','latex','FontSize',14);
% xlim([0,25]);
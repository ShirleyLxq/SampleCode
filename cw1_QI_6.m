%% I.6
figure;
gamma=1; epsilon=0.01; % 0.01;
N=1000; % ntrials
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
%rewardIsPresentInEachTrial=zeros(N, 1);
%pArray=[0:0.25:1];
pArray=[0:0.05:1];
num_p=length(pArray);
peakAtTheTimeOfStimulusArray=zeros(1,num_p);
peakAtTheTimeOfRewardArray=zeros(1,num_p);
timeAroundTheTimeOfStimulus=10+[-2:delta_t:2];
indAroundTheTimeOfStimulus=convertTimeToArrayIndex(timeAroundTheTimeOfStimulus, delta_t);
timeAroundTheTimeOfReward=20+[-2:delta_t:2];
indAroundTheTimeOfReward=convertTimeToArrayIndex(timeAroundTheTimeOfReward, delta_t);

for pInd=1:num_p
p=pArray(pInd); % chance of p that reward is present; 1-p it is not
averageDopamineActivityVector_ForAllTrials=computeAverageDAAcrossTrials(N,p,gamma,epsilon,...
    timeStepsWithinTrial,numTimeStepsWithinTrial,delta_t,S_compositeMatrix,S_compositeMatrix_oneStepBackInTime,weights,TD_learningErrorMatrix);

DA_AroundTheTimeOfStimulus=averageDopamineActivityVector_ForAllTrials(indAroundTheTimeOfStimulus);
[maxValStimulus,maxIndStimulus]=max(DA_AroundTheTimeOfStimulus);
peakAtTheTimeOfStimulusArray(pInd)=maxValStimulus;

DA_AroundTheTimeOfReward=averageDopamineActivityVector_ForAllTrials(indAroundTheTimeOfReward);
[maxValReward,maxIndReward]=max(DA_AroundTheTimeOfReward);
peakAtTheTimeOfRewardArray(pInd)=maxValReward;

if(p==0)
    plot(timeStepsWithinTrial, averageDopamineActivityVector_ForAllTrials,'k', 'LineWidth',2); hold on;
else
    plot(timeStepsWithinTrial, averageDopamineActivityVector_ForAllTrials, 'LineWidth', 1); hold on;
end
end
%%
legend({'p=0', 'p=0.25', 'p=0.5', 'p=0.75', 'p=1'}, 'Interpreter','latex','FontSize',14);
title('The dopamine activity with different reward probabilities $p$','Interpreter','latex','FontSize',14);
ylabel('$DA(t)$','Interpreter','latex','FontSize',14);
xlabel('Time t','Interpreter','latex','FontSize',14);
xlim([0,25]);

%% I.7
figure;
subplot(2,1,1);
plot(pArray, peakAtTheTimeOfStimulusArray,'*-');
title('The maximum dopamine activity around the time of the stimulus','Interpreter','latex','FontSize',14);
ylabel('DA(p)','Interpreter','latex','FontSize',14);

subplot(2,1,2);
plot(pArray, peakAtTheTimeOfRewardArray, '*-'); hold on;
ylabel('DA(p)','Interpreter','latex','FontSize',14);

yyaxis right
pSamples=0:0.01:1;
plot(pSamples, bernoulliEntropy(pSamples), '--');
title('The maximum dopamine activity around the time of the reward','Interpreter','latex','FontSize',14);
ylabel('H(p)','Interpreter','latex','FontSize',14);
legend({'DA(p)', 'The entropy H(p)'}, 'Interpreter','latex','FontSize',14);
xlabel('The probability of reward $p$', 'Interpreter','latex','FontSize',14);
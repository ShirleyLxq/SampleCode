function averageDopamineActivityVector_ForAllTrials=computeAverageDAAcrossTrials(N,p,gamma,epsilon,...
    timeStepsWithinTrial,numTimeStepsWithinTrial,delta_t,S_compositeMatrix,S_compositeMatrix_oneStepBackInTime,weights,TD_learningErrorMatrix)
for n=1:N % N trials
    % determine whether reward is presence or absence
    rewardIsPresent=(rand(1)<p); % uniform distribution between 0 and 1
    %rewardIsPresentInEachTrial(n)=rewardIsPresent;
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
%     valueFnVector=weights*S_compositeMatrix';
%     valueFnMatrix(n,:)=valueFnVector;
     TDVector=weights*S_compositeMatrix'-weights*S_compositeMatrix_oneStepBackInTime';
%     TDMatrix(n,:)=TDVector;
    TD_learningErrorVector=reward_oneStepBackInTime+TDVector;
    TD_learningErrorMatrix(n,:)=TD_learningErrorVector;
end

%% plot three graphs
totalNumTrials=100;
selectedTrials= [N-totalNumTrials+1: N]; % retain only the data for the last 100 trials
%selectedTrials= [1:totalNumTrials]; % retain only the data for the last 100 trials

%rewardIsPresentInEachTrial=logical(rewardIsPresentInEachTrial(selectedTrials));
%valueFnMatrix=valueFnMatrix(selectedTrials, :);
%TDMatrix=TDMatrix(selectedTrials, :);
TD_learningErrorMatrix=TD_learningErrorMatrix(selectedTrials, :);
%numRewardedTrials=sum(rewardIsPresentInEachTrial);
%numUnrewardedTrials=totalNumTrials-numRewardedTrials;

%%
dopamineActivityMatrix=DA(TD_learningErrorMatrix); % this has already gone down to 100 entries rather than 1000
%averageDopamineActivityVector_ForRewardedTrials=sum(dopamineActivityMatrix(rewardIsPresentInEachTrial, :))/numRewardedTrials;
%averageDopamineActivityVector_ForUnrewardedTrials=sum(dopamineActivityMatrix(~rewardIsPresentInEachTrial, :))/numUnrewardedTrials;
averageDopamineActivityVector_ForAllTrials=sum(dopamineActivityMatrix)/totalNumTrials;
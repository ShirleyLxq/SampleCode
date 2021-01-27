function plotEvolutionTripleGraphs(timeStepsWithinTrial, selectedPlotIndices, valueFnMatrix, TDMatrix, TD_learningErrorMatrix)
figure;
subplot(3,1,1);
%selectedPlotIndices=1:10:31;
numSelectedPlots=length(selectedPlotIndices);
colourList=colormap(jet(numSelectedPlots));
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, valueFnMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The value function','Interpreter','latex','FontSize',14);
ylabel('$\hat{V}(t)$','Interpreter','latex','FontSize',14);
xlim([0,25]);
%ylim([0,3]);
    
subplot(3,1,2);
for ind=1:length(selectedPlotIndices)
    plotInd=selectedPlotIndices(ind);
    plot(timeStepsWithinTrial, TDMatrix(plotInd, :), 'Color' ,colourList(ind, :)); hold on;
end
title('The TD of the value','Interpreter','latex','FontSize',14);
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

cb=colorbar; %('Location', 'SouthOutside', 'Limits', [1 numSelectedPlots/(numSelectedPlots+10)]);
caxis([selectedPlotIndices(1), selectedPlotIndices(end)]);
set(cb, 'YTick',[selectedPlotIndices(1), selectedPlotIndices(end)])

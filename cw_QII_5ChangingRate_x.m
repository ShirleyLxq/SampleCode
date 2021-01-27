function cw_QII_5ChangingRate_x(isRunningRemotely)
rate_xArray=5:5:20;
numRate_x=length(rate_xArray);
rate_EArray=zeros(1, numRate_x);
rate_IArray=zeros(1, numRate_x);?

for rate_xInd=1:numRate_x
    rate_x=rate_xArray(rate_xInd);
    fprintf('rate_xInd=%d, rate_x=%d...\n',rate_xInd, rate_x);
    [rate_E,rate_I]=cw_QII_5(rate_x);
    rate_EArray(rate_xInd)=rate_E;
    rate_IArray(rate_xInd)=rate_I;
end
varyingRatesFigHandle=figure;
plot(rate_xArray, rate_xArray); hold on;
plot(rate_xArray, rate_EArray);
plot(rate_xArray, rate_IArray);
title('The firing rates of the three populations of neurons','Interpreter','latex','FontSize',14);
ylabel('The firing rates $r_X, r_E, r_I$','Interpreter','latex','FontSize',14);
xlabel('The firing rate of the X population $r_X$', 'Interpreter','latex','FontSize',14);
legend({'The external population $r_X$', 'The excitatory population $r_E$', 'The inhibitory population $r_I$'}, 'Interpreter','latex','FontSize',10);
varyingRatesFigName=sprintf('varyingRates-%s', datestr(now));
saveas(varyingRatesFigHandle, varyingRatesFigName, 'jpg');
saveas(varyingRatesFigHandle, varyingRatesFigName, 'fig');

if(isRunningRemotely)
    exit;
end
function [indOfTheInputFromXToEachE, indOfTheInputFromEToEachE, indOfTheInputFromIToEachE,...
    indOfTheInputFromXToEachI, indOfTheInputFromEToEachI, indOfTheInputFromIToEachI]= fullNetwork_generateConnectivityMatrices(N,K)
% function [indOfTheInputFromXToEachX, indOfTheInputFromEToEachX, indOfTheInputFromIToEachX,...
%     indOfTheInputFromXToEachE, indOfTheInputFromEToEachE, indOfTheInputFromIToEachE,...
%     indOfTheInputFromXToEachI, indOfTheInputFromEToEachI, indOfTheInputFromIToEachI]= fullNetwork_generateConnectivityMatrices(N,K)

%indOfTheInputFromXToEachX =fullNetwork_generateConnectivityMatrix(N,K, 'same');
%indOfTheInputFromEToEachX =fullNetwork_generateConnectivityMatrix(N,K, 'different');
%indOfTheInputFromIToEachX =fullNetwork_generateConnectivityMatrix(N,K, 'different');

indOfTheInputFromXToEachE =fullNetwork_generateConnectivityMatrix(N,K, 'different');
indOfTheInputFromEToEachE =fullNetwork_generateConnectivityMatrix(N,K, 'same');
indOfTheInputFromIToEachE =fullNetwork_generateConnectivityMatrix(N,K, 'different');

indOfTheInputFromXToEachI =fullNetwork_generateConnectivityMatrix(N,K, 'different');
indOfTheInputFromEToEachI =fullNetwork_generateConnectivityMatrix(N,K, 'different');
indOfTheInputFromIToEachI =fullNetwork_generateConnectivityMatrix(N,K, 'same');
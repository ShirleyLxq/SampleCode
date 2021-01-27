function indOfTheInputFromXToEachE =fullNetwork_generateConnectivityMatrix(N,K, connectionTypeFlag)
indOfTheInputFromXToEachE=zeros(N, K); % each of the N neurons in the pop E is connected to K neurons in pop X
% each row is for one E neuron, containing the indices of its K presyn X partners
for E_neuronInd=1:N
    switch connectionTypeFlag
        case 'different'
            indOfTheInputFromX=randperm(N, K); % a row vector
            indOfTheInputFromXToEachE(E_neuronInd, :)=indOfTheInputFromX; % K presyn X partners
        case 'same'
            % avoid self connections
            usableIndices=1:N;
            usableIndices(usableIndices==E_neuronInd)=[];
            indOfTheInputFromX= datasample(usableIndices,K,'Replace',false);
            indOfTheInputFromXToEachE(E_neuronInd, :)=indOfTheInputFromX;
        otherwise
            error('Invalid connectionTypeFlag! \n');
    end
end
function sequence = generateMarkovSequence(transitionTable, initialState, sequenceLength)
    % Converta a tabela em uma matriz
    transitionMatrix = table2array(transitionTable);
    
    % Verifique o n�mero de estados
    numStates = size(transitionMatrix, 1);
    
    % Inicialize a sequ�ncia de estados
    sequence = zeros(1, sequenceLength);
    
    % Defina o estado inicial
    currentState = initialState;
    
    % Gere a sequ�ncia de estados
    for i = 1:sequenceLength
        % Calcule as probabilidades de transi��o a partir do estado atual
        transitionProbabilities = transitionMatrix(currentState, :);
        
        % Escolha o pr�ximo estado com base nas probabilidades de transi��o
        nextState = find(rand() <= cumsum(transitionProbabilities), 1);
        
        % Atualize o estado atual
        currentState = nextState;
        
        % Armazene o estado na sequ�ncia
        sequence(i) = currentState;
    end
end

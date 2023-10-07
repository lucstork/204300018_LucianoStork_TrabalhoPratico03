function sequence = generateMarkovSequence(transitionTable, initialState, sequenceLength)
    % Converta a tabela em uma matriz
    transitionMatrix = table2array(transitionTable);
    
    % Verifique o número de estados
    numStates = size(transitionMatrix, 1);
    
    % Inicialize a sequência de estados
    sequence = zeros(1, sequenceLength);
    
    % Defina o estado inicial
    currentState = initialState;
    
    % Gere a sequência de estados
    for i = 1:sequenceLength
        % Calcule as probabilidades de transição a partir do estado atual
        transitionProbabilities = transitionMatrix(currentState, :);
        
        % Escolha o próximo estado com base nas probabilidades de transição
        nextState = find(rand() <= cumsum(transitionProbabilities), 1);
        
        % Atualize o estado atual
        currentState = nextState;
        
        % Armazene o estado na sequência
        sequence(i) = currentState;
    end
end

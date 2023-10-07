clear; close all; clc;

% ------------------------------------------------- Parte 01 -----------------------------------------------------

% Leia o arquivo hexokinase.txt como um sequenciamento FASTA
fasta_data = fastaread('hexokinase.txt');

% Acesse o sequenciamento de DNA a partir da estrutura FASTA
sequence = fasta_data.Sequence;

% Defina rótulos para os estados (A, C, G, T)
states = {'A', 'C', 'G', 'T'};

% Inicialize as variáveis para contar as transições
num_states = length(states);
transition_counts = zeros(num_states, num_states); % Matriz NxN

% Analise a sequência de DNA e conte as transições
for i = 1:(length(sequence) - 1)
    current_nucleotide = sequence(i);
    next_nucleotide = sequence(i + 1);
    
    % Mapeie os nucleotídeos para índices
    current_index = find(strcmp(states, current_nucleotide));
    next_index = find(strcmp(states, next_nucleotide));
    
    % Atualize a contagem de transição
    transition_counts(current_index, next_index) = transition_counts(current_index, next_index) + 1;
end

% Calcule as probabilidades de transição
transition_probabilities = transition_counts ./ sum(transition_counts, 2);

% Crie uma tabela para representar a matriz de transição
transition_table = array2table(transition_probabilities, 'RowNames', states, 'VariableNames', states);

% Exiba a matriz de transição como uma tabela
disp('Matriz de Transição da Cadeia de Markov:');

disp(transition_table);

% Gere uma sequência de estados com a função generateMarkovSequence
initialState = 1; % Escolha um estado inicial (1, 2, 3 ou 4)
sequenceLength = 1000; % Comprimento da sequência desejada
sequence = generateMarkovSequence(transition_table, initialState, sequenceLength);

% Mapeie os números de volta para os rótulos de estados
state_labels = {'A', 'C', 'G', 'T'};
sequence_text = cellfun(@(x) state_labels{x}, num2cell(sequence), 'UniformOutput', false);
sequence_text = [sequence_text{:}];

% Exiba a sequência em formato de texto
disp('Sequência de Estados Gerada (Formato de Texto):');
disp(sequence_text);

% Calcule a distribuição estacionária
[V, D] = eig(transition_probabilities');
[~, index] = max(abs(diag(D)) - 1);
stationary_distribution = abs(V(:, index)) / sum(abs(V(:, index)));

% Crie uma tabela para representar a distribuição estacionária
DistribuicaoEstacionaria = array2table(stationary_distribution, 'RowNames', states, 'VariableNames', {'Distribuicao_Estacionaria'});

% Exiba a distribuição estacionária como uma tabela
disp("  ")
disp(DistribuicaoEstacionaria);

% Calcule a distribuição estacionária a partir da matriz de transição
calculated_stationary_distribution = stationary_distribution' * transition_probabilities;

% Verifique a consistência entre as distribuições
is_consistent = ismembertol(stationary_distribution, calculated_stationary_distribution, 1e-6);

% Exiba o resultado da verificação
if all(is_consistent)
    disp('A distribuição estacionária é consistente com a matriz de transição.');
    disp('A verificação é consistente porque a distribuição estacionária é o vetor próprio correspondente ao autovalor 1 da matriz de transição.');
    disp('Segue, a nível de comparação, os resultados que comprovam tal afimação:')
else
    disp('A distribuição estacionária NÃO é consistente com a matriz de transição.');
end

% Exiba as distribuições estacionárias separadamente
disp('Distribuição Estacionária:');
disp(stationary_distribution');
disp('Distribuição Calculada a partir da Matriz de Transição:');
disp(calculated_stationary_distribution);


% ------------------------------------------------- Parte 02 -----------------------------------------------------
disp("  ")
disp("---------------------------------Parte 02-------------------------------")
disp("  ")

%  Cálculo da taxa de entropia da matriz de transição.
entropy_transition = calculateEntropy(transition_probabilities);

% Exiba a taxa de entropia da matriz de transição.
disp('Taxa de Entropia da Matriz de Transição:');
disp(entropy_transition);

    
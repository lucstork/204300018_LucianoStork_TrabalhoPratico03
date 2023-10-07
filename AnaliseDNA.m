clear; close all; clc;

% ------------------------------------------------- Parte 01 -----------------------------------------------------

% Leia o arquivo hexokinase.txt como um sequenciamento FASTA
fasta_data = fastaread('hexokinase.txt');

% Acesse o sequenciamento de DNA a partir da estrutura FASTA
sequence = fasta_data.Sequence;

% Defina r�tulos para os estados (A, C, G, T)
states = {'A', 'C', 'G', 'T'};

% Inicialize as vari�veis para contar as transi��es
num_states = length(states);
transition_counts = zeros(num_states, num_states); % Matriz NxN

% Analise a sequ�ncia de DNA e conte as transi��es
for i = 1:(length(sequence) - 1)
    current_nucleotide = sequence(i);
    next_nucleotide = sequence(i + 1);
    
    % Mapeie os nucleot�deos para �ndices
    current_index = find(strcmp(states, current_nucleotide));
    next_index = find(strcmp(states, next_nucleotide));
    
    % Atualize a contagem de transi��o
    transition_counts(current_index, next_index) = transition_counts(current_index, next_index) + 1;
end

% Calcule as probabilidades de transi��o
transition_probabilities = transition_counts ./ sum(transition_counts, 2);

% Crie uma tabela para representar a matriz de transi��o
transition_table = array2table(transition_probabilities, 'RowNames', states, 'VariableNames', states);

% Exiba a matriz de transi��o como uma tabela
disp('Matriz de Transi��o da Cadeia de Markov:');

disp(transition_table);

% Gere uma sequ�ncia de estados com a fun��o generateMarkovSequence
initialState = 1; % Escolha um estado inicial (1, 2, 3 ou 4)
sequenceLength = 1000; % Comprimento da sequ�ncia desejada
sequence = generateMarkovSequence(transition_table, initialState, sequenceLength);

% Mapeie os n�meros de volta para os r�tulos de estados
state_labels = {'A', 'C', 'G', 'T'};
sequence_text = cellfun(@(x) state_labels{x}, num2cell(sequence), 'UniformOutput', false);
sequence_text = [sequence_text{:}];

% Exiba a sequ�ncia em formato de texto
disp('Sequ�ncia de Estados Gerada (Formato de Texto):');
disp(sequence_text);

% Calcule a distribui��o estacion�ria
[V, D] = eig(transition_probabilities');
[~, index] = max(abs(diag(D)) - 1);
stationary_distribution = abs(V(:, index)) / sum(abs(V(:, index)));

% Crie uma tabela para representar a distribui��o estacion�ria
DistribuicaoEstacionaria = array2table(stationary_distribution, 'RowNames', states, 'VariableNames', {'Distribuicao_Estacionaria'});

% Exiba a distribui��o estacion�ria como uma tabela
disp("  ")
disp(DistribuicaoEstacionaria);

% Calcule a distribui��o estacion�ria a partir da matriz de transi��o
calculated_stationary_distribution = stationary_distribution' * transition_probabilities;

% Verifique a consist�ncia entre as distribui��es
is_consistent = ismembertol(stationary_distribution, calculated_stationary_distribution, 1e-6);

% Exiba o resultado da verifica��o
if all(is_consistent)
    disp('A distribui��o estacion�ria � consistente com a matriz de transi��o.');
    disp('A verifica��o � consistente porque a distribui��o estacion�ria � o vetor pr�prio correspondente ao autovalor 1 da matriz de transi��o.');
    disp('Segue, a n�vel de compara��o, os resultados que comprovam tal afima��o:')
else
    disp('A distribui��o estacion�ria N�O � consistente com a matriz de transi��o.');
end

% Exiba as distribui��es estacion�rias separadamente
disp('Distribui��o Estacion�ria:');
disp(stationary_distribution');
disp('Distribui��o Calculada a partir da Matriz de Transi��o:');
disp(calculated_stationary_distribution);


% ------------------------------------------------- Parte 02 -----------------------------------------------------
disp("  ")
disp("---------------------------------Parte 02-------------------------------")
disp("  ")

%  C�lculo da taxa de entropia da matriz de transi��o.
entropy_transition = calculateEntropy(transition_probabilities);

% Exiba a taxa de entropia da matriz de transi��o.
disp('Taxa de Entropia da Matriz de Transi��o:');
disp(entropy_transition);

    
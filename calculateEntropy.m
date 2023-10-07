function entropy = calculateEntropy(input_data)
    if ischar(input_data)  % Verifica se a entrada é uma sequência de caracteres
        % Se a entrada for uma sequência de caracteres, calcule a entropia de texto
        input_data = upper(input_data);
        input_data = regexprep(input_data, '[^A-Z]', ''); % Remove não-alfabéticos
        freqs = histc(double(input_data) - 64, 1:26);
        freqs = freqs(freqs > 0);
        p = freqs / length(input_data);
        entropy = -sum(p .* log2(p));
    elseif isnumeric(input_data) && ismatrix(input_data) % Verifica se a entrada é uma matriz numérica
        % Se a entrada for uma matriz numérica, calcule a entropia da matriz

        entropy = -sum(input_data(:) .* log2(input_data(:)));
    else
        error('Tipo de entrada não suportado. A entrada deve ser uma sequência de caracteres ou uma matriz numérica.');
    end
end

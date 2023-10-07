function entropy = calculateEntropy(input_data)
    if ischar(input_data)  % Verifica se a entrada � uma sequ�ncia de caracteres
        % Se a entrada for uma sequ�ncia de caracteres, calcule a entropia de texto
        input_data = upper(input_data);
        input_data = regexprep(input_data, '[^A-Z]', ''); % Remove n�o-alfab�ticos
        freqs = histc(double(input_data) - 64, 1:26);
        freqs = freqs(freqs > 0);
        p = freqs / length(input_data);
        entropy = -sum(p .* log2(p));
    elseif isnumeric(input_data) && ismatrix(input_data) % Verifica se a entrada � uma matriz num�rica
        % Se a entrada for uma matriz num�rica, calcule a entropia da matriz

        entropy = -sum(input_data(:) .* log2(input_data(:)));
    else
        error('Tipo de entrada n�o suportado. A entrada deve ser uma sequ�ncia de caracteres ou uma matriz num�rica.');
    end
end

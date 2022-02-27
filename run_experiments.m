%% run_experiments

entrada = [ 
            5,5,179.313
            7,7,351.454
           11,2,448.283
           11,4,491.318
           11,6,563.044
           11,11,867.876
           13,11,1040.017
           13,12,1122.501
                            ];

% entrada = (0.1:0.1:1).';

nexp = size(entrada,1); % numero de experimentos
nrepet = 1; % numero de repeti��es de cada experimentos
for i = 1:nexp
    for j = 1:nrepet
        calc_IU_placa([2,2],[entrada(i,1),entrada(i,2)],[63,63],entrada(i,3))
%         calc_IU_cilindro_com_tampas(entrada(i))
    end
end


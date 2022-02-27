function calc_IU_cilindro(L,modos,Ne,freq_adim)

    %%%%%%%%%%%%% Designando variáveis de entrada %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    L = L(:);
    modos = modos(:);
    Ne = Ne(:);
    
    raio = L(1);      % raio do cilindro
    altura = L(2);      % altura do cilindro
    m = modos(1);  % modo de vibração na direção radial
    n = modos(2);   % modo de vibração na direção z
    Ne_arco = Ne(1);  % número de elementos por arco 
    Ne_z = Ne(2);  % número de elementos na direção z
    
    %%%%%%%%%%%%% Constantes  Físicas  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     altura = 5; % altura desejada para o cilindro
%     raio = 0.5; % raio desejado para o cilindro
    
    c = 343;  % velocidade do som
    rho = 1.29; % massa específicado meio (ar)
        
    %%%%%%% Calculo de variáveis importantes na placa %%%%%%%%
    
    tic % inicia o cronômetro
    
    teta = linspace(0,pi,Ne_arco + 1); %  angulos discretos da malha
    z = linspace(0,altura,Ne_z + 1); % valores discretos em z
    
    
    A = 2*pi*raio*altura;  % área superficial do cilindro
    Ne = Ne_arco*Ne_z; % número total de elementos na malha
    Np = (Ne_arco + 1)*(Ne_z + 1); % número de pontos nodais
   
    kz= n*pi/altura; % número de onda na direção axial 
    kr = m*pi/(2*pi*raio); % número de onda na direção radial
    kf = sqrt(kz^2+kr^2); % nï¿½mero de onda livre
    freq_coinc = kf*c/(2*pi);
    frequencia = freq_adim*freq_coinc; % frequencia de excitacao
    omega = 2*pi*frequencia; %  frequencia angular
    k = omega/c; % número de onda

   
%     kr = ones(Ne_arco+1,1)*sqrt(k^2 - kz.^2); % número de onda na direção radial
    delta = 10^(-7); % tolerancia do truncamento

    if  kz^2 + kr^2 <= k^2
        disp('Modo de Superficie')
    elseif (kr > k) && (kz > k)
        disp('Modo de Canto')
    elseif (kr > k && kz < k)||(kr < k && kz > k)
        disp('Modo de Borda')
    else
        disp('Modo de Transicao')
    end
      
    
    %%%% Calculo do Deslocamento e Velocidade %%%%%%%%%
    
    d = zeros(Ne_arco+1,Ne_z+1);
    for i = 1:Ne_arco + 1
        for j =1:Ne_z + 1
            d(i,j) = 2/sqrt(A)*sin(kz*z(j))*sin(kr*teta(i));   % Calculo do deslocamento na superfície da placa pela formula analítica
        end
    end 
    
    V = sqrt(-1)*omega*d; % calculo da velocidade normal

    v = V(:);   % vetor com os pontos nodais da velocidade 
    r = size(v,1);
    clear d % limpando variáveis
    
    % Calculo da Solucao Analitica via Abordagem de Fourier
    
%     kr1 = ifftshift( (2*pi/Lx)*(-Ne_x/2:Ne_x/2)  ); % nï¿½meros de onda da direï¿½ï¿½o x 
%     kz1 = ifftshift( (2*pi/Ly)*(-Ne_y/2:Ne_y/2)  ); % nï¿½meros de onda da direï¿½ï¿½o y
    kr1 = (2*pi/(2*pi*raio))*[0:((Ne_arco+1)/2-1) (-(Ne_arco+1)/2):-1];
    kz1 = (2*pi/altura)*[0:((Ne_z+1)/2-1) (-(Ne_z+1)/2):-1];
    
    kf1 = zeros(Ne_arco + 1,Ne_z + 1);    
    for i = 1:Ne_arco + 1
        for j = 1:Ne_z + 1
            kf1(i,j) = sqrt(k^2 - kr1(i)^2 - kz1(j)^2);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tV = fft2(V); % velocidade angular
    tP = rho*c*k*tV./kf1;  % pressao angular
    P = ifft2(tP);   % pressao 
    I2 = 1/2*real(P.*conj(V)); % intensidade acï¿½stica
    Pot2 = integral_dupla(2*pi*raio/Ne_arco, altura/Ne_z, I2); % Potï¿½ncia sonora
    NPS2 = 10*log10((Pot2)/(10^(-12)));
    fprintf('Intensidade Acústica Analítica - NPS: %f \n',NPS2)
    
    %%%%%%%%%%%%%%%%% Calculo das Grandezas Supersï¿½nicas %%%%%%%%%%%%%%%%
     
    clear D V P kf1
    tPS = zeros(Ne_arco + 1,Ne_z + 1);  % pressï¿½o angular supersï¿½nica
    tVS = zeros(Ne_arco + 1,Ne_z + 1);  % velocidade angular supersï¿½nica
    for i = 1:Ne_arco + 1
        for j =1:Ne_z + 1
            if  kr1(i)^2 + kz1(j)^2 <= k^2
                tPS(i,j) = tP(i,j);
                tVS(i,j) = tV(i,j);
            end
        end
    end
    VS = ifft2(tVS);   % velocidade supersonica
    PS = ifft2(tPS);   % pressao 
    IS = 1/2*real(PS.*conj(VS)); % intensidade acï¿½stica
    Pot2 = integral_dupla(2*pi*raio/Ne_arco, altura/Ne_z,IS); % Potï¿½ncia sonora
    NPS_IS = 10*log10((Pot2)/(10^(-12)));
    fprintf('Intensidade Supersônica Analitica - NPS: %f \n',NPS_IS)
    
    clear VS tVS  % limpando variï¿½veis


    %%%%%%%%%%%%% Montando Matrizes de Influência %%%%%%%%%%%%%%
    
    [coord,edof,conec,dof] = topologia_cilindro(raio,altura,Ne_arco,Ne_z);
    [Ex,Ey,Ez] = coordxtr(edof,coord,dof,4);
 
    H = zeros(r,r);
    G = H;

    for k=1:r
        for j=1:Ne
            [He,Ge] = bem_infl4q(coord(k,:),Ex(j,:),Ey(j,:),Ez(j,:),[omega c rho],-1);
            H = bem_assem(edof,H,He,k,j);
            G = bem_assem(edof,G,Ge,k,j);
        end
    end

    H = H + 1/2*eye(r);
    
    coord2 = coord; coord2(:,3) = 10*ones(Np,1);
    H_bar = zeros(r,r);
    G_bar = zeros(r,r);
    for k=1:r
        for j=1:Ne
            [He,Ge] = bem_infl4q(coord2(k,:),Ex(j,:),Ey(j,:),Ez(j,:),[omega c rho],-1);
            H_bar = bem_assem(edof,H_bar,He,k,j);
            G_bar = bem_assem(edof,G_bar,Ge,k,j);
        end
    end

    %%%%%%%%%% Calculo do Operador de Superfície e da Intensidade Acústica %%%%%%%%%%%%%%
    
    R = H\G;      % operador de superfície
    R_bar = G_bar - H_bar*R;
    p = R*v;      % calculo da pressao na superficie
    I = 1/2*real(p.*conj(v)); %  Intensidade acustica convencional 
    Pot = integral_dupla(2*pi*raio/Ne_arco, altura/Ne_z, I); % Potï¿½ncia sonora
    NPS = 10*log10((Pot)/(10^(-12)));
    disp(NPS)
    
    clear H G p % limpando variáveis
    
   
    %%%%%%%%%% Calculo do Operador de Potencia %%%%%%%%%%%
    
    %%%%%%%%%%%%%%% Pontos e Pesos de Gauss %%%%%%%%%%%%%%%%%%%%%%%
 
    xi=[-0.8611363116 -0.3399810436 0.3399810436 0.8611363116];
    eta=[-0.8611363116 -0.3399810436 0.3399810436 0.8611363116];
    wxi=[0.3478548451 0.6521451549 0.6521451549 0.3478548451];
    weta=[0.3478548451 0.6521451549 0.6521451549 0.3478548451];
    
    Q = zeros(r,r);
    for e = 1:Ne
        Qe = zeros(Np,4);
        for i=1:4
            for j=1:4
                [N,A,Area] = Matrizelemento(e,edof,R,Ex(e,:),Ey(e,:),Ez(e,:),xi(i),eta(j));
                Qe = Qe + wxi(i)*weta(j)*(A*(N'*N))*(Area/4);
            end
        end
        Q(:,edof(e,2))= Q(:,edof(e,2))+ Qe(:,1);
        Q(:,edof(e,3))= Q(:,edof(e,3))+ Qe(:,2);
        Q(:,edof(e,4))= Q(:,edof(e,4))+ Qe(:,3);
        Q(:,edof(e,5))= Q(:,edof(e,5))+ Qe(:,4);
    end
%     Q = 1/4*(transp(Q) + conj(Q));  % Operador de potência
    Q = 1/2*real(Q);
    Pot = v'*Q*v;  % Cálculo da potência sonora sem o uso da Eigen-Decomposition
    NPS = 10*log10((Pot)/(10^(-12)));
    fprintf('Intensidade Acústica Numérica - NPS: %f \n',NPS)
        
    %%%%%%%%%% Decomposição em autovalores e sua ordenação %%%%%%%%%%%
    
    [V,D] = eig(Q);
    lambda = diag(D);
    
    clear Q Qe D dof % limpando variáveis
    
    [~,order] = sort(abs(lambda),'descend');
    V = V(:,order);
    lambda = lambda(order);
    
    %%%%%%%%%% Aplicação do Critério de Truncamento e Calculo da Velocidade útil %%%%%%%%%%%

    vu = zeros(r,1); 
    soma = 0;
    pot = zeros(r,1);
    somas_pot = zeros(r,1);
    aux = 0;
    cond = 0;
    traco = sum(lambda);
    
    for i=1:r
            soma = soma + lambda(i);
            pot(i) = real(lambda(i)*abs(v'*V(:,i))^2);
            somas_pot(i) = aux + pot(i);
            aux = somas_pot(i);
            if cond == 0
                vu = vu + V(:,i)'*v*V(:,i);  % Velocidade útil
                if abs((1 - soma/traco)*(1 - somas_pot(i)/Pot)) < delta  
                    rc = i;  % posição truncada
                    cond = 1;
                end
            end
    end
    
    clear V % limpando variáveis

    %%%%%%%%%% Calculo da Pressão útil e da Intensidade útil %%%%%%%%%%%

    pu = R*vu; % pressão útil
    IU = 1/2*real(pu.*conj(vu)); % intensidade útil
    Pot_IU = somas_pot(rc); % Potência útil
    NPS_IU = 10*log10((Pot_IU)/(10^(-12))); % nível de potência sonora útil
    fprintf('Intensidade Útil - NPS: %f \n',NPS_IU) 
%     disp(rc)
    
    clear vu pu  % limpando variáveis
        
    %%%%%%%%%% Decomposiï¿½ï¿½o em valores singulares e sua ordenaï¿½ï¿½o %%%%%%%%%%%
    
    [~,S,M] = svd(R_bar);
    alpha = diag(S);
    clear S
            
    vu = zeros(r,1); 
    soma = 0;
    pot = zeros(r,1);
    somas_pot2 = zeros(r,1);
    aux = 0;
    cond = 0;
    delta = 10^(-3);
    traco = sum(alpha);
    for i=1:r
            soma = soma + alpha(i);
            pot(i) = alpha(i)^2*abs(v'*M(:,i))^2;
            somas_pot2(i) = aux + pot(i);
            aux = somas_pot2(i);
            if abs(1 - soma/traco) >= delta
                vu = vu + M(:,i)'*v*M(:,i);  % Velocidade util   
            end
            if cond == 0 && abs(1 - soma/traco) < delta  
                    rc2 = i;  % posicao truncada
                    cond = 1;
            end
    end
    Pot_IU2 = somas_pot2(rc2); % Potï¿½ncia ï¿½til
     
    clear M alpha  % limpando variï¿½veis

    %%%%%%%%%% Calculo da Pressï¿½o ï¿½til e da Intensidade ï¿½til %%%%%%%%%%%

    pu = R*vu; % pressï¿½o ï¿½til
    IU2 = 1/2*real(pu.*vu'.'); % intensidade ï¿½til
    NPS_IU2 = 10*log10((Pot_IU2)/(10^(-12))); % nï¿½vel de potï¿½ncia sonora ï¿½til
    fprintf('Intensidade Supersônica  segundo Magalhães - NPS: %f \n',NPS_IU2)

    
    % Calculo da Intensidade Nao-Negativa
    
    Teta = zeros(r,r);
    for e = 1:Ne
        Tetae = zeros(4,4);
        for i=1:size(xi,2)
            for j=1:size(eta,2)
                [N,~,Area] = Matrizelemento(e,edof,R,Ex(e,:),Ey(e,:),Ez(e,:),xi(i),eta(j));
                Tetae = Tetae + wxi(i)*weta(j)*(N'*N)*(Area/4);
            end
        end
        Teta(edof(e,2:5),edof(e,2:5))= Teta(edof(e,2:5),edof(e,2:5))+ Tetae;
    end
    Q = 0.5*real(R'*Teta);
    [Psi,Lambda]=eig(Q,Teta);
    beta= Psi*sqrt(Lambda)*Psi'*Teta*v;
    INN = abs(beta).^2;
    Pot_INN = beta'*Teta*conj(beta);
    NPS_INN = 10*log10((Pot_INN)/(10^(-12)));
    fprintf('Intensidade Não-Negativa - NPS: %f \n',real(NPS_INN))
    
    toc % para o cronômetro e exibe o tempo de operação
    
    %%%%%%%%%%%%%%%%%% Plotando Resultados %%%%%%%%%%%%%%%%%%%
    
    close all;
    figure
    
    subplot(3,3,1)
    plotagrandeza3(coord,conec,imag(v))% plota a intensidade acústica
    title('velocidade normal')
    
    subplot(3,3,2)
    plotagrandeza3(coord,conec,I2(:)) % plota a intensidade acustica analitica
    title('intensidade analitica')
%   
    subplot(3,3,3)
    plotagrandeza3(coord,conec,I) % plota a intensidade acústica
    title('intensidade acustica numerica')
%     saveas(gcf,['intensidade acustica convencional - ',num2str(raio),'-',num2str(altura),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_arco),'-',num2str(Ne_z),'- ',num2str(frequencia),'.png'])
    subplot(3,3,4)
    plotagrandeza3(coord,conec,IS(:)) % plota a intensidade acustica supersonica analitica
    title('intensidade supersonica analitica')
% 
    subplot(3,3,5)
    plotagrandeza3(coord,conec,IU) % plota a intensidade acústica útil
    title('intensidade acustica util')
%     saveas(gcf,['intensidade acustica util - ',num2str(raio),'-',num2str(altura),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_arco),'-',num2str(Ne_z),'- ',num2str(frequencia),'.png'])
    subplot(3,3,6)
    plotagrandeza3(coord,conec,INN) % plota a intensidade acústica útil
    title('intensidade acustica nao negativa')
    
    subplot (3,3,7)
    plotagrandeza3(coord,conec,IU2) % plota a intensidade acustica supersonica via Magalhaes
    title('intensidade supersonica Magalhaes')
    colorbar
%     saveas(gcf,['intensidade acustica nao negativa - ',num2str(raio),'-',num2str(altura),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_arco),'-',num2str(Ne_z),'- ',num2str(frequencia),'.png'])
    saveas(gcf,['cilindro grandezas - ',num2str(raio),'-',num2str(altura),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_arco),'-',num2str(Ne_z),'- ',num2str(frequencia),'.png'])
%     dados = [I, IU, I2(:), IS(:), INN];
%     save( ['cilindro grandezas - ',num2str(raio),'-',num2str(altura),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_arco),'-',num2str(Ne_z),'- ',num2str(frequencia)],'dados','-ascii','-double')
%  
%     sigma_IU = Pearson(IS(:),IU);
%     sigma_INN = Pearson(IS(:),INN);
% 
%     eta_IU = Spearman(IS(:),IU);
%     eta_INN = Spearman(IS(:),INN);
%     
%     if exist(['cilindro medidas - ',num2str(raio),'-',num2str(altura),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_arco),'-',num2str(Ne_z),'.txt'],'file')
%          dados2 = load(['cilindro medidas - ',num2str(raio),'-',num2str(altura),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_arco),'-',num2str(Ne_z),'.txt']);
%          dados2 = [dados2;[frequencia,NPS_IS,NPS_IU,NPS_INN,sigma_IU,sigma_INN,eta_IU,eta_INN]];
%     else
%          dados2 = [frequencia,NPS_IS,NPS_IU,NPS_INN,sigma_IU,sigma_INN, eta_IU,eta_INN];
%     end
%     save(['cilindro medidas - ',num2str(raio),'-',num2str(altura),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_arco),'-',num2str(Ne_z),'.txt'],'dados2','-ascii','-double')
    
%     % Plota a contribuição dos valores próprios 
%     
%     figure
%     plot(pot)
%     grid
%     xlabel('Valores próprios de velocidade')
%     ylabel('Potência Sonora (Watts)')
% 
%     % Plota a potência de acordo com os valores ordenados considerados
% 
    figure
    plot( (1:r)',somas_pot, rc*ones(1000, 1), linspace(min(somas_pot),somas_pot(rc), 1000 ),'b--')
    grid
    xlabel('Valores próprios de velocidade considerados')
    ylabel('Potência Sonora (Watts)')
    saveas(gcf,['cilindro convergencia potencia - ',num2str(raio),'-',num2str(altura),' - ',num2str(frequencia),'.png'])

end
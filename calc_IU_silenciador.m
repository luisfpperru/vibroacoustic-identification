function calc_IU_silenciador(frequencia)

    %%%%%%%%%%%%% Constantes  Físicas  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    c = 343;  % velocidade do som
    rho = 1.29; % massa específicado meio (ar)
        
    %%%%%%% Calculo de variáveis importantes na placa %%%%%%%%
    
    tic % inicia o cronômetro
    
    load('silenciador.mat','desl','coord','conec')
    Ne = size(conec,1); % número total de elementos na malha
    Np = size(coord,1); % número de pontos nodais
    d = desl(:);
    dof = (1:Ne)';
    edof= [dof conec];
        
%     frequencia = 525
%     freq_coinc = 240;
%     frequencia = freq_coinc*freq_adim;
    omega = 2*pi*frequencia; %  frequencia angular
    k = omega/c; % número de onda

%     kz = n*pi/altura; % número de onda na direção axial 
%     kr = m*pi/raio; % número de onda na direção radial
%     kr = ones(Ne_arco+1,1)*sqrt(k^2 - kz.^2); % número de onda na direção radial
    delta = 10^(-7); % tolerancia do truncamento   
    
    %%%% Calculo do Deslocamento e Velocidade %%%%%%%%%

    V = sqrt(-1)*omega*d; % calculo da velocidade normal

    v = V(:);   % vetor com os pontos nodais da velocidade 
    r = size(v,1);
    clear d V % limpando variáveis

     %%%%%%%%%%%%% Montando Matrizes de Influencia %%%%%%%%%%%%%%
    
    tic
%     [Ex,Ey,Ez] = coordxtr(edof,coord,dof,4);
    Ex = zeros(Ne,4);
    Ey = Ex; Ez = Ex;
    for e = 1:Ne
        Ex(e,:) = coord(conec(e,:),1);
        Ey(e,:) = coord(conec(e,:),2);
        Ez(e,:) = coord(conec(e,:),3);
    end
    
    H = zeros(r,r);
    G = zeros(r,r);
    for k=1:r
        for j=1:Ne
            [He,Ge] = bem_infl4q(coord(k,:),Ex(j,:),Ey(j,:),Ez(j,:),[omega c rho],-1);
            H = bem_assem(edof,H,He,k,j);
            G = bem_assem(edof,G,Ge,k,j);
        end
    end

    H = H + 1/2*eye(r);
    
    time = [];
    time = [time;toc];
%     coord2 = coord; coord2(:,3) = 10*ones(Np,1);
%     H_bar = zeros(r,r);
%     G_bar = zeros(r,r);
%     for k=1:r
%         for j=1:Ne
%             [He,Ge] = bem_infl4q(coord2(k,:),Ex(j,:),Ey(j,:),Ez(j,:),[omega c rho],-1);
%             H_bar = bem_assem(edof,H_bar,He,k,j);
%             G_bar = bem_assem(edof,G_bar,Ge,k,j);
%         end
%     end
  
    %%%%%%%%%% Calculo do Operador de Superfï¿½cie e da Intensidade Acï¿½stica %%%%%%%%%%%%%%
    
    
    tic
    R = H\G;      % operador de superfï¿½cie
%     R_bar = G_bar - H_bar*R;
    p = R*v;      % calculo da pressao na superficie
    I = 1/2*real(p.*conj(v)); %  Intensidade acustica convencional
    time = [time;toc];
    clear p % limpando variï¿½veis
    
   
    %%%%%%%%%% Calculo do Operador de Potencia %%%%%%%%%%%
    
    %%%%%%%%%%%%%%% Pontos e Pesos de Gauss %%%%%%%%%%%%%%%%%%%%%%%
 
    tic
    
    xi=[-0.8611363116 -0.3399810436 0.3399810436 0.8611363116];
    eta=[-0.8611363116 -0.3399810436 0.3399810436 0.8611363116];
    wxi=[0.3478548451 0.6521451549 0.6521451549 0.3478548451];
    weta=[0.3478548451 0.6521451549 0.6521451549 0.3478548451];

%     xi = [-sqrt(0.6) 0 -sqrt(0.6)];
%     eta = [-sqrt(0.6) 0 -sqrt(0.6)];;
%     wxi = [5/9 8/9 5/9];
%     weta = [5/9 8/9 5/9];
%     
%     xi = [-0.577350269189626 0.577350269189626];
%     eta = [-0.577350269189626 0.577350269189626];
%     wxi = [1 1];
%     weta = [1 1];
    
    Q = zeros(r,r);
    for e = 1:Ne
        Qe = zeros(Np,4);
        for i=1:size(xi,2)
            for j=1:size(eta,2)
                [N,A,Area] = Matrizelemento(e,edof,R,Ex(e,:),Ey(e,:),Ez(e,:),xi(i),eta(j));
                Qe = Qe + wxi(i)*weta(j)*(A*(N'*N))*(Area/4);
            end
        end
        Q(:,edof(e,2))= Q(:,edof(e,2))+ Qe(:,1);
        Q(:,edof(e,3))= Q(:,edof(e,3))+ Qe(:,2);
        Q(:,edof(e,4))= Q(:,edof(e,4))+ Qe(:,3);
        Q(:,edof(e,5))= Q(:,edof(e,5))+ Qe(:,4);
    end
    Q = 1/4*(transp(Q) + conj(Q));  % Operador de potï¿½ncia

%     Q = 1/2*real(Q);

%     Teta = zeros(r,r);
%     for e = 1:Ne
%         Tetae = zeros(4,4);
%         for i=1:size(xi,2)
%             for j=1:size(eta,2)
%                 [N,~,Area] = Matrizelemento(e,edof,R,Ex(e,:),Ey(e,:),Ez(e,:),xi(i),eta(j));
%                 Tetae = Tetae + wxi(i)*weta(j)*(N'*N)*(Area/4);
%             end
%         end
%         Teta(edof(e,2:5),edof(e,2:5))= Teta(edof(e,2:5),edof(e,2:5))+ Tetae;
%     end
%     Q = 0.5*real(R'*Teta);

    Pot = v'*Q*v;  % Cï¿½lculo da potï¿½ncia sonora sem o uso da Eigen-Decomposition
    NPS = 10*log10((Pot)/(10^(-12)));
    fprintf('Intensidade Acústica Numérica - NPS: %f \n',NPS)
    time = [time;toc];
    
    %%%%%%%%%% Decomposiï¿½ï¿½o em autovalores e sua ordenaï¿½ï¿½o %%%%%%%%%%%
    
    tic
    [V,D] = eig(Q);
    lambda = diag(D);
    
    clear Qe D  % limpando variï¿½veis
    
    [~,order] = sort(abs(lambda),'descend');    
    V = V(:,order);
    lambda = lambda(order);
    
    %%%%%%%%%% Aplicaï¿½ï¿½o do Critï¿½rio de Truncamento e Calculo da Velocidade ï¿½til %%%%%%%%%%%

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
            if abs((1 - soma/traco)*(1 - somas_pot(i)/Pot)) >= delta
                vu = vu + V(:,i)'*v*V(:,i);  % Velocidade util   
            end
            if cond == 0 && abs((1 - soma/traco)*(1 - somas_pot(i)/Pot)) < delta  
                    rc = i;  % posicao truncada
                    cond = 1;
            end
    end
    Pot_IU = somas_pot(rc); % Potï¿½ncia ï¿½til
     
    clear V lambda  % limpando variï¿½veis

    %%%%%%%%%% Calculo da Pressï¿½o ï¿½til e da Intensidade ï¿½til %%%%%%%%%%%

    pu = R*vu; % pressao util
    IU = 1/2*real(pu.*conj(vu)); % intensidade ï¿½til
    NPS_IU = 10*log10((Pot_IU)/(10^(-12))); % nï¿½vel de potï¿½ncia sonora ï¿½til
    fprintf('Intensidade Útil - NPS: %f \n',NPS_IU) 
    time = [time;toc];
    
   %%%%%%%%%% Decomposiï¿½ï¿½o em valores singulares e sua ordenaï¿½ï¿½o %%%%%%%%%%%
    
%     tic
%     coord2 = coord; coord2(:,3) = 10*ones(Np,1);
%     H_bar = zeros(r,r);
%     G_bar = zeros(r,r);
%     for k=1:r
%         for j=1:Ne
%             [He,Ge] = bem_infl4q(coord2(k,:),Ex(j,:),Ey(j,:),Ez(j,:),[omega c rho],-1);
%             H_bar = bem_assem(edof,H_bar,He,k,j);
%             G_bar = bem_assem(edof,G_bar,Ge,k,j);
%         end
%     end
%     R_bar = G_bar - H_bar*R;
%     [S,D,V] = svd(R_bar);
%     lambda = diag(D);
%     clear D
%             
%     vs = zeros(r,1); 
%     pot = zeros(r,1);
%     somas_pot2 = zeros(r,1);
%     aux = 0;
%     cond = 0;
%     delta = 10^(-7);
%     traco = sum(lambda);
%     sum_traco = 0;
%     for i=1:r
%             pot(i) = real(lambda(i)*(v.'*S(:,i))*(V(:,i)'*conj(v)));
%             somas_pot2(i) = aux + pot(i);
%             aux = somas_pot2(i);
%             sum_traco = sum_traco + lambda(i);
%             E = abs((traco - sum_traco)/traco)*abs((Pot - somas_pot2(i))/Pot);
%             if E >= delta
%                 vs = vs + V(:,i)'*v*V(:,i);  % Velocidade util   
%             end
%             if cond == 0 && (E < delta)
%                     rc2 = i;  % posicao truncada
%                     cond = 1;
%             end
%             
%     end
%     Pot_IS3 = somas_pot2(rc2); % Potencia supersônica numérica
%      
%     clear S V lambda  % limpando variï¿½veis

    %%%%%%%%%% Calculo da Pressï¿½o ï¿½til e da Intensidade ï¿½til %%%%%%%%%%%

%     ps = R*vs; % pressao supersonica
%     IS3 = 1/2*real(ps.*vs'.'); % intensidade supersonica segundo Magalhaes
%     NPS_IS3 = 10*log10((Pot_IS3)/(10^(-12))); % nivel de potencia sonora supersonica
%     fprintf('Intensidade Supersônica  segundo Magalhães - NPS: %f \n',NPS_IS3)
%     time = [time;toc];
    
    % Calculo da Intensidade Nao-Negativa Simétrica
    tic
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
    Z = R'*Teta;
    Zr = real(Z);
    [Psi,Lambda]=eig(Zr,Teta);
    beta = Psi*sqrt(Lambda)*Psi'*Teta*v;
    INN_Assim = 1/2*abs(beta).^2;
    Pot = 1/2*beta'*Teta*conj(beta);
    NPS_INN_Assim = 10*log10((Pot)/(10^(-12)));
    fprintf('Intensidade Não-Negativa Assimétrica - NPS: %f \n',NPS_INN_Assim)
    time = [time;toc];
    
    % Calculo da Intensidade Nao-Negativa Assimétrica
    tic
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
    Z = R'*Teta;
    Zr = real(Z);
    Zr_simm = (Zr' + Zr)/2;
    [Psi,Lambda]=eig(Zr_simm,Teta);
    beta = Psi*sqrt(Lambda)*Psi'*Teta*v;
    INN_Sim = 1/2*abs(beta).^2;
    Pot = 1/2*beta'*Teta*conj(beta);
    NPS_INN_Sim = 10*log10((Pot)/(10^(-12)));
    fprintf('Intensidade Não-Negativa Simétrica - NPS: %f \n',NPS_INN_Sim)
    time = [time;toc];
    
%     if exist(['placa times.txt'],'file')
%          dados = load('placa times.txt');
%          dados = [dados;frequencia,time'];
%     else
%          dados = [frequencia,time'];
%     end
%     save('resultados/silenciador times.txt','dados','-ascii','-double')
    
    clear R vu pu ps vs % limpando variï¿½veis
 
%%    Elementos de contorno de Galerkin 
%
%     H = zeros(r,r);
%     G = H;
%     tic
%     for e=1:Ne
%         for f=1:Ne
%             [He,Ge] = galerkin_bem_infl4q(Ex(e,:),Ey(e,:),Ez(e,:),Ex(f,:),Ey(f,:),Ez(f,:),[omega c rho],-1);
%             H(edof(e,2:5),edof(f,2:5))=H(edof(e,2:5),edof(f,2:5))+ He;
%             G(edof(e,2:5),edof(f,2:5))=G(edof(e,2:5),edof(f,2:5))+ Ge;
%         end
%     end
%     toc
%     H = H - 1/2*Teta;
%     figure
%     surf(isnan(H))
%     figure
%     surf(isnan(G))
%     R  = pinv(H)*G;      % operador de superfï¿½cie
%     p = R*v;      % calculo da pressao na superficie
%     I_GAL = 1/2*real(p.*conj(v)); %  Intensidade acustica convencional 
%     Pot = integral_dupla(Lx/Ne_x, Ly/Ne_y, I_GAL);
%     NPS = 10*log10((Pot)/(10^(-12)));
%     disp(NPS)
%     
%     toc % para o cronï¿½metro e exibe o tempo de operaï¿½ï¿½o
    
    %%%%%%%%%%%%%%%%%% Plotando Resultados %%%%%%%%%%%%%%%%%%%
    
    close all
    figure
    subplot(2,2,1)
%     plotagrandeza3(coord,conec,IS3) % plota a intensidade acustica numerica
%     title('IS via Magalhaes')
    plotagrandeza3(coord,conec,I) % plota a intensidade acustica numerica
    title('I')
    subplot(2,2,2)
    plotagrandeza3(coord,conec,IU) % plota a intensidade acustica util
    title('IU')
    subplot(2,2,3)
    plotagrandeza3(coord,conec,INN_Assim) % plota a intensidade nao negativa assimetrica
    title('INN Assim')
    subplot(2,2,4)
    plotagrandeza3(coord,conec,INN_Sim) % plota a intensidade nao negativa simetrica
    title('INN Sim')
    
    saveas(gcf,['resultados/silenciador grandezas - ',num2str(frequencia),'.png'])


    figure
    plotagrandeza3(coord,conec,imag(v))% plota a intensidade acústica
    saveas(gcf,['resultados/silenciador Velos - ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza3(coord,conec,I) % plota a intensidade acustica numerica
    saveas(gcf,['resultados/silenciador IA - ',num2str(frequencia),'.png'])
    
  
    figure
    plotagrandeza3(coord,conec,IU) % plota a intensidade acustica util
    saveas(gcf,['resultados/silenciador IU - ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza3(coord,conec,INN_Assim) % plota a intensidade nao negativa assimetrica
    saveas(gcf,['resultados/silenciador INN Assim - ',num2str(frequencia),'.png'])
   
     figure
    plotagrandeza3(coord,conec,INN_Sim) % plota a intensidade nao negativa simetrica
    saveas(gcf,['resultados/silenciador INN Sim - ',num2str(frequencia),'.png'])
    
%     figure
%     plotagrandeza3(coord,conec,IS3) % plota a intensidade acustica supersonica numerica
%     saveas(gcf,['resultados/silenciador IS via Magalhaes - ',num2str(frequencia),'.png'])

     fprintf('\n\n Operadores H e G: %f \n',time(1))
     fprintf('\n Operador R: %f \n',time(2))
     fprintf('\n Operador Q: %f \n',time(3))
     fprintf('\n IU: %f \n',time(4))
%      fprintf('\n IS via Magalhaes: %f \n',time(5))
     fprintf('\n INN ASim: %f \n',time(5))
     fprintf('\n INN Sim: %f \n',time(6))

%     plotagrandeza(coord,conec,I_GAL) % plota a intensidade acústica numérica via Galerkin BEM
%     title('intensidade acustica numerica via Galerkin BEM')
%     colorbar
%     saveas(gcf,['intensidade acustica numerica via Galerkin BEM - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
   
%     subplot (3,3,8)
%     plotagrandeza(coord,conec,IU2) % plota a intensidade acustica supersonica via Magalhaes
%     title('IS Magalhaes')
%     colorbar
%     saveas(gcf,['intensidade acustica nao negativa - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
%     P1 = integral_dupla(Lx/Ne_x, Ly/Ne_y, reshape(I,Ne_x+1,Ne_y+1));
%     P1 = 10*log10((P1)/(10^(-12)));
%     P2 = integral_dupla(Lx/Ne_x, Ly/Ne_y, reshape(IU,Ne_x+1,Ne_y+1));
%     P2 = 10*log10((P2)/(10^(-12)));
%     P3 = integral_dupla(Lx/Ne_x, Ly/Ne_y, reshape(INN,Ne_x+1,Ne_y+1));
%     P3 = 10*log10((P3)/(10^(-12)));
%     
%     P4 = sum(I)*A/Ne;
%     P4 = 10*log10((P4)/(10^(-12)));
%     P5 = sum(IU)*A/Ne;
%     P5 = 10*log10((P5)/(10^(-12)));
%     P6 = sum(INN)*A/Ne;
%     P6 = 10*log10((P6)/(10^(-12)));
    
%     fprintf('\n %f \n %f \n %f \n\n %f \n %f \n %f \n',P1,P2,P3,P4,P5,P6)
    
%     dados = [I, IU, IS, I2(:), IS2(:), INN];
%     save( ['resultados/placa grandezas - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.txt'],'dados','-ascii','-double')
%     
%     sigma_IU = Pearson(IS2(:),IU);
%     fprintf('\n Pearson IU-IS: %f',sigma_IU)
%     sigma_INN = Pearson(IS2(:),INN);
%     fprintf('\n Pearson INN-IS: %f',sigma_INN)
%     sigma = Pearson(IU,INN);
%     fprintf('\n Pearson IU-INN: %f \n',sigma)
% 
%     
%     if exist(['resultados/placa medidas - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'.txt'],'file')
%          dados2 = load( ['placa medidas - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'.txt']);
%          dados2 = [dados2;[frequencia,NPS_IS2,NPS_IU,NPS_INN,sigma_IU,sigma_INN]];
%     else
%         dados2 = [frequencia,ipsilon,NPS2,NPS_IS2,NPS,NPS_IU,NPS_INN,sigma_IU,sigma_INN];
%     end
%     save(['resultados/placa medidas - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'.txt'],'dados2','-ascii','-double')
% 
%     % Plota a contribuiï¿½ï¿½o dos valores prï¿½prios   

    figure
    plot( (1:r)',somas_pot, rc*ones(1000, 1), linspace(min(somas_pot),somas_pot(rc), 1000 ),'b--')
    grid
    xlabel('Valores próprios de velocidade considerados')
    ylabel('Potência Sonora (Watts)')
    saveas(gcf,['resultados/silenciador convergencia potencia - ',num2str(frequencia),'.png'])
    aux = [rc;somas_pot];
    save(['resultados/silenciador convergencia potencia - ',num2str(frequencia),'.txt'],'aux','-ascii','-double')

%     figure
%     plot( (1:r)',somas_pot2,'r', rc2*ones(1000, 1), linspace(min(somas_pot2),somas_pot2(rc2), 1000 ),'r--')
%     grid
%     xlabel('Valores singulares de velocidade considerados')
%     ylabel('Potência Sonora (Watts)')
    
%     close all
    
end
function calc_IU_placa(L,modos,Ne,frequencia)
    %%%%%%%%%%%%% Designando vari�veis de entrada %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all;
    L = L(:);
    modos = modos(:);
    Ne = Ne(:);
    
    Lx = L(1);      % Comprimento da placa na dire��o x
    Ly = L(2);      % Comprimento da placa na dire��o y
%     Lz = L(3);
    m = modos(1);   % modo de vibra��o na dire��o x
    n = modos(2);   % modo de vibra��o na dire��o y
    Ne_x = Ne(1);  % n�mero de elementos na dire��o x
    Ne_y = Ne(2);  % n�mero de elementos na dire��o y
    
    %%%%%%%%%%%%% Constantes  F�sicas  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = 343;  % velocidade do som
    rho = 1.29; % massa espec�ficado meio (ar)
        
    %%%%%%% Calculo de vari�veis importantes na placa %%%%%%%%
    
%     tic % inicia o cronometro
    
    x = linspace(0,Lx,Ne_x + 1); % abcissas dos pontos  na malha
    y = linspace(0,Ly,Ne_y + 1); % ordenadas dos pontos na malha
    
    
    A = Lx*Ly;  % �rea superficial da placa
    Ne = Ne_x*Ne_y; % n�mero total de elementos na malha
    Np = (Ne_x + 1)*(Ne_y + 1);
   
    kx = m*pi/Lx; % n�mero de onda na dire��o x
    ky = n*pi/Ly; % n�mero de onda na dire��o y
    kf = sqrt(kx^2+ky^2); % n�mero de onda livre
%     frequencia = ipsilon*kf*c/(2*pi); % frequencia de excita��o
    omega = 2*pi*frequencia; %  frequencia angular
    k = omega/c; % n�mero de onda
    ipsilon = k/kf; % frequencia adimensional
    fprintf('A frequ�ncia adimensional ser�: %f \n',ipsilon)
    delta = 10^(-7); % tolerancia do truncamento

    
    %%%%%%%%%%%%% Determina��o do modo de vibra��o %%%%%%%%%%%%%%
    fprintf('A placa possui um Modo de ')
    if  kx^2 + ky^2 <= k^2
        fprintf('Superficie \n')
    elseif (kx > k) && (ky > k)
        fprintf('Canto \n')
    elseif (kx > k && ky < k)||(kx < k && ky > k)
        fprintf('Borda \n')
    else
        fprintf('Transicao \n')
    end
    
    
    %%%% Calculo do Deslocamento e Velocidade %%%%%%%%%
    
    D = zeros(Ne_x + 1,Ne_y + 1);
    for i = 1:Ne_x + 1
        for j =1:Ne_y + 1
            D(i,j) = 2/sqrt(A)*sin(kx*x(i))*sin(ky*y(j))*10^(-6);   % Calculo do deslocamento na superf�cie da placa pela formula anal�tica
        end
    end 
    V = sqrt(-1)*omega*D; % calculo da velocidade normal

    v = V(:);   % vetor com os pontos nodais da velocidade 
    r = size(v,1);
    
    kx1 = ifftshift( (2*pi/Lx)*(-Ne_x/2:Ne_x/2)  ); % n�meros de onda da dire��o x 
    ky1 = ifftshift( (2*pi/Ly)*(-Ne_y/2:Ne_y/2)  ); % n�meros de onda da dire��o y
    kx1 = (2*pi/Lx)*[0:((Ne_x+1)/2-1) (-(Ne_x+1)/2):-1];
    ky1 = (2*pi/Ly)*[0:((Ne_y+1)/2-1) (-(Ne_y+1)/2):-1];

    kz = zeros(Ne_x + 1,Ne_y + 1);    
    for i = 1:Ne_x + 1
        for j = 1:Ne_y + 1
            kz(i,j) = sqrt(k^2 - kx1(i)^2 - ky1(j)^2);
        end
    end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    tV = fft2(V); % velocidade angular
    tP = rho*c*k*tV./kz;  % pressao angular
    P = ifft2(tP);   % pressao 
    I2 = 1/2*real(P.*V'.'); % intensidade ac�stica
    Pot2 = integral_dupla(Lx/Ne_x, Ly/Ne_y, I2); % Pot�ncia sonora
    NPS2 = 10*log10((Pot2)/(10^(-12)));
    fprintf('Intensidade Ac�stica Convencional - NPS: %f \n',NPS2)
    
    %%%%%%%%%%%%%%%%% Calculo das Grandezas Supers�nicas %%%%%%%%%%%%%%%%
     
    clear D V P kz
    tPS = zeros(Ne_x + 1,Ne_y + 1);  % press�o angular supers�nica
    tVS = zeros(Ne_x + 1,Ne_y + 1);  % velocidade angular supers�nica
    for i = 1:Ne_x + 1
        for j =1:Ne_y + 1
            if  kx1(i)^2 + ky1(j)^2 <= k^2
                tPS(i,j) = tP(i,j);
                tVS(i,j) = tV(i,j);
            end
        end
    end
    VS = ifft2(tVS);   % velocidade supersonica
    PS = ifft2(tPS);   % pressao 
    IS2 = 1/2*real(PS.*conj(VS)); % intensidade ac�stica
    Pot2 = integral_dupla(Lx/Ne_x, Ly/Ne_y, IS2); % Pot�ncia sonora
    NPS_IS2 = 10*log10((Pot2)/(10^(-12)));
    fprintf('Intensidade Supers�nica - NPS: %f \n',NPS_IS2)
    
%     vs = VS(:);
%     clear VS tVS  % limpando vari�veis


    %%%%%%%%%%%%% Montando Matrizes de Influencia %%%%%%%%%%%%%%
    
    tic
    [coord,edof,conec,dof] = topologia_placa(Lx,Ly,Ne_x,Ne_y);
    [Ex,Ey,Ez] = coordxtr(edof,coord,dof,4);
%     Ex = zeros(Ne,4);
%     Ey = Ex; Ez = Ex;
%     for e = 1:Ne
%         Ex(e,:) = coord(conec(e,:),1);
%         Ey(e,:) = coord(conec(e,:),2);
%         Ez(e,:) = coord(conec(e,:),3);
%     end
    
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

  
    %%%%%%%%%% Calculo do Operador de Superf�cie e da Intensidade Ac�stica %%%%%%%%%%%%%%
    
    
    tic
    R = H\G;      % operador de superf�cie
    p = R*v;      % calculo da pressao na superficie
    I = 1/2*real(p.*conj(v)); %  Intensidade acustica convencional
    time = [time;toc];
    clear p % limpando vari�veis
    
   
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
    Q = 1/4*(transp(Q) + conj(Q));  % Operador de pot�ncia

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

    Pot = v'*Q*v;  % C�lculo da pot�ncia sonora sem o uso da Eigen-Decomposition
    NPS = 10*log10((Pot)/(10^(-12)));
    fprintf('Intensidade Ac�stica Num�rica - NPS: %f \n',NPS)
    time = [time;toc];
    
    %%%%%%%%%% Decomposi��o em autovalores e sua ordena��o %%%%%%%%%%%
    
    tic
    [V,D] = eig(Q);
    lambda = diag(D);
    
%     pot = zeros(r,1);
%     for i = 1:r
%         pot(i) = real(lambda(i)*abs(v'*V(:,i))^2);  
%     end
%     [pot,order] = sort(pot);
%     V = V(:,order);
%     sum_pot_neg = abs(sum(pot(pot < 0)));
%     pot_pos = pot(pot >= 0);
%     psum_pot_pos = cumsum(pot_pos);
%     V = V(:,pot >= 0);
%     ind = find(abs(psum_pot_pos >= sum_pot_neg));
%     V = V(:,ind);
%     rc = size(ind,1);
%     vu = zeros(r,1); 
%     Pot_IU = sum(pot_pos(ind));
%     for i = 1:rc
%         vu = vu + V(:,i)'*v*V(:,i);
%     end

    %  
    clear Qe D  % limpando vari�veis
    
    [~,order] = sort(abs(lambda),'descend');    
    V = V(:,order);
    lambda = lambda(order);
    
    %%%%%%%%%% Aplica��o do Crit�rio de Truncamento e Calculo da Velocidade �til %%%%%%%%%%%

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
    Pot_IU = somas_pot(rc); % Pot�ncia �til
    % 
    clear V lambda  % limpando vari�veis

    %%%%%%%%%% Calculo da Press�o �til e da Intensidade �til %%%%%%%%%%%

    pu = R*vu; % pressao util
    IU = 1/2*real(pu.*conj(vu)); % intensidade �til
    NPS_IU = 10*log10((Pot_IU)/(10^(-12))); % n�vel de pot�ncia sonora �til
    fprintf('Intensidade �til - NPS: %f \n',NPS_IU) 
    time = [time;toc];
    
    %%%%%%%%%% Decomposi��o em valores singulares e sua ordena��o %%%%%%%%%%%
    
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
%     Pot_IS3 = somas_pot2(rc2); % Potencia supers�nica num�rica
%      
%     clear S V lambda  % limpando vari�veis

    %%%%%%%%%% Calculo da Press�o �til e da Intensidade �til %%%%%%%%%%%

%     ps = R*vs; % pressao supersonica
%     IS3 = 1/2*real(ps.*vs'.'); % intensidade supersonica segundo Magalhaes
%     NPS_IS3 = 10*log10((Pot_IS3)/(10^(-12))); % nivel de potencia sonora supersonica
%     fprintf('Intensidade Supers�nica  segundo Magalh�es - NPS: %f \n',NPS_IS3)
%     time = [time;toc];
    
    % Calculo da Intensidade Nao-Negativa Sim�trica
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
    fprintf('Intensidade N�o-Negativa Assim�trica - NPS: %f \n',NPS_INN_Assim)
    time = [time;toc];
    
    % Calculo da Intensidade Nao-Negativa Assim�trica
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
    fprintf('Intensidade N�o-Negativa Sim�trica - NPS: %f \n',NPS_INN_Sim)
    time = [time;toc];
    
%     if exist(['placa times.txt'],'file')
%          dados = load('placa times.txt');
%          dados = [dados;frequencia,time'];
%     else
%          dados = [frequencia,time'];
%     end
%     save('resultados/placa times.txt','dados','-ascii','-double')
    
    clear R vu pu ps vs % limpando vari�veis
 
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
%     R  = pinv(H)*G;      % operador de superf�cie
%     p = R*v;      % calculo da pressao na superficie
%     I_GAL = 1/2*real(p.*conj(v)); %  Intensidade acustica convencional 
%     Pot = integral_dupla(Lx/Ne_x, Ly/Ne_y, I_GAL);
%     NPS = 10*log10((Pot)/(10^(-12)));
%     disp(NPS)
%     
%     toc % para o cron�metro e exibe o tempo de opera��o
    
    %%%%%%%%%%%%%%%%%% Plotando Resultados %%%%%%%%%%%%%%%%%%%
    
    close all
    figure
    
    subplot(3,3,1)
    plotagrandeza(coord,conec,imag(v))% plota a intensidade ac�stica
    title('velocidade normal')
    
    subplot(3,3,3)
    plotagrandeza(coord,conec,I) % plota a intensidade acustica numerica
    title('IA numerica')
    subplot(3,3,5)
    plotagrandeza(coord,conec,IU) % plota a intensidade acustica util
    title('IU')
    subplot(3,3,2)
    plotagrandeza(coord,conec,I2(:)) % plota a intensidade acustica analitica
    title('IA analitica')
    subplot(3,3,4)
    plotagrandeza(coord,conec,IS2(:)) % plota a intensidade acustica supersonica analitica
    title('IS analitica')
    subplot(3,3,6)
    plotagrandeza(coord,conec,INN_Assim) % plota a intensidade acustica supersonica analitica
    title('INN Assim')
    subplot(3,3,7)
    plotagrandeza(coord,conec,INN_Sim) % plota a intensidade acustica supersonica analitica
    title('INN Sim')
%     plotagrandeza(coord,conec,IS3) % plota a intensidade acustica supersonica analitica
%     title('IS via Magalhaes')
    saveas(gcf,['resultados/placa grandezas - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])


    figure
    plotagrandeza(coord,conec,imag(v))% plota a intensidade ac�stica
    saveas(gcf,['resultados/placa Velos - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza(coord,conec,I) % plota a intensidade acustica numerica
    saveas(gcf,['resultados/placa IA - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza(coord,conec,IS2(:)) % plota a intensidade acustica supersonica analitica
    saveas(gcf,['resultados/placa IS analitica - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza(coord,conec,IU) % plota a intensidade acustica util
    saveas(gcf,['resultados/placa IU - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza(coord,conec,INN_Assim) % plota a intensidade acustica util
    saveas(gcf,['resultados/placa INN Assim - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
   
     figure
    plotagrandeza(coord,conec,INN_Sim) % plota a intensidade acustica util
    saveas(gcf,['resultados/placa INN Sim - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
   
    figure
    plotagrandeza(coord,conec,I2(:)) % plota a intensidade acustica analitica
    saveas(gcf,['resultados/placa IA analitica - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
%     figure
%     plotagrandeza(coord,conec,IS3) % plota a intensidade acustica analitica
%     saveas(gcf,['resultados/placa IS via Magalhaes - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
%     3D figures
    
    figure
    plotagrandeza2(x,y,imag(reshape(v,Ne_x+1,Ne_y+1)))% plota a intensidade ac�stica
    saveas(gcf,['resultados/placa 3D Velos - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza2(x,y,reshape(I,Ne_x+1,Ne_y+1)) % plota a intensidade acustica numerica
    saveas(gcf,['resultados/placa 3D IA - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza2(x,y,IS2) % plota a intensidade acustica supersonica analitica
    saveas(gcf,['resultados/placa 3D IS analitica - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza2(x,y,reshape(IU,Ne_x+1,Ne_y+1)) % plota a intensidade acustica util
    saveas(gcf,['resultados/placa 3D IU - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
    figure
    plotagrandeza2(x,y,reshape(INN_Assim,Ne_x+1,Ne_y+1)) % plota a intensidade acustica nao negativa
    saveas(gcf,['resultados/placa 3D INN Assim - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])

    figure
    plotagrandeza2(x,y,reshape(INN_Sim,Ne_x+1,Ne_y+1)) % plota a intensidade acustica nao negativa
    saveas(gcf,['resultados/placa 3D INN Sim - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])

%     figure
%     plotagrandeza2(x,y,reshape(IS3,Ne_x+1,Ne_y+1)) % plota a intensidade acustica nao negativa
%     saveas(gcf,['resultados/placa 3D IS via Magalhaes - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
%     
    figure
    plotagrandeza2(x,y,I2) % plota a intensidade acustica analitica
    saveas(gcf,['placa 3D IA analitica - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
    
     fprintf('\n\n Operadores H e G: %f \n',time(1))
     fprintf('\n Operador R: %f \n',time(2))
     fprintf('\n Operador Q: %f \n',time(3))
     fprintf('\n IU: %f \n',time(4))
%      fprintf('\n IS via Magalhaes: %f \n',time(5))
     fprintf('\n INN ASim: %f \n',time(5))
     fprintf('\n INN Sim: %f \n',time(6))
     

%     plotagrandeza(coord,conec,I_GAL) % plota a intensidade ac�stica num�rica via Galerkin BEM
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
%     % Plota a contribui��o dos valores pr�prios   

%     figure
%     plot( (1:r)',somas_pot, rc*ones(1000, 1), linspace(min(somas_pot),somas_pot(rc), 1000 ),'b--')
%     grid
%     xlabel('Valores pr�prios de velocidade considerados')
%     ylabel('Pot�ncia Sonora (Watts)')
%     saveas(gcf,['resultados/placa convergencia potencia - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.png'])
%     aux = [rc;somas_pot];
%     save(['resultados/placa convergencia potencia - ',num2str(Lx),'-',num2str(Ly),' - ',num2str(m),'-',num2str(n),' - ',num2str(Ne_x),'-',num2str(Ne_y),'- ',num2str(frequencia),'.txt'],'aux','-ascii','-double')

%     figure
%     plot( (1:r)',somas_pot2,'r', rc2*ones(1000, 1), linspace(min(somas_pot2),somas_pot2(rc2), 1000 ),'r--')
%     grid
%     xlabel('Valores singulares de velocidade considerados')
%     ylabel('Pot�ncia Sonora (Watts)')
    
%     close all

end
% L'utilité de chaque section est précisée en commentaire de chaque partie.
% pour etre utilisée, chaque partie doit etre décommentée suivant les
% besoins.

%% Données
% tfin    = 179.428  % temps total de simualation (s)
% mass    = 1.5e-1   % masse en (kg)
% q       = 1.e-4    % charge en (C)
% L       = 2.e-1    % longueur du pendule (m)
% g       = 9.81e0   % acceleration de gravite (m/s^2)
% k       = 0        % constante de la force du frottement (kg/m)
% E0      = 3.e3     % intensite du champ electrique (V/m)
% omega   = 7.0035   % vitesse angulaire champ electrique (r/s)
% phi     = 0.e0     % phase champ electrique
% theta0  = 3        % angle initiale (r)
% thetaDot0 = 0      % vitesse angulaire initiale (r/s)
% nsteps = 10000 

%tfin = 9395621.2017;     % recherchée
tfin = 9400753.64319;      % marcelo 
%tfin = 9400964.925;


%% Omega : si je veux faire varier omega (b,ii) avec linespace:


% repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
% executable = 'Exercice3_2020_student'; % Nom de l'executable (NB: ajouter .exe sous Windows)
% input = 'config.in'; % Nom du fichier d'entree de base MODIFIER SELON VOS BESOINS
% 
% nsimul = 25; % Nombre de simulations a faire
% 
% omega = linspace(sqrt(g/L)-0.27,sqrt(g/L)+0.22,nsimul); % Nombre d'iterations entier de 10^2 a 10^4  MODIFIER SELON VOS BESOINS
% 
% paramstr = 'omega'; % Nom du parametre a scanner  MODIFIER SELON VOS BESOINS
% param = omega; % Valeurs du parametre a scanner  MODIFIER SELON VOS BESOINS
% 
% output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
% for i = 1:nsimul
%     output{i} = [paramstr, '=', num2str(param(i)), '.out']
%     % Execution du programme en lui envoyant la valeur a scanner en argument
%     cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i})
%     disp(cmd)
%     system(cmd);
% end


%% theta0 : si je veux faire varier theta0 (d,iv) avec linespace:


% repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
% executable = 'Exercice3_2020_student'; % Nom de l'executable (NB: ajouter .exe sous Windows)
% input = 'config.in'; % Nom du fichier d'entree de base MODIFIER SELON VOS BESOINS
% 
% nsimul = 2; % Nombre de simulations a faire
% 
% theta0 = [2,2+1.0e-8];
% %theta0 = linspace(3,3-1.0e-8,nsimul);
% 
% %theta0 = logspace(2.356194-4.39999999998e-5,2.356194,nsimul); % Nombre d'iterations entier de 10^2 a 10^4  MODIFIER SELON VOS BESOINS
% 
% paramstr = 'theta0'; % Nom du parametre a scanner  MODIFIER SELON VOS BESOINS
% param = theta0; % Valeurs du parametre a scanner  MODIFIER SELON VOS BESOINS
% 
% output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
% for i = 1:nsimul
%     output{i} = [paramstr, '=', num2str(param(i)), num2str(i), '.out']
%     % Execution du programme en lui envoyant la valeur a scanner en argument
%     cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i})
%     disp(cmd)
%     system(cmd);
% end



%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice4_Harfouche_Frybes'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'config.in'; % Nom du fichier d'entree de base MODIFIER SELON VOS BESOINS

nsimul = 10; % Nombre de simulations a faire

nsteps = round(logspace(4,5,nsimul)); % Nombre d'iterations entier de 10^2 a 10^4  MODIFIER SELON VOS BESOINS

paramstr = 'nsteps'; % Nom du parametre a scanner  MODIFIER SELON VOS BESOINS
param = nsteps; % Valeurs du parametre a scanner  MODIFIER SELON VOS BESOINS

% Simulations %% 
%%%%%%%%%%%%%%%%
% Lance une serie de simulations (= executions du code C++)
% Normalement, on ne devrait pas avoir besoin de modifier cette partie

output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out']
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i})
    disp(cmd)
    system(cmd);
end

%% 1ere partie: tracage des trajectoires (a-i) : 
%%%%%%%%%%%%%

% La position 

for i=1:nsimul
    data= load(output{i});
    t = [zeros,nsteps(i)];
    
    x1= [zeros,nsteps(i)];
    y1= [zeros,nsteps(i)];
    
    x2= [zeros,nsteps(i)];
    y2= [zeros,nsteps(i)];
    
    x3= [zeros,nsteps(i)];
    y3= [zeros,nsteps(i)];
    
    energy=[zeros,nsteps(i)];
    
    for j=1:size(data)
        t(j)=data(j,1);
        
        x1(j)=data(j,2);    
        y1(j)=data(j,3);
        
        x2(j)=data(j,4);    
        y2(j)=data(j,5);
        
        x3(j)=data(j,6);    
        y3(j)=data(j,7);
        
        energy(j)=abs(data(j,14)-data(1,14));
        
    end        
end



    figure
    plot(x1, y1)
    set(gca,'FontSize',20)
    xlabel('x [m]','FontSize',18)
    ylabel('y [m]','FontSize',18)
    
    figure
    plot(x2, y2)
    set(gca,'FontSize',20)
    xlabel('x [m]','FontSize',18)
    ylabel('y [m]','FontSize',18)
    
    figure
    plot(t, energy)
    set(gca,'FontSize',20)
    xlabel('t [s]','FontSize',18)
    ylabel('Mechanical energy \Delta E_{mec} [J]','FontSize',18)
 
    
%     figure
%     plot(x3,y3)







% dessin de l'évloutuin théorique de l'angle theta (correspondant au dernier simul de nsimul):
% theorique = [zeros,nsteps(3)];
% for k=1:size(data)
%    theta_theorique(k)=theta_ana(g,L,data(k,1),theta0,thetaDot0); 
% end
% 
% 
% figure
% plot(t,theta,t,theta_theorique)


% La vitesse  
% 
% for i=1:nsimul
%     data= load(output{i});
%     t= [zeros,nsteps(i)];
%     thetaDot= [zeros,nsteps(i)];
%     for j=1:size(data)
%         t(j)=data(j,1);    % theta en fct du temps
%         thetaDot(j)=data(j,3);
%     end    
%  
%     
% end
% 
% figure
% plot(t,thetaDot)

% dessin de l'évloutuin théorique de la vitesse angulaire theta point (correspondant au dernier simul de nsimul):
% thetaDot_theorique = [zeros,nsteps(nsimul)];
% for k=1:size(data)
%    thetaDot_theorique(k)=thetaDot_ana(g,L,data(k,1),theta0,thetaDot0); 
% end


% figure
% plot(t,thetaDot)

% figure
% plot(t,thetaDot_theorique)



%% Etude de convergence de la position finale (a.ii) en connaissant la solution analytique : 

% pour runge kutta avec nsteps

erreurpos = [zeros,nsimul];
erreurvit = [zeros,nsimul];
Deltat = [zeros,nsimul];
fin = [zeros,nsimul];
for i=1:nsimul
    data= load(output{i});
    Deltat(i)= tfin/nsteps(i);
    %erreur(i)= abs( data(nsteps(i)+1,4)-data(1,4) ); %x
    %erreur(i)= abs( data(nsteps(i)+1,5)- data(1,5)); %y
    %erreur(i)= abs( sqrt( data(nsteps(i)+1,4)*data(nsteps(i)+1,4)+data(nsteps(i)+1,5)*data(nsteps(i)+1,5) )-sqrt( data(1,4)*data(1,4)+data(1,5)*data(1,5) ) ); %norme
    erreurpos(i)= sqrt ( (data(end,4)-data(1,4))*(data(end,4)-data(1,4)) + (data(end,5)-data(1,5))*(data(end,5)-data(1,5)));
    %erreurvit(i)= sqrt ( (data(nsteps(i)+1,10)-data(1,10))*(data(nsteps(i)+1,10)-data(1,10)) + (data(nsteps(i)+1,11)-data(1,11))*(data(nsteps(i)+1,11)-data(1,11)));
    
     %x = interp1(data(end-4:end,1),data(end-4:end,4),tfin);
     %y = interp1(data(end-4:end,1),data(end-4:end,5),tfin);
     %erreurpos(i)= sqrt ( (x-data(1,4))*(x-data(1,4)) + (y-data(1,5))*(y-data(1,5)));

end
% 
% figure
% scatter(Deltat.^4,erreur,'r+')

figure('Name','fin')
loglog(Deltat, erreurpos, 'r+', 'linewidth' , 1.5 )
set(gca,'FontSize',16)
xlabel('\Delta t [s]','FontSize',18)
ylabel('Error on final position [m]','FontSize',18)
ppos = polyfit (log10(Deltat),log10(erreurpos),1)
hold on
plot(Deltat,10.^(ppos(2)).*Deltat.^(ppos(1)),'r')
grid on

% 
% 
% figure('Name','fin')
% loglog(Deltat, erreurvit, 'r+', 'linewidth' , 1.5 )
% set(gca,'FontSize',16)
% xlabel('\Delta t [s]','FontSize',18)
% ylabel('Error on final velocity [m.s^{-1}]','FontSize',18)
% pvit = polyfit (log10(Deltat),log10(erreurvit),1)
% hold on
% plot(Deltat,10.^(pvit(2)).*Deltat.^(pvit(1)),'r')
% grid on


% pour runge kutta avec un while t<tfin

% erreurpos = [zeros,nsimul];
% erreurvit = [zeros,nsimul];
% Deltat = [zeros,nsimul];
% fin = [zeros,nsimul];
% for i=1:nsimul
%     data= load(output{i});
%     
%     Deltat(i)= tfin/nsteps(i);
%     %erreur(i)= abs( data(nsteps(i)+1         ,4)-data(1,4) ); %x
%     %erreur(i)= abs( data(nsteps(i)+1,5)- data(1,5)); %y
%     %erreur(i)= abs( sqrt( data(nsteps(i)+1,4)*data(nsteps(i)+1,4)+data(nsteps(i)+1,5)*data(nsteps(i)+1,5) )-sqrt( data(1,4)*data(1,4)+data(1,5)*data(1,5) ) ); %norme
%     erreurpos(i)= sqrt ( (data(end,4)-data(1,4))*(data(end,4)-data(1,4)) + (data(end,5)-data(1,5))*(data(end,5)-data(1,5)));
%     erreurvit(i)= sqrt ( (data(end,10)-data(1,10))*(data(end,10)-data(1,10)) + (data(end,11)-data(1,11))*(data(end,11)-data(1,11)));
%     
% end
% 
% 
% figure('Name','fin')
% loglog(Deltat, erreurpos, 'r+', 'linewidth' , 1.5 )
% set(gca,'FontSize',16)
% xlabel('\Delta t [s]','FontSize',18)
% ylabel('Error on final position [m]','FontSize',18)
% ppos = polyfit (log10(Deltat),log10(erreurpos),1)
% hold on
% plot(Deltat,10.^(ppos(2)).*Deltat.^(ppos(1)),'r')
% grid on

% 
% 
% figure('Name','fin')
% loglog(Deltat, erreurvit, 'r+', 'linewidth' , 1.5 )
% set(gca,'FontSize',16)
% xlabel('\Delta t [s]','FontSize',18)
% ylabel('Error on final velocity [m.s^{-1}]','FontSize',18)
% pvit = polyfit (log10(Deltat),log10(erreurvit),1)
% hold on
% plot(Deltat,10.^(pvit(2)).*Deltat.^(pvit(1)),'r')
% grid on



%% Etude de convergence sur les position et vitesse maximales (question a) :

% calcul des erreurs maximales :

% errtheta = [zeros,nsimul];
% errthetadot = [zeros,nsimul];
% dt = [zeros,nsimul]
% for i = 1:nsimul % Parcours des resultats de toutes les simulations
%     data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
%     dt(i) = data(2,1)-data(1,1);
%     err1 = 0;
%     err2 = 0;
% %     theta = [zeros,nsteps(i)];
% %     t = [zeros,nsteps(i)];
%     for  j = 1:nsteps(i)
%         if abs(theta_ana(g,L,data(j,1),theta0,thetaDot0)- data(j,2))> err1
%            err1 = abs(theta_ana(g,L,data(j,1),theta0,thetaDot0)- data(j,2));
%         end
%         
%         if abs(thetaDot_ana(g,L,data(j,1),theta0,thetaDot0)- data(j,3)) > err2
%            err2 = abs(thetaDot_ana(g,L,data(j,1),theta0,thetaDot0)- data(j,3));
%         end
%     end 
% errtheta(i) = err1;
% errthetadot(i) = err2;
% end
% 
% 
% % plot des erreurs maximales : 
% 
% % max de la position
% figure
% loglog(dt, errtheta, 'r+')
% xlabel('\Delta t [s]')
% ylabel('Maximum de l''erreur sur la position (\Deltaz) [m\cdots^{-1}]')
% p= polyfit (log10(dt),log10(errtheta),1);
% hold on
% p = plot(dt,10.^(p(2)).*dt.^(p(1)),'r')
% grid on
% 
% % max de la vitesse
% figure
% loglog(dt, errthetadot, 'r+')
% xlabel('\Delta t [s]')
% ylabel('Maximum de l''erreur sur la vitesse (\Deltav_z) [m\cdots^{-1}]')
% p= polyfit (log10(dt),log10(errthetadot),1);
% hold on
% p = plot(dt,10.^(p(2)).*dt.^(p(1)),'r')
% grid on




%% Plot de l'énergie mécanique (question a): 

% for i=1:nsimul
%     data= load(output{i});
%     t= [zeros,nsteps(i)];
%     Delta_energie= [zeros,nsteps(i)];
%     for j=1:size(data)
%         t(j)=data(j,1);    
%         Delta_energie(j)=data(j,4)-data(1,4);    % energie en fct du temps
%     end    
%  
%     figure
%     plot(t,Delta_energie)
%     
% end


%% %% Etude de convergence de la position finale (en ignorant la solution analytique) : 

% erreurposfin = [zeros,nsimul]
% Deltat = [zeros,nsimul]
% for i=1:nsimul
%     data= load(output{i});
%     Deltat(i)= tfin*(1/nsteps(i));
%     erreurposfin(i)= data(nsteps(i)+1,2);
% end
% Deltat2 = Deltat.^2;
% 
% figure
% plot(Deltat2, erreurposfin, 'r+')
% xlabel('\Delta t [s]')
% ylabel('Maximum de l''erreur sur la vitesse (\Deltav_z) [m\cdots^{-1}]')


%% Vérification du théorème de l'energie mécanique b)i:  (le nsteps pris est le dernier nsteps du nsimul)
  
% Vérification du théorème de l'énergie mécanique:  b-i
% dt= [zeros,nsteps(nsimul)];
% 
% DEmec = [zeros,nsteps(nsimul)];
% data= load(output{nsimul});
% for i=1:nsteps(nsimul)
%    dt(i)= data(i,1);
%    DEmec(i) = (data(i+1,4)-data(i,4))/(data(i+1,1)-data(i,1));
% end
% 
% figure
% plot(dt, DEmec, 'r+')
% xlabel('t [s]')
% ylabel('Energie mecanique [J]')



% Vérification du théorème de l'énergie mécanique:  b-i
% Puiss = [zeros,nsteps(nsimul)];
% for i=1:nsteps(nsimul)
%    Puiss(i) = data(i,5); % P_ana(g,L,k,q,E0,omega,data(i,3),data(i,1));
% end
% 
% figure
% plot(dt, Puiss, 'r+')
% xlabel('t [s]')
% ylabel('Puissance [J]')



%% maximum de lenergie mécanique b-ii:


% Emec_max =  [zeros,nsimul];
% Delt =  [zeros,nsimul];
% for i = 1:nsimul % Parcours des resultats de toutes les simulations
%     data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
%     val= 0;
%     Delt(i)= tfin/nsteps(i);
%     for  j = 1:nsteps
%         if abs(data(j,14)-data(1,14)) > val
%            val = abs(data(j,14)-data(1,14));
%         end
%     end 
%     Emec_max(i)= val;
% end


% figure
% plot(Delt, Emec_max , 'r+')
% xlabel('\Delta t [s]')
% ylabel('Energie maximale \Delta E_{mec,max} [J]')

% max de la position
% figure('Name','Emec')
% loglog(Delt, Emec_max, 'r+', 'linewidth' , 1.5 )
% set(gca,'FontSize',18)
% xlabel('\Delta t [s]','FontSize',20)
% ylabel('Maximal mechanical energy \Delta E_{mec,max} [J]','FontSize',18)
% pmec = polyfit (log10(Delt),log10(Emec_max),1)
% hold on
% plot(Delt,10.^(pmec(2)).*Delt.^(pmec(1)),'r')
% grid on




%% diagramme de phase (question d-ii):
% 
% data1 = load(output{1});
% data2 = load(output{2});
% 
% 
% 
% 
% 
% % figure
% % plot(data(:,2), data(:,3) , 'r.')
% figure
% scatter (mod2pi(data1(:,2)),data1(:,3),'.')
% xlabel('\theta [rad]')
% ylabel('d\theta/dt [rad.s^{-1}]')
% 
% figure
% scatter (mod2pi(data2(:,2)),data2(:,3),'.')
% xlabel('\theta [rad]')
% ylabel('d\theta/dt [rad.s^{-1}]')




%% Lyapounov (question d-iv):



% data1 = load(output{1});
% data2 = load(output{2});
% data3 = data2 - data1;
% figure
% plot(data1(:,1),log(sqrt(data3(:,3).^2 + omega^2 .* (mod2pi(data3(:,2))).^2)))
% set(gca, 'fontsize', 24)
% xlabel('$t$', 'interpreter', 'latex')
% ylabel('$log(d(t))$', 'interpreter', 'latex')


%% question (d.v): sert a faire les diagrammes de phase 

% data = load(output{nsimul});
% 
% 
% 
% 
% figure
% scatter (mod2pi(data(:,2)),data(:,3),'.')
% xlabel('\theta [rad]')
% ylabel('d\theta/dt [rad.s^{-1}]')


%% Définition de la fonction de calcul de la norme : 






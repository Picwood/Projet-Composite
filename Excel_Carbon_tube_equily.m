syms teta Aa Bb Dd
%%%%%%% Description matériaux

%prop elast: Ex,Ey,Es(GPa) NUyx (coef de poisson)
E = [181 85 17 220; 10.3 5.6 3.5 220; 7.2 2.1 2.4 80].*10^9;
CPl = [0.28 0.34 0.28 0.3];
%Contrainte Ultime: X,X',Y,Y',S(MPa)
CU = [92.25 1410 540 1500;1500 280 210 1800;40 28 25 1500;246 141 105 1800;68 45 45 500];

%Ceof de dilatation: CTEx,CTEy(/MK)
CD = [1.57 -4 15.5 8;58.57 79 78.6 8];

%Température
T = [20,0];
dT = T(2)-T(1);

%Coeff de poisson trans
CPt = [0.0159 0.0224 0.0576 0.3];

nbmat = length(CPt);

alpha = 1./(1-CPl.*CPt); 
%Constante elastique
Qxx = E(1,:).*alpha; Qyy = E(2,:).*alpha; Qxy = CPl.*E(2,:).*alpha; Qss = E(3,:);
Q = zeros(3,3,nbmat); %GPa
Q(1,1,:) = Qxx; Q(1,2,:) = Qxy; Q(2,1,:) = Qxy; Q(2,2,:) = Qyy; Q(3,3,:) = Qss;

%Complaisance
Sxx = 1000./E(1,:); Syy = 1000./E(2,:); Sxy = -1000.*CPt./E(2,:); Sss = 1000./E(3,:);
S = zeros(3,3,nbmat); %1/TPa
S(1,1,:) = Sxx; S(1,2,:) = Sxy; S(2,1,:) = Sxy; S(2,2,:) = Syy; S(3,3,:) = Sss;

%Invariant elastique
U1 = (3.*Qxx + 3.*Qyy + 2*Qxy + 4.*Qss)./8;
U2 = (Qxx - Qyy)./2;
U3 = (Qxx + Qyy - 2.*Qxy - 4.*Qss)/8;
U4 = (Qxx + Qyy + 6.*Qxy - 4.*Qss)/8;

u1 = (3.*Sxx + 3.*Syy + 2.*Sxy + Sss)./8;
u2 = (Sxx - Syy)./2;
u3 = (Sxx + Syy - 2.*Sxy - Sss)./8;
u4 = (Sxx + Syy + 6.*Sxy - Sss)./8;

%Tsai Wu
Fxx = 1./(CU(1,:).*CU(2,:));
Fx = 1./CU(1,:) - 1./CU(2,:);
Fyy = 1./(CU(3,:).*CU(4,:));
Fy = 1./CU(3,:) - 1./CU(4,:);
Fss = 1./(CU(5,:).*CU(5,:));
Fxy = -0.5./(CU(1,:).*CU(2,:).*CU(3,:).*CU(4,:)).^0.5;

%%%%%%%%%% Stratifié non symétrique

%angleInt = linspace(10,80);
angleInt = [15 30 45 60 75];
angleExt = linspace(10,80);
deform = zeros(6,length(angleInt));
minus = zeros(6,6,length(angleInt));
d_11 = zeros(length(angleInt)^2,1);
d_66 = d_11;
ite = 1;
for j=1:length(angleInt)
    for k=1:length(angleExt)
        angl1 = angleExt(k);
        angl2 = angleInt(j);
        %strat: code(=1,2,3ou4 si mat et 0 si rien) epaisseur (mm) angle (°)
        e_ply=1;
        strat = [
            1 e_ply -angl1;
            1 e_ply angl2;
            ];
        
        strat(:,2) = strat(:,2).*0.001;
        nbstrat = length(strat(:,1));
        e_coeur = 10; %demi épaisseur du coeur vide
        htot = sum(strat(:,2));
        R = (e_coeur*0.001+sum(strat(find(strat(:,1)),2))./2);
        
        %position plan médian
        Z_index=0;
        for i=1:nbstrat
            if sum(strat(1:i,2))>=htot/2 && Z_index == 0
                Zpm_rel = sum(strat(1:i,2));
                Z_index = i;
                break
            end
        end
        coteZ = zeros(nbstrat,2);
        
        %calcul des côtes
        for i=1:nbstrat
            if i<=Z_index
                coteZ(i,1) = sum(strat(i:Z_index,2));
                coteZ(i,2) = sum(strat(i+1:Z_index,2));
            else
                coteZ(i,2) = -sum(strat(Z_index+1:i,2));
                coteZ(i,1) = -sum(strat(Z_index+1:i-1,2));
            end
        end
        
        %calcul de Q_ dans le repère local
        Q_ = zeros(3,3,nbstrat);
        for n=1:nbstrat
            i = strat(n,1);
            c = cos(strat(n,3)*pi/180);
            s = sin(strat(n,3)*pi/180);
            if i~=0
                Q_(1,1,n) = Q(1,1,i).*c^4+2.*(Q(1,2,i)+2*Q(3,3,i)).*s^2.*c^2+Q(2,2,i).*s^4;
                Q_(1,2,n) = Q(1,2,i).*(s^4+c^4)+(Q(1,1,i)+Q(2,2,i)-4.*Q(3,3,i)).*s^2.*c^2;
                Q_(2,1,n) = Q_(1,2,n);
                Q_(2,2,n) = Q(1,1,i).*s^4+2.*(Q(1,2,i)+2.*Q(3,3,i)).*s^2.*c^2+Q(2,2,i).*c^4;
                Q_(1,3,n) = (Q(1,1,i)-Q(1,2,i)-2.*Q(3,3,i)).*s.*c^3+(Q(1,2,i)-Q(2,2,i)+2.*Q(3,3,i)).*c.*s^3;
                Q_(3,1,n) = Q_(1,3,n);
                Q_(2,3,n) = (Q(1,1,i)-Q(1,2,i)-2.*Q(3,3,i)).*c.*s^3-(Q(1,2,i)-Q(2,2,i)+2.*Q(3,3,i)).*s.*c^3;
                Q_(3,2,n) = Q_(2,3,n);
                Q_(3,3,n) = (Q(1,1,i)+Q(2,2,i)-2.*Q(1,2,i)-2.*Q(3,3,i)).*s^2.*c^2+Q(3,3,i).*(s^4+c^4);
            end
        end
        
        A = sum(Q_.*reshape((coteZ(:,1)-coteZ(:,2)),1,1,nbstrat),3)/1; %GPa.m
        B = sum(Q_.*reshape((coteZ(:,1).^2-coteZ(:,2).^2),1,1,nbstrat),3)./2; %GPa.m²
        D = sum(Q_.*reshape((coteZ(:,1).^3-coteZ(:,2).^3),1,1,nbstrat),3)./3; %GPa.m^3
        E = sum(Q_.*reshape((coteZ(:,1).^4-coteZ(:,2).^4),1,1,nbstrat),3)./4; %GPa.m^3
        
        Aa = A;
        Bb = B+R.*cos(teta).*A;
        Dd = D + 2.*R.*cos(teta).*B+(R*cos(teta))^2.*A;
        Ee = E + 3.*R.*cos(teta).*D+3*(R*cos(teta))^2.*B+(R*cos(teta)).*A;
        
        A_ = eval(R.*int(Aa,teta,0,2*pi));
        B_ = eval(R.*int(Bb,teta,0,2*pi));
        D_ = eval(R.*int(Dd,teta,0,2*pi));
        E_ = eval(R.*int(Ee,teta,0,2*pi));
        
        MAJ = [A_ B_;B_ D_];
        minus(:,:,ite) = inv(MAJ);
        
        ite = ite + 1;
    end
end

figure
for i=1:length(angleInt)
hold on
plot(angleExt,squeeze(minus(4,4,(i-1)*length(angleExt)+1:(i)*length(angleExt))).*1000)
hold off
end
xlabel('angle $\theta_{ext}$ du stratifie extérieur','interpreter','latex')
ylabel('complaisance en flexion $\overline{d_{11}}$ $kN^{-1}.m^{-1}$','interpreter','latex')
legend('\theta_{int}='+string(angleInt)+'°')
grid on

figure
for i=1:length(angleInt)
hold on
plot(angleExt,squeeze(minus(6,6,(i-1)*length(angleExt)+1:(i)*length(angleExt))).*1000)
hold off
end
xlabel('angle $\theta_{ext}$ du stratifie extérieur','interpreter','latex')
ylabel('complaisance en torsion $\overline{d_{66}}$ $kN^{-1}.m^{-1}$','interpreter','latex')
legend('\theta_{int}='+string(angleInt)+'°')
grid on

%centre de masse
% Rin = linspace(0,11);
% vb = pi.*(15^2-Rin.^2)*1560;
% cm = (-198.31*53835+780.*vb)./(53835+vb);
% cm_mod = (-198.31*53835+780.*vb+1560*7500)./(53835+vb+7500);
% 
% figure
% plot(Rin,(780-cm)/10)
% hold on
% plot(Rin,(780-cm_mod)/10)
% xlabel('Rayon interne (mm)')
% ylabel('d_m-d_g(cm)')
% legend('Bâton nu','Bâton avec capuchon de 15g')
% grid on
% hold off



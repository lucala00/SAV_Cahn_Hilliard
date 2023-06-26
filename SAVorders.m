clear all; clf;
%%%%%%%%%

L1 = [1,0;0,1];
L2 = [0.5,-0.4;-0.4,0.5];


L = chol(L1'*L1-L2'*L2);
A = [0.25,0;0,2];
b = [0;0];

E = @(u)  sum((A*u -b).^2)/2;
nabla_E = @(u) A'*(A*u -b);
J2 = @(Mu)  (L2*Mu)'*(L2*Mu)/2;
J_relax1 =@(Mu,R) 0.5*(L1*Mu)'*(L1*Mu) - (R)^2;


T = 5;
U0=[0.1;2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% implicit/exact solver
Un_e = expm(-T*(L'*L)*(A'*A))*U0;
rel_err_u=@(U)  norm(U-Un_e,2);

dt_vec = [1e-4,5e-4,0.001,0.01,0.02];


for i=1:length(dt_vec)

    Un = U0; Mun = nabla_E(U0); Rn = sqrt(J2(Mun));

    Un_i = Un; Mun_i=Mun; Rn_i=Rn;
    Un_2 = Un; Mun_2 = Mun; Rn_2 = Rn;
    Un_2i = Un; Mun_2i = Mun; Rn_2i = Rn;
    Un_cvx = Un; Mun_cvx= Mun;

    dt=dt_vec(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAV order 1

    MA = [eye(2),dt*L1'*L1; -A'*A, eye(2)];

    for n=1:T/dt
        energie(n) = E(Un);
        pos1(n) = Un(1);
        pos2(n) = Un(2);
        t(n) = (n-1)*dt;

        %%%%%%%%%%%%% ordre 1 %%%%%%%%%%%
        temp = MA\[Un;-A'*b]; Un1 = temp(1:2); Mun1 = temp(3:4);
        temp = (MA\[dt*L2'*L2*Mun;[0;0]])/sqrt(J2(Mun));
        Un2 = temp(1:2); Mun2 = temp(3:4);
        hn = @(Mu) 1/(2*sqrt(J2(Mun)))*(L2*Mun)'*(L2*Mu);
        Rn = (Rn - hn(Mun) + hn(Mun1))/(1 - hn(Mun2));
        Rn = max(Rn,0);
        Un = Un1 + Rn*Un2;
        Mun = Mun1 + Rn*Mun2;

    end

    err_Un_SAV1(i)=rel_err_u(Un);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAV order 1 improved


    for n=1:T/dt
        energie_i(n) = E(Un_i);
        Rn_i = sqrt(J2(Mun_i));

        pos1_i(n) = Un_i(1);
        pos2_i(n) = Un_i(2);
        t_i(n) = (n-1)*dt;

        %%%%%%%%%%%%% ordre 1 %%%%%%%%%%%
        temp = MA\[Un_i;-A'*b]; Un1_i = temp(1:2); Mun1_i = temp(3:4);
        temp = (MA\[dt*L2'*L2*Mun_i;[0;0]])/sqrt(J2(Mun_i));
        Un2_i = temp(1:2); Mun2_i = temp(3:4);
        hn_1i = @(Mu) 1/(2*sqrt(J2(Mun_i)))*(L2*Mun_i)'*(L2*Mu);
        Rn_i = (Rn_i - hn_1i(Mun_i) + hn_1i(Mun1_i))/(1 - hn_1i(Mun2_i));
        Rn_i = max(Rn_i,0);
        Un_i = Un1_i + Rn_i*Un2_i;
        Mun_i = Mun1_i + Rn_i*Mun2_i;

    end

    err_Un_SAV1i(i)=rel_err_u(Un_i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAV order 2

    MA2 = [eye(2),dt*L1'*L1/2; -A'*A, eye(2)];

    for n=1:T/dt
        energie_2(n) = E(Un_2);
        ratio_2(n) = Rn_2/sqrt(J2(Mun_2));
        pos1_2(n) = Un_2(1); pos2_2(n) = Un_2(2); t_2(n) = (n-1)*dt;


        %%%%%%%%%%%%% ordre 1 // calcul de Mu^{n+1/2} %%%%%%%%%%%
        temp = MA\[Un_2;-A'*b]; Un1 = temp(1:2); Mun1 = temp(3:4);
        temp = (MA\[dt*L2'*L2*Mun_2;[0;0]])/sqrt(J2(Mun_2)); Un2 = temp(1:2); Mun2 = temp(3:4);
        hn_2 = @(Mu)1/(2*sqrt(J2(Mun_2)))*(L2*Mun_2)'*(L2*Mu);
        Rn_ = (Rn_2 - hn_2(Mun_2) + hn_2(Mun1))/(1 - hn_2(Mun2));  Rn_ = max(Rn_,0);
        Mun_1_2 = ((Mun1 + Rn_*Mun2) + Mun_2)/2;


        %%%%%%%%%%% ordre 2 // calcul de Mu^{n+1} %%%%%%%%%%%%
        temp = MA2\[Un_2 - dt*(L1'*L1)/2*Mun_2;-A'*b]; Un1 = temp(1:2); Mun1 = temp(3:4);
        temp = (MA2\[dt*L2'*L2*Mun_1_2;[0;0]])/sqrt(J2(Mun_1_2)); Un2 = temp(1:2); Mun2 = temp(3:4);
        hn_21 = @(Mu)1/(2*sqrt(J2(Mun_1_2)))*(L2*Mun_1_2)'*(L2*Mu);
        Rn_p =  max((Rn_2 - hn_21(Mun_2) + hn_21(Mun1)  + 1/2*hn_21(Mun2)*Rn_2)/(1 - hn_21(Mun2)/2),0);
        Un_2 = Un1 + (Rn_2 + Rn_p)/2*Un2; Mun_2 = Mun1 + (Rn_2 + Rn_p)/2*Mun2;  Rn_2= Rn_p;

        J_rel2(n) = J_relax1(Mun_2,Rn_2);

    end

    err_Un_SAV2(i)=rel_err_u(Un_2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAV order 2 (improved)

    MA2 = [eye(2),dt*L1'*L1/2; -A'*A, eye(2)];

    for n=1:T/dt
        energie_2i(n) = E(Un_2i);
        ratio_2i(n) = Rn_2i/sqrt(J2(Mun_2i));
        Rn_2i = sqrt(J2(Mun_2i));
        pos1_2i(n) = Un_2i(1); pos2_2i(n) = Un_2i(2); t_2i(n) = (n-1)*dt;


        %%%%%%%%%%%%% ordre 1 // calcul de Mu^{n+1/2} %%%%%%%%%%%
        temp = MA\[Un_2i;-A'*b]; Un1_2i = temp(1:2); Mun1_2i = temp(3:4);
        temp = (MA\[dt*L2'*L2*Mun_2i;[0;0]])/sqrt(J2(Mun_2i)); Un2_2i = temp(1:2); Mun2_2i = temp(3:4);
        hn_2i = @(Mu)1/(2*sqrt(J2(Mun_2i)))*(L2*Mun_2i)'*(L2*Mu);
        Rn_ = (Rn_2i - hn_2i(Mun_2i) + hn_2i(Mun1_2i))/(1 - hn_2i(Mun2_2i));  Rn_ = max(Rn_,0);
        Mun_1_2i = ((Mun1_2i + Rn_*Mun2_2i) + Mun_2i)/2;


        %%%%%%%%%%% ordre 2 // calcul de Mu^{n+1} %%%%%%%%%%%%
        temp = MA2\[Un_2i - dt*(L1'*L1)/2*Mun_2i;-A'*b]; Un1_2i = temp(1:2); Mun1_2i = temp(3:4);
        temp = (MA2\[dt*L2'*L2*Mun_1_2i;[0;0]])/sqrt(J2(Mun_1_2i)); Un2_2i = temp(1:2); Mun2_2i = temp(3:4);
        hn_2i1 = @(Mu)1/(2*sqrt(J2(Mun_1_2i)))*(L2*Mun_1_2i)'*(L2*Mu);
        Rn_p_2i =  max((Rn_2i - hn_2i1(Mun_2i) + hn_2i1(Mun1_2i)  + 1/2*hn_2i1(Mun2_2i)*Rn_2i)/(1 - hn_2i1(Mun2_2i)/2),0);
        Un_2i = Un1_2i + (Rn_2i + Rn_p_2i)/2*Un2_2i; Mun_2i = Mun1_2i + (Rn_2i + Rn_p_2i)/2*Mun2_2i;  Rn_2i= Rn_p_2i;

        J_rel2i(n) = J_relax1(Mun_2i,Rn_2i);

    end

    err_Un_SAV2i(i)=rel_err_u(Un_2i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Convex splitting

    MA = [eye(2),dt*L1'*L1; -A'*A, eye(2)];

    for n=1:T/dt
        energie_cvx(n) = E(Un_cvx);
        pos1_cvx(n) = Un_cvx(1);
        pos2_cvx(n) = Un_cvx(2);
        J_cvx(n) = 0.5*(L1*Mun_cvx)'*(L1*Mun_cvx) - (L2*Mun_cvx)'*(L2*Mun_cvx)/2;
        t_cvx(n) = (n-1)*dt;
        temp = MA\[Un_cvx + dt*(L2'*L2)*Mun_cvx;-A'*b];
        Un_cvx = temp(1:2); Mun_cvx = temp(3:4);


    end

    err_Un_cvx(i)=rel_err_u(Un_cvx);


end

figure,
loglog(dt_vec,err_Un_SAV1,'b',dt_vec,err_Un_SAV1i,'b--', dt_vec, err_Un_SAV2,'g', dt_vec,err_Un_SAV2i,'g--',...
    dt_vec,err_Un_cvx,'c','LineWidth', 2.5)
hold on
loglog(dt_vec,1e-7*dt_vec,'r',dt_vec,1e-7*dt_vec.^2,'r--','LineWidth', 1)
legend('Mb-SAV 1', 'Mb-SAV 1 improved', 'Mb-SAV 2','Mb-SAV 2 improved','Cvx split','$O(\delta t)$','$O(\delta t^2)$','Location','Best', 'FontSize',20,...
    'interpreter','latex')
xlabel('$\delta t$','FontSize',20,'interpreter','latex'), ylabel('$\|u^N(\delta t)-u_{ex}\|_2$','FontSize',20,'interpreter','latex')
title('Order of methods','FontSize',24,'interpreter','latex')

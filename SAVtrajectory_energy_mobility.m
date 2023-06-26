clear all; clf;
close all
%%%%%%%%%

L1 = [1,0;0,1];
L2 = [0.5,-0.4;-0.4,0.5];

L = chol(L1'*L1-L2'*L2);
A = [0.25,0;0,2];
b = [0;0];

E = @(u)  sum((A*u -b).^2)/2;
nabla_E = @(u) A'*(A*u -b);

J_total = @(Mu)  (L*Mu)'*(L*Mu)/2;

J2 = @(Mu)  (L2*Mu)'*(L2*Mu)/2;
J_relax1 =@(Mu,R) 0.5*(L1*Mu)'*(L1*Mu) - (R)^2;


T = 200;
U0=[0.1;2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% implicit/exact solver
exact_solve = @(t,U) expm(-t*(L'*L)*(A'*A))*U;

Un_e_fin = exact_solve(T,U0);
pos1_e_fin=Un_e_fin(1); pos2_e_fin=Un_e_fin(2);

dt_vec = [0.1, 1, 4];


for i=1:length(dt_vec)

    Un = U0; Mun = nabla_E(U0); Rn = sqrt(J2(Mun));


    Un_e = Un; Mun_e = Mun; Rn_e=Rn;
    Un_i = Un; Mun_i=Mun; Rn_i=Rn;
    Un_2 = Un; Mun_2 = Mun; Rn_2 = Rn;
    Un_2i = Un; Mun_2i = Mun; Rn_2i = Rn;
    Un_cvx = Un; Mun_cvx= Mun;

    dt=dt_vec(i);

    for n=1:T/dt
        energie_e(n,i)=E(Un_e);
        Mun_e = nabla_E(Un_e);
        J_e(n,i)=J_total(Mun_e);
        J2e_fun(n,i)=J2(Mun_e);

        t_ex(n,i)=(n-1)*dt;
        Un_e= exact_solve(dt,Un_e);
        pos1_e(n,i) = Un_e(1); pos2_e(n,i) =Un_e(2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAV order 1

    MA = [eye(2),dt*L1'*L1; -A'*A, eye(2)];

    for n=1:T/dt
        energie(n,i) = E(Un);
        J(n,i)=J_total(Mun);
        J_relax(n,i)=J_relax1(Mun,Rn);
        pos1(n,i) = Un(1);
        pos2(n,i) = Un(2);
        rn(n,i)=Rn;
        t(n,i) = (n-1)*dt;

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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAV order 1 improved


    for n=1:T/dt
        energie_i(n,i) = E(Un_i);
        J_i(n,i)=J_total(Mun_i);
        J_relax_i(n,i)=J_relax1(Mun_i,Rn_i);
        rn_i(n,i)=Rn_i;
        Rn_i = sqrt(J2(Mun_i));
        pos1_i(n,i) = Un_i(1);
        pos2_i(n,i) = Un_i(2);
        t_i(n,i) = (n-1)*dt;

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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAV order 2

    MA2 = [eye(2),dt*L1'*L1/2; -A'*A, eye(2)];

    for n=1:T/dt
        energie_2(n,i) = E(Un_2);
        J_2(n,i)=J_total(Mun_2);
        J_relax_2(n,i)=J_relax1(Mun_2,Rn_2);
        pos1_2(n,i) = Un_2(1); pos2_2(n,i) = Un_2(2); t_2(n,i) = (n-1)*dt;
        rn2(n,i)=Rn_2;

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

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAV order 2 (improved)

    MA2 = [eye(2),dt*L1'*L1/2; -A'*A, eye(2)];

    for n=1:T/dt
        energie_2i(n,i) = E(Un_2i);
        J_2i(n,i)=J_total(Mun_2i);
        J_relax_2i(n,i)=J_relax1(Mun_2i,Rn_2i);
        ratio_2i(n,i) = Rn_2i/sqrt(J2(Mun_2i));
        rn2i(n,i)=Rn_2i;
        Rn_2i = sqrt(J2(Mun_2i));
        pos1_2i(n,i) = Un_2i(1); pos2_2i(n,i) = Un_2i(2); t_2i(n,i) = (n-1)*dt;


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

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Convex splitting

    MA = [eye(2),dt*L1'*L1; -A'*A, eye(2)];

    for n=1:T/dt
        energie_cvx(n,i) = E(Un_cvx);
        Mun_cvx = nabla_E(Un_cvx);
        J_cvx(n,i)=J_total(Mun_cvx);
        pos1_cvx(n,i) = Un_cvx(1);
        pos2_cvx(n,i) = Un_cvx(2);
        t_cvx(n,i) = (n-1)*dt;
        temp = MA\[Un_cvx + dt*(L2'*L2)*Mun_cvx;-A'*b];
        Un_cvx = temp(1:2); Mun_cvx = temp(3:4);

    end


    xmin = min(min([pos1(:,i),pos1_i(:,i),pos1_2(:,i),pos1_2i(:,i),pos1_cvx(:,i)]));
    xmax = max(max([pos1(:,i),pos1_i(:,i),pos1_2(:,i),pos1_2i(:,i),pos1_cvx(:,i)]));


    ymin = min(min([pos2(:,i),pos2_i(:,i),pos2_2(:,i),pos2_2i(:,i),pos2_cvx(:,i)]));
    ymax = max(max([pos2(:,i),pos2_i(:,i),pos2_2(:,i),pos2_2i(:,i),pos2_cvx(:,i)]));

    x1=linspace(xmin-0.2,xmax+0.2);
    x2 = linspace(ymin-0.2,ymax+0.2);
    [X1,X2] = meshgrid(x1,x2);
    E_X = (A(1,1)*X1.^2 + A(2,2)*X2.^2)/2;


    figure(i),
    plot(pos1(:,i),pos2(:,i),'b', 'LineWidth', 2)
    hold on;
    plot(pos1_i(:,i),pos2_i(:,i),'b--', 'LineWidth', 2)
    plot(pos1_2(:,i),pos2_2(:,i),'g', 'LineWidth', 2)
    plot(pos1_2i(:,i),pos2_2i(:,i),'g--', 'LineWidth', 2)
    plot(pos1_cvx(:,i),pos2_cvx(:,i),'c', 'LineWidth', 2)
    plot(pos1_e(:,i),pos2_e(:,i),'r', 'LineWidth', 1.5)
    contour(x1,x2,E_X,'--')
    legend('Mb-SAV 1', 'Mb-SAV 1 improved', 'Mb-SAV 2','Mb-AV 2 improved','Cvx split','Exact','Location','Best', ...
        'FontSize', 14, 'interpreter', 'latex')
    axis equal
    title(['Trajectory of $u^n$, $\delta_t$ =', num2str(dt)], 'FontSize', 25,'interpreter','latex')
    xlabel('$x$','FontSize',20,'interpreter','latex'), ylabel('$y$','FontSize',20,'interpreter','latex')

    figure(i+length(dt_vec))
    loglog(t(:,i),energie(:,i),'b','LineWidth',2.5);
    hold on,
    loglog(t_i(:,i),energie_i(:,i),'b--','LineWidth',2.5);
    loglog(t_2(:,i),energie_2(:,i),'g','LineWidth',2.5);
    loglog(t_2(:,i),energie_2i(:,i),'g--','LineWidth',2.5);
    loglog(t_cvx(:,i),energie_cvx(:,i),'c','LineWidth',2.5)
    loglog(t_ex(:,i),energie_e(:,i), 'r', 'LineWidth', 2.5)
    loglog(t_cvx(:,i),energie_cvx(2,i)*ones(size(t_cvx(:,i))),'k:','LineWidth',1)
    legend('Mb-SAV 1', 'Mb-SAV 1 improved', 'Mb-SAV 2','Mb-AV 2 improved','Cvx split','Exact','Location','Best','FontSize', 14, 'interpreter', 'latex')
    title(['Energy $E$ decay, $\delta_t$ =', num2str(dt)], 'FontSize', 25, 'interpreter', 'latex')
    xlabel('Iterations','FontSize',20,'interpreter','latex'), ylabel('$E(u^n)$','FontSize',20,'interpreter','latex')

    figure(i+2*length(dt_vec))
    loglog(t(:,i),J(:,i),'b','LineWidth',2.5);
    hold on,
    loglog(t_i(:,i),J_i(:,i),'b--','LineWidth',2.5);
    loglog(t_2(:,i),J_2(:,i),'g','LineWidth',2.5);
    loglog(t_2(:,i),J_2i(:,i),'g--','LineWidth',2.5);
    loglog(t_cvx(:,i),J_cvx(:,i),'c','LineWidth',2.5)
    loglog(t_ex(:,i),J_e(:,i), 'r', 'LineWidth', 2.5)
    legend('Mb-SAV 1', 'Mb-SAV 1 improved', 'Mb-SAV 2','Mb-AV 2 improved','Cvx split','Exact','Location','Best', 'FontSize', 14, 'interpreter', 'latex')
    title(['Mobility $J$ decay, $\delta_t$ =', num2str(dt)], 'FontSize', 25, 'interpreter', 'latex')
    xlabel('Iterations','FontSize',20,'interpreter','latex'), ylabel('$J(u^n)$','FontSize',20,'interpreter','latex')

    figure(i+3*length(dt_vec))
    loglog(t_i(:,i),J_relax(:,i),'b','LineWidth',2.5);
    hold on,
    loglog(t_2(:,i),J_relax_i(:,i),'b--','LineWidth',2.5);
    loglog(t_2(:,i),J_relax_2(:,i),'g','LineWidth',2.5);
    loglog(t_cvx(:,i),J_relax_2i(:,i),'g--','LineWidth',2.5)
    loglog(t(:,i),J_e(:,i),'r','LineWidth',2.5);
    legend('Mb-SAV 1', 'Mb-SAV 1 improved', 'Mb-SAV 2','Mb-AV 2 improved','Exact',...
        'Location','Best', 'FontSize', 14, 'interpreter', 'latex')
    title(['Relaxed mobility $J$ decay, $\delta_t$ =', num2str(dt)], 'FontSize', 25, 'interpreter', 'latex')
    xlabel('Iterations','FontSize',20,'interpreter','latex'), ylabel('$\tilde{J}(u^n)$','FontSize',20,'interpreter','latex')

    figure(i+4*length(dt_vec))
    semilogx(t_i(:,i),rn(:,i),'b','LineWidth',2.5);
    hold on
    semilogx(t_2(:,i),rn_i(:,i),'b--','LineWidth',2.5);
    semilogx(t_2(:,i),rn2(:,i),'g','LineWidth',2.5);
    semilogx(t_2(:,i),rn2i(:,i),'g--','LineWidth',2.5)
    semilogx(t(:,i),sqrt(J2e_fun(:,i)),'r','LineWidth',2.5);
    legend('Mb-SAV 1', 'Mb-SAV 1 improved', 'Mb-SAV 2','Mb-AV 2 improved','Exact',...
        'Location','Best', 'FontSize', 14, 'interpreter', 'latex')
    title(['$r^n$ approximation of $\sqrt{J_2(\mu^n)}$, $\delta_t$ =', num2str(dt)], 'FontSize', 25, 'interpreter', 'latex')
    xlabel('Iterations','FontSize',20,'interpreter','latex'), ylabel('$r^n$','FontSize',20,'interpreter','latex')

    hold off

end

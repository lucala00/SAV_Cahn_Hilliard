clear all; clf; close all;   colormap('jet');   
N = 2^8; x1 = linspace(-1/2,1/2,N); x2 = linspace(-1/2,1/2,N); [X1,X2] = meshgrid(x1,x2);
%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%
epsilon =2/N; dt = epsilon^4; T =1*10^(-4);
alpha = 2; W_prim = @(U) (U.*(U-1).*(2*U-1));W = @(U) 1/2*(U.*(U-1)).^2;
k =  [0:N/2,-N/2+1:-1]; [K1,K2] = meshgrid(k,k);
Delta = -4*pi^2*((K1.^2 + (K2).^2));
T_vec = linspace(0,0.99*T,10);
T_vec = [T_vec,2*T];

%%%%%%%%%%%%%%%%%% condition initiale %%%%%%%%%%%
R_0 = 0.25;
rho =  sqrt(X1.^2 + X2.^2);
theta = angle(X1 + 1i*X2);
d1 = rho - R_0*( 1 + 0.5*cos(4*theta)); 
U_init = (1-tanh(d1/(2*epsilon)))/2; Uinit_fourier = fft2(U_init);
Muinit_fourier = -Delta.*Uinit_fourier + fft2(1/epsilon^2*W_prim(U_init)); Mu_init = real(ifft2(Muinit_fourier));


%%%%%%%%%%%%%%%%%%%%% modèle M %%%%%%%%%%
U = U_init; U_fourier = Uinit_fourier;Mu_fourier = Muinit_fourier; Mu= Mu_init;
nabla1_Mu= ifft2(2*1i*pi*K1.*Mu_fourier);nabla2_Mu = ifft2(2*1i*pi*K2.*Mu_fourier);
MobM = @(U) 36*((((U).*(1-U)).^2));
s = linspace(0,1,N); m=max(MobM(s));
M_L = 1./(1 + dt*m*Delta.*(Delta - alpha/epsilon^2));
j_sauvegarde = 1;

for i=1:T/dt,
temps(i) = i*dt;    
mobU = MobM(U);
div_mob_laplacien_fourier =  2*1i*pi*K1.*fft2((mobU-m).*nabla1_Mu ) + 2*1i*pi*K2.*fft2((mobU-m).*nabla2_Mu);
J2mu = real(sum( (m-mobU(:)).*(nabla1_Mu(:).^2 + nabla2_Mu(:).^2) + eps^2)/N^2)/2;
R = sqrt(J2mu);
%%%%%%%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1_1 = U_fourier; B1_2 = dt*div_mob_laplacien_fourier/sqrt(J2mu);
B2 = fft2(W_prim(U)/epsilon^2 - alpha/epsilon^2*U);
Mu1_fourier = M_L.*((alpha/epsilon^2 - Delta).*B1_1 + B2);
Mu2_fourier = M_L.*((alpha/epsilon^2 - Delta).*B1_2); 
nabla1_Mu1 = real(ifft2(2*1i*pi*K1.*Mu1_fourier));nabla2_Mu1 = real(ifft2(2*1i*pi*K2.*Mu1_fourier));
nabla1_Mu2 = real(ifft2(2*1i*pi*K1.*Mu2_fourier));nabla2_Mu2 = real(ifft2(2*1i*pi*K2.*Mu2_fourier));

%%%%%%%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_mu1_mu = real((sum( (m-mobU(:)).*(nabla1_Mu1(:).*nabla1_Mu(:) + nabla2_Mu1(:).*nabla2_Mu(:))+eps^2)/N^2)/(2*sqrt(J2mu)));
h_mu2_mu = real((sum( (m-mobU(:)).*(nabla1_Mu2(:).*nabla1_Mu(:) + nabla2_Mu2(:).*nabla2_Mu(:))+eps^2)/N^2)/(2*sqrt(J2mu)));
R = (h_mu1_mu)/(1-  h_mu2_mu);
%%%%%%%% Final step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mu_fourier = Mu1_fourier + R*Mu2_fourier; nabla1_Mu = nabla1_Mu1 + R*nabla1_Mu2;  nabla2_Mu = nabla2_Mu1 + R*nabla2_Mu2;
U_fourier = M_L.*(B1_1 + R*B1_2 + dt*m*Delta.*B2); U = real(ifft2(U_fourier));

EnergieM(i) = epsilon*sum(sum( (4*pi^2*(K1.^2 + K2.^2)).*abs(U_fourier).^2 ))/N^4 + 1/epsilon*sum(W(U(:)))/N^2    ;
 
if (i*dt > T_vec(j_sauvegarde))
       
           clf;
        colormap('jet')

       imagesc(x1,x2,U)    
       colorbar
       name_title = ['t = ',num2str(i*dt)];
       title(name_title,'linewidth',2)
       axis tight
       
       name_fig = ['Test1_M_',num2str( j_sauvegarde),'.eps'];
       print('-depsc', name_fig)
       j_sauvegarde = j_sauvegarde+1;

 

   end


end
U_M = U;

%%%%%%%%%%%%%%%%%%%%% modèle NMM %%%%%%%%%%
U = U_init; U_fourier = Uinit_fourier;Mu_fourier = Muinit_fourier; Mu= Mu_init;
nabla1_mu= ifft2(2*1i*pi*K1.*Mu_fourier);nabla2_mu = ifft2(2*1i*pi*K2.*Mu_fourier);
MobM = @(U) ((((U).*(1-U)).^2+epsilon^2)); MobN = @(U) 1./sqrt(MobM(U) );
%%%%%%%%%%%%%%%%% Diffusion Fourier %%%%%%%%%%%%%%%
m=1; beta = 2/epsilon^2;
M_L = 1./(1 + dt*(m*Delta - beta).*(Delta - alpha/epsilon^2));
j_sauvegarde = 1;


for i=1:T/dt,
mobMU = MobM(U); mobNU = MobN(U); sqrtM = sqrt(mobMU); sqrtM_fourier = fft2(sqrtM); 
nabla1_sqrtM=  real(ifft2(2*pi*1i*K1.*sqrtM_fourier )); nabla2_sqrtM=  real(ifft2(2*pi*1i*K2.*sqrtM_fourier ));
muN_fourier = fft2(Mu.*mobNU); muN = real(ifft2(muN_fourier));
nabla1_muN =  real(ifft2(2*pi*1i*K1.*muN_fourier)); nabla2_muN =  real(ifft2(2*pi*1i*K2.*muN_fourier));
laplacien_muN =  real(ifft2(Delta.*muN_fourier)); NdivMgradNMu = sqrtM.*laplacien_muN + 2*(nabla1_sqrtM.*nabla1_muN +nabla2_sqrtM.*nabla2_muN);
J2mu = (sum( -mobMU(:).*(nabla1_muN(:).^2 + nabla2_muN(:).^2 ) +  m*((nabla1_mu(:).^2 + nabla2_mu(:).^2 )) + beta*Mu(:).^2 + eps^2)/N^2)/2;
R = sqrt(J2mu);
%%%%%%%%%%%%%%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1_1 = U_fourier; B1_2 = dt*(fft2(NdivMgradNMu) - (m*Delta-beta).*Mu_fourier)/sqrt(J2mu);
B2 = fft2(W_prim(U)/epsilon^2 - alpha/epsilon^2*U);
Mu1_fourier = M_L.*((alpha/epsilon^2 - Delta).*B1_1 + B2);   Mu1 = real(ifft2(Mu1_fourier));
Mu2_fourier = M_L.*((alpha/epsilon^2 - Delta).*B1_2);   Mu2 = real(ifft2(Mu2_fourier));
%%%%%%%%%%%%%%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu1N_fourier = fft2(Mu1.*mobNU); nabla1_mu1N =  real(ifft2(2*pi*1i*K1.*mu1N_fourier)); nabla2_mu1N =  real(ifft2(2*pi*1i*K2.*mu1N_fourier));
mu2N_fourier = fft2(Mu2.*mobNU); nabla1_mu2N =  real(ifft2(2*pi*1i*K1.*mu2N_fourier)); nabla2_mu2N =  real(ifft2(2*pi*1i*K2.*mu2N_fourier));
nabla1_mu1 =  real(ifft2(2*pi*1i*K1.*Mu1_fourier)); nabla2_mu1 =  real(ifft2(2*pi*1i*K2.*Mu1_fourier));
nabla1_mu2 =  real(ifft2(2*pi*1i*K1.*Mu2_fourier)); nabla2_mu2 =  real(ifft2(2*pi*1i*K2.*Mu2_fourier));
h_mu1_mu = (sum( -mobMU(:).*(nabla1_muN(:).*nabla1_mu1N(:) + nabla2_muN(:).*nabla2_mu1N(:) ) ...
    +  m*((nabla1_mu(:).*nabla1_mu1(:) + nabla2_mu(:).*nabla2_mu1(:) )) + beta*Mu(:).*Mu1(:) + eps^2)/N^2)/(2*sqrt(J2mu));
h_mu2_mu = (sum( -mobMU(:).*(nabla1_muN(:).*nabla1_mu2N(:) + nabla2_muN(:).*nabla2_mu2N(:) ) ...
    +  m*((nabla1_mu(:).*nabla1_mu2(:) + nabla2_mu(:).*nabla2_mu2(:) )) + beta*Mu(:).*Mu2(:) + eps^2)/N^2)/(2*sqrt(J2mu));
R = h_mu1_mu/(1-  h_mu2_mu);
%%%%%%%%%%%%%%%% Final Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Mu_fourier = Mu1_fourier + R*Mu2_fourier; nabla1_mu = nabla1_mu1 + R*nabla1_mu2;  nabla2_mu = nabla2_mu1 + R*nabla2_mu2; Mu = Mu1 + R*Mu2;
U_fourier = M_L.*(B1_1 + R*B1_2 + dt*(m*Delta-beta).*B2); U = real(ifft2(U_fourier));


EnergieNMN(i) = epsilon*sum(sum( (4*pi^2*(K1.^2 + K2.^2)).*abs(U_fourier).^2 ))/N^4 + 1/epsilon*sum(W(U(:)))/N^2    ;
 
if (i*dt > T_vec(j_sauvegarde))
       
           clf;
        colormap('jet')

       imagesc(x1,x2,U)    
       colorbar
       name_title = ['t = ',num2str(i*dt)];
       title(name_title,'linewidth',2)
       axis tight
       
       name_fig = ['Test1_NMN_',num2str( j_sauvegarde),'.eps'];
       print('-depsc', name_fig)
       j_sauvegarde = j_sauvegarde+1;
 

   end


end
U_NMN = U;








clf
plot(temps,EnergieM,'r');
hold on;
plot(temps,EnergieNMN);
name_title = 'Cahn Hilliard energy';     
title(name_title,'linewidth',2)
axis tight
name_fig = ['Test1_energy.eps'];
legend('MCH model','NMNCH model')
print('-depsc', name_fig)
    

%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
plot(x1,U_M(:,N/2),'r');
hold on;
plot(x1,U_NMN(:,N/2));
name_title = 'Profile of the phase field function';     
title(name_title,'linewidth',2)
axis([-0.5,0.5,0-1.5*epsilon,1+epsilon])
legend('MCH model','NMNCH model')
name_fig = ['Test1_profile.eps'];
print('-depsc', name_fig)


clf
plot(x1,U_M(:,N/2),'r');
hold on;
plot(x1,U_NMN(:,N/2));
name_title = 'Profile of the phase field function';     
title(name_title,'linewidth',2)
axis([-0.5,0.5,1-1.5*epsilon,1+epsilon])
legend('MCH model','NMNCH model')
name_fig = ['Test1_zoom_profile.eps'];
print('-depsc', name_fig)



clf
plot(x1,U_M(:,N/2),'r');
hold on;
plot(x1,U_NMN(:,N/2));
name_title = 'Profile of the phase field function';     
title(name_title,'linewidth',2)
axis([0.21-4*epsilon,0.21+4*epsilon,1-1.5*epsilon,1+epsilon])
legend('MCH model','NMNCH model')
name_fig = ['Test1_zoom2_profile.eps'];
print('-depsc', name_fig)






clear all;
 
N1 = 2^8; N2 = N1/4; N3 = N1/4;
L1 = 1; L2 = 1/4; L3 = 1/4;

x1 = linspace(-L1/2,L1/2,N1);
x2 = linspace(-L2/2,L2/2,N2);
x3 = linspace(-L3/2,L3/2,N3);

[X1,X2,X3] = meshgrid(x1,x2,x3);

%%%%%%%%%%%%%%%%% parametre %%%%%%%%%%%%%%%%%%%%%%%%
epsilon =2/N1;
dt = epsilon^4;
T =5*10^(-6);

alpha = 2;
W = @(U) 1/2*(U.*(U-1)).^2;
W_prim = @(U) (U.*(U-1).*(2*U-1));
MobM = @(U) 36*((((U).*(1-U)).^2));

k1 =  [0:N1/2,-N1/2+1:-1]/L1; k2 =  [0:N2/2,-N2/2+1:-1]/L2; k3 = [0:N3/2,-N3/2+1:-1]/L3;
[K1,K2,K3] = meshgrid(k1,k2,k3);
Delta = -4*pi^2*((K1.^2 + K2.^2 + K3.^2));



%%%%%%%%%%%%%%%%%% condition initiale %%%%%%%%%%%

d = max(max(max(max(max(X3-0.2*L3/2,-X3-0.2*L3/2),X1-0.9*L1/2),X2-0.2*L2/2),-X1-0.9*L1/2),-X2-0.2*L2/2); 
U = (1-tanh(d/(epsilon)))/2; U_fourier = fftn(U);
Mu_fourier = -Delta.*U_fourier + fftn(1/epsilon^2*W_prim(U)); Mu = real(ifftn(Mu_fourier));
nabla1_Mu= ifftn(2*1i*pi*K1.*Mu_fourier);nabla2_Mu = ifftn(2*1i*pi*K2.*Mu_fourier);nabla3_Mu = ifftn(2*1i*pi*K3.*Mu_fourier);


s = linspace(0,1,100); m=max(MobM(s)); 
M_L = 1./(1 + dt*m*Delta.*(Delta - alpha/epsilon^2));

%%%%%%%%%%%%%%%%
k = 1;

T_vec = linspace(0,0.99*T,10);
T_vec = [T_vec(1),2*10^(-7),T_vec(2:10)];
T_vec = [T_vec,2*T];
j_sauvegarde = 1;


 clf
  v = U;
 p = patch(isosurface(x1,x2,x3,v,0.5));
 isonormals(x1,x2,x3,v,p)
 set(p,'FaceColor',[1,0,0],'EdgeColor','none');
 
 daspect([1 1 1]);
 camlight headlight;
 lighting gouraud;
 axis([-L1/2,L1/2,-L2/2,L2/2,-L3/2,L3/2])

 
for i=1:T/dt,

mobU = MobM(U);
div_mob_laplacien_fourier =  2*1i*pi*K1.*fftn((mobU-m).*nabla1_Mu ) ...
    + 2*1i*pi*K2.*fftn((mobU-m).*nabla2_Mu) ... 
    + 2*1i*pi*K3.*fftn((mobU-m).*nabla3_Mu);
J2mu = real(sum( (m-mobU(:)).*(nabla1_Mu(:).^2 + nabla2_Mu(:).^2 + nabla3_Mu(:).^2 ) + eps^2)/(N1*N2*N3))/2;
R = sqrt(J2mu);

%%%%%%%% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1_1 = U_fourier; B1_2 = dt*div_mob_laplacien_fourier/sqrt(J2mu);
B2 = fftn(W_prim(U)/epsilon^2 - alpha/epsilon^2*U);
Mu1_fourier = M_L.*((alpha/epsilon^2 - Delta).*B1_1 + B2);
Mu2_fourier = M_L.*((alpha/epsilon^2 - Delta).*B1_2); 
nabla1_Mu1 = real(ifftn(2*1i*pi*K1.*Mu1_fourier));nabla2_Mu1 = real(ifftn(2*1i*pi*K2.*Mu1_fourier));nabla3_Mu1 = real(ifftn(2*1i*pi*K3.*Mu1_fourier));
nabla1_Mu2 = real(ifftn(2*1i*pi*K1.*Mu2_fourier));nabla2_Mu2 = real(ifftn(2*1i*pi*K2.*Mu2_fourier));nabla3_Mu2 = real(ifftn(2*1i*pi*K3.*Mu2_fourier));
%%%%%%%% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_mu1_mu = real((sum( (m-mobU(:)).*(nabla1_Mu1(:).*nabla1_Mu(:) + nabla2_Mu1(:).*nabla2_Mu(:) + nabla3_Mu1(:).*nabla3_Mu(:)  )+eps^2)/(N1*N2*N3))/(2*sqrt(J2mu)));
h_mu2_mu = real((sum( (m-mobU(:)).*(nabla1_Mu2(:).*nabla1_Mu(:) + nabla2_Mu2(:).*nabla2_Mu(:) + nabla3_Mu2(:).*nabla3_Mu(:)  )+eps^2)/(N1*N2*N3))/(2*sqrt(J2mu)));
R = (h_mu1_mu)/(1-  h_mu2_mu);    
%%%%%%%% Final step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mu_fourier = Mu1_fourier + R*Mu2_fourier; nabla1_Mu = nabla1_Mu1 + R*nabla1_Mu2;  nabla2_Mu = nabla2_Mu1 + R*nabla2_Mu2; nabla3_Mu = nabla3_Mu1 + R*nabla3_Mu2;
U_fourier = M_L.*(B1_1 + R*B1_2 + dt*m*Delta.*B2); U = real(ifftn(U_fourier));


EnergieNMN(i) = epsilon*sum(sum(sum( (4*pi^2*(K1.^2 + K2.^2 + K3.^2)).*abs(U_fourier).^2 )))/(N1*N2*N3)^2 + 1/epsilon*sum(W(U(:)))/(N1*N2*N3)    ;
 
    


 if (i*dt > T_vec(j_sauvegarde))
              
       
        clf
       v = U;
       p = patch(isosurface(x1,x2,x3,v,0.5));
       isonormals(x1,x2,x3,v,p)
       set(p,'FaceColor',[1,0,0],'EdgeColor','none');
 
       daspect([1 1 1]);
       camlight headlight;
       lighting gouraud;
      axis([-L1/2,L1/2,-L2/2,L2/2,-L3/2,L3/2])

       name_title = ['t = ',num2str(i*dt)];
       title(name_title,'linewidth',2)
      
       
       name_fig = ['Test3D1_UM_',num2str( j_sauvegarde),'.eps'];
       print('-depsc', name_fig)
%      
       
       j_sauvegarde = j_sauvegarde +1;
 

       
    end


% 
% 
% if (mod(i,10)==1)
%        clf
%        v = U;
%        p = patch(isosurface(x1,x2,x3,v,0.5));
%        isonormals(x1,x2,x3,v,p)
%        set(p,'FaceColor',[1,0,0],'EdgeColor','none');
%  
%        daspect([1 1 1]);
%        camlight headlight;
%        lighting gouraud;
%       axis([-L1/2,L1/2,-L2/2,L2/2,-L3/2,L3/2])
% 
%     
%     
%     drawnow
%    
% end
    
end






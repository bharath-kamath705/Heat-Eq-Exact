%=======================================================
%Solves the heat equation for mixed boundary condition
%BC: phi(0,t)=phi_w, dphi/dx(L,t)=0
%IC: phi=phi_w for 0<x<L/10 and 0 for L/10<x<L
%Modified for diffusion of oxygen in water
%=======================================================

close all, clear all

L = 0.2;                %domain length mm
x = linspace(0,L,3001); %domain discretization
t=0:10000:(100*86400);  %time steps note: 86400s in a day
alpha = 3.24e-9;        %m^2/s diffusivity of O2 in water at 40 deg C

phi_w = 5.931e-3;       %kg/m^3 Left boundary condition - concentration
n = 1:1:5000;           %number of fourier terms
k = (2*n - 1)/2;        %definition for convinience

%temporal part of the solution
lam_sq = (k*pi/L).^2;            %eigen val squared
H = exp(-alpha * (t' * lam_sq)); %unsteady term of solution

%spatial part of the solution
Bn = (-2*phi_w./(k*pi)) .* (-cos(k*pi)+cos(k*pi/20)); %fourier co-eff
theta = (k*pi/L)' * x;

f_x = phi_w + Bn*sin(theta);           %initial condition
f_xt = phi_w+(H*diag(Bn)*sin(theta));  %final solution


%plot Tav_w vs time
for i=1:size(f_xt,1)
phi_av_w(i) = trapz(x,f_xt(i,:))/(L*phi_w);
end
figure('Name','Tav_w')
plot(t/3600,phi_av_w,'LineWidth',2)
xlabel('time (hrs)','FontSize',20)
ylabel('\phi_{av} / \phi_w','FontSize',20)
ylim([0 1])

%plot transient solution
figure('Name','Solution')
for i=1:size(f_xt,1)

    plot(x,f_xt(i,:),'LineWidth',2)
    ylim([0 1.1*phi_w]);
    txt1 = ['t = ',num2str(t(i)/3600), ' hrs'];
    txt2 = ['\phi_{av}/\phi_w = ',num2str(phi_av_w(i))];
    tx1 = text(x(end)/2,0.9*phi_w,txt1);
    tx2 = text(x(end)/2,0.8*phi_w,txt2);
    xlabel('x (m)','FontSize',20)
    ylabel('\phi (kg/m^3)','FontSize',20)
    
    if i==1
        pause(1.5)
    else
        pause(0.01)
    end

end



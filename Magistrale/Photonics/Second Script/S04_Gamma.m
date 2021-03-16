clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIELD AND SYSTEM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% mostrare l'effetto delle perdite --> z_eff


T=1500e-15;              %T_max dell'intervallo campionato in s
T0=80e-15;              %durata dell'impulso in s
s0=T0/2.355*sqrt(2);    %conversion FWHM to sigma... why sqrt(2)?!?
tp0=0;                  %posizione iniziale dell'impulso in s

% alpha=0*1/22000;    
alpha=0;              % esagerated absortpion
C=1*-2;%-2;             %initial chirp 
gamma=3e3;
I0=3e-2;

SaveVideo=0;            % controllo per salvare un video
Lnl=1/(gamma*I0);    %lunghezza nonlineare
N=1024*4;               % numero di tempi campionati
nsteps=40;              % numero di passi per la propagazione in z
zmax=0.08;              %distanza massima all'interno del materiale dispersivo 
dz=zmax/nsteps;
Leff=@(z) (1-exp(-alpha*z))/alpha;

%main program
t=linspace(-T,T,N);     % creo il vettore dei tempi


f=zeros(size(t));       % creo vettore delle frequenze nullo da riempire
dF=1/(2*T);
for i=1:N/2
    f(i)=(i-1)*dF;
end
for i=(N/2+1):N
    f(i)=(i-N-1)*dF;
end


Omega=2*pi*f;
all_intensities=zeros(N,nsteps+1);
zplot=zeros(1,nsteps+1);

A=sqrt(I0)*exp(-(1+1i*C).*((t-tp0).^2/(2*s0^2)));
% A=sqrt(I0)*exp(-(1+1i*C).*((t-tp0).^4/(2*s0^4)));
%figure,plot(t,real(A.*exp(-1i*1e14.*t)))
% A=1..*(abs(t-tp0)/(2*s0)));
% A=sech(t/T0).*exp(-1i*C*(t-tp0).^2/(2*T0^2));
% w=T0; %width of rectangle
% A=rectpuls(t,2*w);
all_intensities(1:N,1)=abs(A).^2;

I_max=max(abs(A).^2);
tt2=t(1:N/2);
tt1=t(N/2:end);
[nt2 mt2]=min(abs(abs(A(1:N/2)).^2-0.5*I_max));
[nt1 mt1]=min(abs(abs(A(N/2:end)).^2-0.5*I_max));
B=fft(A);[n_sortf m_sortf]=sort(f);
all_spectra(1:N,1)=abs(B(m_sortf)).^2;


figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

if SaveVideo
    aviobj2=VideoWriter('Video_Dispersion2.avi');
    aviobj2.FrameRate =10;
    open(aviobj2)
end

FWHM0=tt1(mt1)-tt2(mt2);
Broadening(1)=1;

for iz=1:nsteps
    A=A.*exp(-1/2*alpha*0.5*dz);
    nonlinear_propagator=exp(1i.*(abs(A).^2).*gamma*dz);
    A=A.*nonlinear_propagator;
    A=A.*exp(-1/2*alpha*0.5*dz);
    z=(iz-1)*dz;
    
    fase=unwrap(angle(A));
    chirp=-ifft(fft(fase).*1i.*Omega);

        
    
    
    figure(1)
    subplot(4,1,1)
    yyaxis left
    plot(t,abs(A).^2);
    xlabel('T [s]')
    ylabel('Intensity')
    ylim([0 I0])
    yyaxis right
    plot(t,abs(A).^2);
    ylim([0,I0])
    xlim([-5*T0, 5*T0])
    plot(t,real(A.*exp(-2*pi*1i*4e13.*t)))
    ylim([-sqrt(I0),sqrt(I0)])
    ylabel('Real(A)')
    %title(dz*iz)
    title(['z=' num2str(zplot(iz)) ', $\frac{z}{L_{nl}\pi}$=' num2str(zplot(iz)/Lnl/pi)],'interpreter','latex')
    %title(dz*iz/Ld)
    subplot(4,1,2)
    B=fft(A);
    plot(f(m_sortf),abs(B(m_sortf)).^2);
    xlabel('f [Hz]')
    ylabel('Intensity')
    xlim([-3e13,3e13])
    ylim([0,3e4*I0])
    subplot(4,1,3)
    plot(t,fase);
    xlim([-5*T0, 5*T0])
    
    xlabel('T [s]')
    ylabel('Phase')
    subplot(4,1,4)
    plot(t,real(chirp));
    xlabel('T [s]')
    ylabel('Chirp')
     xlim([-5*T0, 5*T0])
     ylim([-1.4e14 1.4e14])

    set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
    set(findobj(gcf,'type','line'),'LineWidth',2)
    pause(.1)
    drawnow
    
    I_max=max(abs(A).^2);
    [nt2 mt2]=min(abs(abs(A(1:N/2)).^2-0.5*I_max));
    [nt1 mt1]=min(abs(abs(A(N/2:end)).^2-0.5*I_max));
    FWHMs=tt1(mt1)-tt2(mt2);
    Broadening(iz+1)=FWHMs/FWHM0;
    
    if SaveVideo
        F=getframe(gcf);
        writeVideo(aviobj2,F);
    end
    
    %   pause
    x=dz*iz;
    all_intensities(:,iz+1)=abs(A).^2;
    all_spectra(:,iz+1)=abs(B(m_sortf)).^2;
    zplot(iz+1)=zplot(iz)+dz;
    %     text(1100,1000,['t=' num2str(x_position(i)) ' ps'])
    
    
end
% close(fig);
if SaveVideo
    close(aviobj2)
end


z=zplot./Lnl;
figure(3)
    pcolor(zplot/Lnl, t*1e12, all_intensities)
    xlabel('z/L_{nl}')

shading interp
ylabel('t (ps)')
box on
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
set(findobj(gcf,'type','line'),'LineWidth',2)


Col=jet(length(zplot));
figure(4),hold on
for k=1:length(zplot)
    plot(t*1e12, all_intensities(:,k),'linewidth',2,'color',Col(k,:))
end
xlabel('t (ps)')
ylabel('Intensity (AU)')
legend(num2str(round(zplot'*1e2)/1e2))
box on
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
set(findobj(gcf,'type','line'),'LineWidth',2)


Col=jet(length(zplot));
figure(5),hold on
for k=1:length(zplot)
    plot(f(m_sortf)/1e12, all_spectra(:,k),'linewidth',2,'color',Col(k,:))
end
xlabel('\nu (THz)')
xlim([-50 50])
ylabel('Intensity (AU)')
legend(num2str(round(z'*1e2/pi)/1e2))
box on
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
set(findobj(gcf,'type','line'),'LineWidth',2)


figure(6)
    pcolor(zplot/Lnl, f(m_sortf)/1e12, all_spectra)
    xlabel('z/L_{nl}')
ylim([-60 60])
shading interp
ylabel('\nu (THz)')
box on
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',14)
set(findobj(gcf,'type','line'),'LineWidth',2)


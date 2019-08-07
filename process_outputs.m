N = 32;
tfinal = 30000;
tau = 1;
hbar = 1;

state = 18;
    
H = zeros(N,N,tfinal/tau + 1);
for k = 1:N
    for l = 1:N
        filename = ['H' num2str(k) ',' num2str(l) '.mat'];
        filename = ['outputs/' filename];
        try
            load(filename);
            H(k,l,:) = Expression1;
            fprintf([filename ' found.\n']);
        catch
%             fprintf([filename ' not found.\n']);
        end
    end
end

eigs = zeros(N,tfinal/tau + 1);
vecs = zeros(N, N,tfinal/tau + 1);
for k = 1:(tfinal/tau + 1)
    [V D] = eig(H(:,:,k));
%     [V D] = eig(H(:,:,k));
%     eigsreal = real(diag(D));
%     [eigs(:,k),order] = sort(eigsreal);
    eigs(:,k) = diag(D);
%     vecs(:,:,k) = V(:,order);
    vecs(:,:,k) = V;
end


%EIGENENERGIES not zoomed
figure('rend','painters','pos',[0 0 750 500])
trange = 0:tau:tfinal;
theplot = plot(trange, eigs,'LineWidth',1.5);
xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
ylabel('$E_{n,\mathbf{k}_{\mathbf{\Gamma}}}(t)$ (Ryd)','Interpreter','latex','FontName','CMU Serif');
xlim([0 tfinal/3]);
ylim([-0.65 0.65]);
title('Eigenenergies at $\mathbf{\Gamma}$ Point','Interpreter','latex')
set(gca,'FontSize',20)
set(gcf,'PaperSize',[9.7 7])
% legend(theplot(26:31), {'26','27','28','29','30','31'},'Orientation','horizontal');


%EIGENENERGIES zoomed
figure('rend','painters','pos',[0 0 750 500])
trange = 0:tau:tfinal;
theplot = plot(trange, eigs,'LineWidth',1.5);
xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
ylabel('$E_{n,\mathbf{k}_{\mathbf{\Gamma}}}(t)$ (Ryd)','Interpreter','latex','FontName','CMU Serif');
xlim([0 tfinal/6]);
ylim([0.18 0.22]);
title('Cluster of Eigenenergies at $\mathbf{\Gamma}$ Point','Interpreter','latex')
set(gca,'FontSize',20)
set(gcf,'PaperSize',[10.5 7])
legend(theplot(17:25), {'17','18','19','20','21','22','23','24','25'},'Orientation','horizontal');


psi0 = vecs(:,state,1);

save('H.mat','H');
save('psi0.mat','psi0');

load('td.mat');
td = permute(Expression1, [2 1]);

vecstr = permute(vecs, [2 1 3]);

aboprob = zeros(tfinal/tau + 1, N);
for k = 1:(tfinal/tau) + 1
    aboprob(k,:) = abs(conj(vecstr(:,:,k))*td(:,k)).^2;
end

figure('rend','painters','pos',[0 0 750 500])
trange = 0:tau:tfinal;
theplot = plot(trange, aboprob,'LineWidth',1.5);
xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
ylabel('$P_{n,\mathbf{k}_{\mathbf{\Gamma}}}(t)$','Interpreter','latex');
title('Overlap Probabilties at $\mathbf{\Gamma}$ Point, Nearly Degenerate State','Interpreter','latex','FontName','CMU Serif')
xlim([0 2000]);
ylim([0 1]);
set(gca,'FontSize',20)
set(gcf,'PaperSize',[9.7 7])
legend(theplot(17:25), {'17','18','19','20','21','22','23','24','25'},'Orientation','horizontal');
% legend(theplot(1), {'1'},'Orientation','horizontal');

% figure('rend','painters','pos',[0 0 750 500])
% trange = 0:tau:tfinal;
% theplot = plot(trange, aboprob,'LineWidth',1.5);
% xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
% ylabel('$P_{n,\mathbf{k}_{\mathbf{\Gamma}}}(t)$','Interpreter','latex');
% title('Overlap Probabilties at $\mathbf{\Gamma}$ Point, Ground State','Interpreter','latex','FontName','CMU Serif')
% xlim([0 2000]);
% ylim([0 1.05]);
% set(gca,'FontSize',20)
% set(gcf,'PaperSize',[9.7 7])
% % legend(theplot(17:25), {'17','18','19','20','21','22','23','24','25'},'Orientation','horizontal');
% legend(theplot(1), {'1'},'Orientation','horizontal');

load('high_sym.mat');
U = Expression1;

dboprob = zeros(tfinal/tau + 1, N);
for k = 1:(tfinal/tau) + 1
    dboprob(k,:) = abs(conj(U(:,:))*td(:,k)).^2;
end

figure('rend','painters','pos',[0 0 750 500])
trange = 0:tau:tfinal;
theplot = plot(trange, dboprob,'LineWidth',1.5);
xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
ylabel('$P_{n,\mathbf{k}_{\mathbf{\Gamma}}}(t)$','Interpreter','latex');
title('Diabatic Overlap Probs. at $\mathbf{\Gamma}$ Point, Nearly Degenerate State','Interpreter','latex','FontName','CMU Serif')
xlim([0 tfinal]);
ylim([0 1]);
set(gca,'FontSize',20)
set(gcf,'PaperSize',[9.7 7])
legend(theplot(17:25), {'17','18','19','20','21','22','23','24','25'},'Orientation','horizontal');


DAboprob = zeros(tfinal/tau + 1, N);
for k = 1:(tfinal/tau) + 1
    DAboprob(k,:) = abs(conj(vecstr(:,:,1))*td(:,k)).^2;
end

figure('rend','painters','pos',[0 0 750 500])
trange = 0:tau:tfinal;
theplot = plot(trange, DAboprob,'LineWidth',1.5);
xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
ylabel('$P_{n,\mathbf{k}_{\mathbf{\Gamma}}}(t)$','Interpreter','latex');
title('Diabatic Overlap Probs. at $\mathbf{\Gamma}$ Point, Nearly Degenerate State','Interpreter','latex','FontName','CMU Serif')
xlim([0 tfinal]);
ylim([0 1]);
set(gca,'FontSize',20)
set(gcf,'PaperSize',[9.7 7])
legend(theplot(17:25), {'17','18','19','20','21','22','23','24','25'},'Orientation','horizontal');




% figure('rend','painters','pos',[0 0 900 450])
% trange = 0:tau:tfinal;
% theplot = plot(trange, abs(td(2,:)').^2,'LineWidth',1.5);
% hold on
% plot(trange, squeeze(abs(vecs(2,state,:))).^2, 'LineWidth', 1.5);
% xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
% ylabel('$|\langle \phi_2(t)|\psi(t) \rangle|^2$','Interpreter','latex');
% title('Atomic Projection Probabilities','Interpreter','latex','FontName','CMU Serif')
% xlim([0 10000]);
% ylim([0 0.2]);
% set(gca,'FontSize',20)
% set(gcf,'PaperSize',[9.7 7])
% legend({'$|\langle \phi_{2}(t)|\psi_{TDSE}(t) \rangle|^2$','$|\langle \phi_2(t)|\psi_{ABO}(t) \rangle|^2$'},'Interpreter','latex','Orientation','horizontal');



% corrabo = zeros(tfinal/tau + 1, N);
% for k = 1:(tfinal/tau) + 1
%     corrabo(k,:) = abs(conj(vecstr(state,:,1))*vecs(:,state,k)).^2;
% end

% figure
% plot(trange,corrabo,'LineWidth',1.5)
% xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
% ylabel('$|\langle \psi_{ABO}(0)|\psi_{ABO}(t) \rangle|^2$','Interpreter','latex');
% title('ABO validity probability for $N = 36$ Atoms','Interpreter','latex')

% normtd = zeros(tfinal/tau + 1, 1);
% for k = 1:(tfinal/tau) + 1
%     normtd(k) = abs(td(:,k)'*td(:,k));
% end
% 
% figure
% plot(trange,normtd,'LineWidth',1.5)
% xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
% ylabel('$|\langle \psi_{TDSE}(t)|\psi_{TDSE}(t) \rangle|^2$','Interpreter','latex');
% title('Normalization $N = 36$ Atoms','Interpreter','latex')

corrtd = zeros(tfinal/tau + 1, 1);
for k = 1:(tfinal/tau) + 1
    corrtd(k) = abs(td(:,1)'*td(:,k));
end

figure('rend','painters','pos',[0 0 750 500])
trange = 0:tau:tfinal;
theplot = plot(trange, corrtd,'LineWidth',1.5);
hold on
xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
ylabel('$A_{18,\mathbf{k}_{\mathbf{\Gamma}}}(t)$','Interpreter','latex');
title('$A_{18,\mathbf{k}_{\mathbf{\Gamma}}}(t)$, at $\mathbf{\Gamma}$ Point, Nearly Degenerate State','Interpreter','latex')
xlim([0 tfinal]);
ylim([0 1.05]);
set(gca,'FontSize',20)
set(gcf,'PaperSize',[9.7 7])


figure('rend','painters','pos',[0 0 750 500])
trange = 0:tau:tfinal;
theplot = plot(trange, corrtd,'LineWidth',1.5);
hold on
xlabel('Time ($\hbar$/Ryd)', 'Interpreter', 'latex');
ylabel('$A_{1,\mathbf{k}_{\mathbf{\Gamma}}}(t)$','Interpreter','latex');
title('$A_{1,\mathbf{k}_{\mathbf{\Gamma}}}(t)$, at $\mathbf{\Gamma}$ Point, Ground State','Interpreter','latex')
xlim([0 tfinal]);
ylim([0 1.05]);
set(gca,'FontSize',20)
set(gcf,'PaperSize',[9.7 7])


% 
% % state = 1;
% auto = zeros(tfinal/tau + 1, 1);
% for k = 1:(tfinal/tau) + 1
%     auto(k) = td(:,1)'*td(:,k);
% %     auto(k) = vecs(:,state,1)'*S(:,:,k)*vecs(:,state,k);
% end
% auto = auto;%./normtd
% 
% 
% startsmooth = 9967;
% 
% smoother = atan(-(trange - tfinal)/100)./(pi/2)+1e-9;
% smoother = smoother(startsmooth:end)';
% 
% % smoother = -(1/tfinal)*trange + 1;
% % smoother = smoother(startsmooth:end)';
% 
% auto(startsmooth:end) = smoother;
% 
% 
% Y = ifft(auto);
% L = length(auto);
% energy = (2*pi/tau)*(0:(L/2))/L;
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% figure('rend','painters','pos',[0 0 900 450])
% theplot = plot(energy, P1,'LineWidth',1.5);
% hold on
% xlabel('Energy (Ryd)', 'Interpreter', 'latex');
% ylabel('Spectrum (arb. units)','Interpreter','latex');
% title('Amplitude Spectral Density','Interpreter','latex')
% xlim([0 1.5]);
% % ylim([0 1]);
% set(gca,'FontSize',20)
% set(gcf,'PaperSize',[9.7 7])

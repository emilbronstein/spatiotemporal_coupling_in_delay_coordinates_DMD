clear; close all; clc

% ----- time ----- %
dt = 0.1; % timestep delta t
t = (0:dt:10); % time vector

% ----- frequencies ----- %
w = [3 , 5]; % system's natural frequencies

% ----- signal ----- %
x = sin(w(1)*t) + sin(w(2)*t); % the two-mode sine signal in Eq. (47)
x = x + 1e-12*randn(size(x)); % addition of a white standard Gaussian noise
    
% ------- data augmentation and DMD as in step 1 of Algorithm 1 ------- %
s = 11; % augmentation number
[X,Xprime] = data_augmentation(x,s);
[lambda,phi] = exact_DMD(X,Xprime);

[wl,sortind] = w_from_lambda(lambda,dt); % the oscillation frequencies as in Eq. (13) and according to step 2 of Algorithm 1

% sort the eigenvalues and eigenvectors according to the order of the oscillation
% frequencies
lambda = lambda(sortind); 
phi = phi(:,sortind);

% ----- Figure 4 in the paper ----- %

figure('Color',[1 1 1]);
font = 14;
row = 1; col = 2;

subplot(row,col,1)
th = linspace(0,2*pi,1e5);
ucx = cos(th); ucy = sin(th);
plot(real(lambda(1:4)),imag(lambda(1:4)),'ob','markersize',10); hold on
plot(real(lambda(5:end)),imag(lambda(5:end)),'xr','markersize',10)
hold on
plot(ucx,ucy,'-k','linewidth',2)
xlabel('$\Re(\hat{\lambda}_l$)','interpreter','latex','FontSize',font)
ylabel('$\Im(\hat{\lambda}_l$)','interpreter','latex','FontSize',font)
set(gca,'TickLabelInterpreter','latex','FontSize',font);
axis([-1.1 1.1 -1.1 1.1])
axis square
text(-0.1,1.1,'(a)','Units','normalized','interpreter','latex','FontSize',font);
legend('True','Spurious','Unit circle','interpreter','latex','FontSize',font,'location','northoutside','box','off')

subplot(row,col,2)
plot((1:4),wl(1:4),'ob','markersize',10); hold on
plot((5:length(wl)),wl(5:end),'xr','markersize',10);
xlabel('$l$','interpreter','latex','FontSize',font)
ylabel('$\omega_l$','interpreter','latex','FontSize',font)
set(gca,'TickLabelInterpreter','latex','FontSize',font);
set(gca,'XTick',(1:2:length(wl)-1));
set(gca,'YTick', [3,(5:5:30)] );
text(1,1,'(b)','Units','normalized','interpreter','latex','FontSize',font);
axis([0 13 2 31])

%%%%%

ws = (2*pi)/dt; % sampling frequency
Ny = floor( (0.5*ws)./wl ); % Nyquist-Shannon sampling criterion, as in Eq. (40)
% and according to step 3 of Algorithm 1

% repeatedly computed eigenvalues as in Eq. (41)
for l=1:s
    lambda_tilde(l,:) = ( phi(l+1,:) ./ phi(1,:) ) .^ (1/l);
end

% ----- Figure 5 in the paper ----- %

lambda1 = lambda_tilde(:,1);
lambda3 = lambda_tilde(:,3);
lambda9 = lambda_tilde(:,9);


figure('Color',[1 1 1])
font = 12;

% unit circle
th = linspace(0,2*pi,1e5);
ucx = cos(th); ucy = sin(th);

% lambda1
subplot('Position',[0.075 0.4 0.25 0.25]);
plot(real(lambda1(1:Ny(1))),imag(lambda1(1:Ny(1))),'ob','markersize',8); hold on
plot(real(lambda1(Ny(1)+1:end)),imag(lambda1(Ny(1)+1:end)),'.b','markersize',15); hold on
plot(ucx,ucy,'-k','linewidth',2)
title('True','interpreter','latex','FontSize',font)
xlabel('$\Re(\tilde{\lambda}_1$)','interpreter','latex','FontSize',font)
ylabel('$\Im(\tilde{\lambda}_1$)','interpreter','latex','FontSize',font)
set(gca,'TickLabelInterpreter','latex','FontSize',font);
axis([-1.1 1.1 -1.1 1.1])
axis square
text(1,1,'(a)','Units','normalized','interpreter','latex','FontSize',font);

% lambda3
subplot('Position',[0.4 0.4 0.25 0.25]);
plot(real(lambda3(1:Ny(3))),imag(lambda3(1:Ny(3))),'ob','markersize',8); hold on
plot(real(lambda3(Ny(3)+1:end)),imag(lambda3(Ny(3)+1:end)),'.b','markersize',15); hold on
plot(ucx,ucy,'-k','linewidth',2)
title('True','interpreter','latex','FontSize',font)
xlabel('$\Re(\tilde{\lambda}_3$)','interpreter','latex','FontSize',font)
ylabel('$\Im(\tilde{\lambda}_3$)','interpreter','latex','FontSize',font)
set(gca,'TickLabelInterpreter','latex','FontSize',font);
axis([-1.1 1.1 -1.1 1.1])
axis square
text(1,1,'(b)','Units','normalized','interpreter','latex','FontSize',font);

% lambda9
subplot('Position',[0.725 0.4 0.25 0.25]);
plot(real(lambda9),imag(lambda9),'xr','markersize',8); hold on
hold on
plot(ucx,ucy,'-k','linewidth',2)
title('Spurious','interpreter','latex','FontSize',font)
xlabel('$\Re(\tilde{\lambda}_{9}$)','interpreter','latex','FontSize',font)
ylabel('$\Im(\tilde{\lambda}_{9}$)','interpreter','latex','FontSize',font)
set(gca,'TickLabelInterpreter','latex','FontSize',font);
axis([-1.1 1.1 -1.1 1.1])
axis square
text(1,1,'(c)','Units','normalized','interpreter','latex','FontSize',font);

%%%%%

% compute the averaged eigenvalues as in Eq. (42) and according to step 4
% in Algorithm 1 in the paper
for i=1:length(Ny)
    if Ny(i) > s
        lambda_tilde_avg(i,1) = mean(lambda_tilde(1:end,i));
    else    
        lambda_tilde_avg(i,1) = mean(lambda_tilde(1:Ny(i),i));
    end
end

% compute the absolute errors as in Eq. (43) and according to step 5 in
% Algorhitm 1 in the paper
epsilon = abs(lambda_tilde_avg - lambda);

% step 6 in Algorhitm 1 in the paper
[idx,idx_true] = true_index(epsilon);

% ----- Figure 6 in the paper ----- %

figure('Color',[1 1 1]); font = 14;
ind = (1:s+1);
semilogy(ind(idx == idx_true),epsilon(idx == idx_true),'ob','markersize',10); hold on
semilogy(ind(idx ~= idx_true),epsilon(idx ~= idx_true),'xr','markersize',10);
xlabel('$l$','interpreter','latex','FontSize',font)
ylabel('$\varepsilon_l$','interpreter','latex','FontSize',font)
set(gca,'TickLabelInterpreter','latex','FontSize',font);
axis([0 13 1e-16 1e-3])
set(gca,'xtick',(1:12))
set(gca,'ytick',[1e-15, 1e-10, 1e-5, 1])
legend('True','Spurious','interpreter','latex','FontSize',font,'location','east','box','off')

% ----- Figure 3 in the paper ----- %

l = diag(lambda_tilde_avg(ind(idx == idx_true))); % the true averaged computed eigenvalues in Eq. (46) in the paper
p = phi(:,ind(idx == idx_true)); % the true zeroth sub-modes in Eq. (46) in the paper
x0 = X(:,1); % the relevant initial observations
sigma0 = p\x0; % \sigma_l in Eq. (46) in the paper

% reconstruction of the observations according to Eq. (46) in the paper
for k=1:length(t)
    xr(:,k) = p*(l.^(t(k)/dt))*sigma0;
end

figure('Color',[1 1 1]); font = 14;
line2 = plot(t(2:end),real(xr(1,2:end)),'b--','LineWidth',3); hold on
line1 = plot(t,x,'k-','LineWidth',1); hold on
xlabel('$t$','interpreter','latex','FontSize',font)
ylabel('$x$','interpreter','latex','FontSize',font)
set(gca,'TickLabelInterpreter','latex','FontSize',font);
legend([line1,line2],'Original','Reconstructed','interpreter','latex','FontSize',font,'location','southwest','box','off')
axis([-0.1 10.1 -2.5 2])



% ----- functions used in this script ----- %

function [X,Xprime] = data_augmentation(q,s)

% data_augmentation constructs the augmented observation matrices \hat{X}
% and \hat{X}' as in Eq. (14) in the paper

%%% input %%%
% q: input signal
% s: augmentation number

%%% output %%%
% X,Xprime: augmented observation matrices \hat{X} and \hat{X}' as in Eq. (14) in the paper

q_aug = zeros( size(q,1)*(s+1) , size(q,2)-s );

for i=1:size(q,2)-s
    v = q(:,i:i+s);
    v = v(:);
    q_aug(:,i) = v;
end

X = q_aug(: , 1:end-1);
Xprime = q_aug(: , 2:end);

end

function [lambda_DMD,phi_DMD] = exact_DMD(X,Xprime)

% exact_DMD executes the exact DMD algorithm on the input matrices X and
% Xprime.

%%% input %%%
% X,Xprime: augmented observation matrices \hat{X} and \hat{X}' as in Eq. (14) in the paper

%%% output %%%
% lambda_DMD: The exact DMD eigenvalues
% phi_DMD: The exact DMD modes

[U,S,V] = svd(X,'econ');
r = rank(X);
U = U(:, 1:r); % truncate to rank r
S = S(1:r, 1:r);
V = V(:, 1:r);
A_tilde = U' * Xprime * V * S^(-1);
[omega_DMD,lambda_DMD] = eig(A_tilde);
phi_DMD = Xprime * V * S^(-1) * omega_DMD * lambda_DMD^(-1);
lambda_DMD = diag(lambda_DMD);

end

function [w,sortind] = w_from_lambda(lambda,dt)

% w_from_lambda extracts the natural frequencies that correspond to the DMD
% eigenvalues according to Eq. (13) in the paper.
% Additionally, the natural frequencies are sorted in ascending
% order for convenience.

%%% input %%%
% lambda: DMD eigenvalues
% dt: timestep

%%% output %%%
% extracted and sorted natural frequencies and their sorting indices (sortind)

w = abs(log(lambda)) / dt; % Eq. (13) in the paper
[w,sortind] = sort(w);

end

function [idx,idx_true] = true_index(epsilon)

% true_index returns the true indices corresponding to the DMD components

%%% input %%%
% epsilon: the absoulte errors as in Eq. (43) in the paper

%%% output %%%
% idx: the indices of the two subsets of application of k-means with k=2
% idx_true: true index


idx = kmeans(log10(epsilon),2); % indices of the two subsets

% the average absolute errors of the two subsets
epsilon_avg = zeros(1,2);
for i=1:length(epsilon_avg)
    epsilon_avg(i) = mean( epsilon( idx == i ) );
end

idx_true = find(epsilon_avg == min(epsilon_avg)); % true index

end

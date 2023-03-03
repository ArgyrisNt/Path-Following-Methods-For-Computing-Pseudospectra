%% Cobra algorithm to compute pseudospectra

% lamda0  : Initial eigenvalue
% epsilon : Perturbation value
% A       : Matrix
% tol2    : Tolerance for needed points
% h1      : Step of prediction 
% H       : Length of Cobra neck
% flag    : 1 : Plot errors
%           2 : Plot boundary of pseudospectrum with errors in 3D
%           3 : Plot boundary of pseudospectrum

function Cobra(lamda0,epsilon,A,tol2,h1,H,flag)

% INPUTS
K = 6000; % The number of points on boundary of pseudpospectrum
d0 = 1; % Direction to look for the first point on boundary
tol = 1e-7; % Tolerance value to terminate while loop
show = 0; % Set equals 1 to see results while solving

spectrum = eig(A);

% STEP 0: Compute the first point z_piv of the boundary
hold on
[n,m] = size(A);
I = eye(n);
theta0 = epsilon; % Set initial step 
z_piv_new = lamda0+theta0*d0; % Initialize first pivot point
s_new = svd(z_piv_new*I-A); % Find singular values of matrix z_piv*I-A
s_new_min = min(s_new); % Find minimum singular value
s_old_min = s_new_min;
z_piv_all = [];
z_piv_all = [z_piv_all z_piv_new];

while (abs(s_old_min-epsilon)/epsilon)>tol
	z_piv_old = z_piv_new;
    s_old = svd(z_piv_old*I-A); % Find singular values of matrix z_piv*I-A
    s_old_min = min(s_old); % Find minimum singular value, that is g(z)
	[U,~,V] = svd(z_piv_old*I-A); % Find minimum singular vectors
	u_min = U(:,n); % Minimum singular left vector
	v_min = V(:,m); % Minimum singular right vector
    der_h = real(v_min'*u_min); % h'(theta)
	z_piv_new = z_piv_old-(s_old_min-epsilon)/der_h; % Update z_piv
    z_piv_all = [z_piv_all z_piv_new];
end
z_piv = z_piv_new;
% plot(z_piv_all,'b')

% STEP 1,2: Prediction/Correction and Cobra
r = (1i*(v_min')*u_min)/abs((v_min')*u_min); %i*grad(g(z))/|grad(g(z))|
mi = 8;
h = H/mi;
z_sup = zeros(1,K); % Here I'll store the support points
z = zeros(1,mi);
z_final = zeros(1,K); % Here I'll store the points of boundary
z_final(1) = z_piv; 
k = 2;
error1 = abs((s_old_min-epsilon))/epsilon;
all = [];
all = [all error1];

ic = 1;
while k < K
	z_pr(k) = z_piv+h1*r; % Update prediction
	s_pr = svd(z_pr(k)*I-A); % Find singular values of matrix z_pr*I-A
	s_pr_min = min(s_pr); % Find minimum singular value
	[U_pr,~,V_pr] = svd(z_pr(k)*I-A); % Find minimum singular vectors
	u_min = U_pr(:,n); % Minimum singular left vector
	v_min = V_pr(:,m); % Minimum singular right vector
	z_sup(k) = z_pr(k)-(s_pr_min-epsilon)/((u_min')*v_min); % Update z1
    % Calculate Cobra Neck 
	for j = 1:mi
		ic = ic + 1;
        a1 = real(z_sup(k));
        b1 = imag(z_sup(k));
        a2 = real(z_piv);
        b2 = imag(z_piv);
		zeta_real(j) = a2 +j*h*(a1-a2);
        zeta_imag(j) = b2 +j*h*(b1-b2);
        zeta(j) = zeta_real(j) + 1i*zeta_imag(j);
        s_zeta{j} = svd(zeta(j)*I-A); % Find singular values of zeta(mi)*I-A
        s_zeta_min(j) = min(s_zeta{j}); % Find minimum singular value
        [U_zeta{j},~,V_zeta{j}] = svd(zeta(j)*I-A); % Find minimum singular vectors
        u_zeta_min{j} = U_zeta{j}(:,n); % Minimum singular left vector
        v_zeta_min{j} = V_zeta{j}(:,m); % Minimum singular right vector
        z(j) = zeta(j)-(s_zeta_min(j)-epsilon)/((u_zeta_min{j}')*v_zeta_min{j});
		z_final(ic) = z(j);
        s_final{j} = svd(z(j)*I-A);
		s_final_min(j) = min(s_final{j});
        error(j)=abs(s_final_min(j)-epsilon)/epsilon;
        all = [all error(j)];
    end
    
    z_piv = z(mi); % Choose the next pivot point 
    if (abs(z_final(ic)-z_final(1))<tol2)
        break;
    end
	r = (1i*(v_zeta_min{mi}')*u_zeta_min{mi})/abs((v_zeta_min{mi}')*u_zeta_min{mi}); % Update direction r
    k = k+1;
    if show==1
       if ceil(k/10)*10==k
          clf 
          hold on
          z_final = z_final(z_final~=0);
          plot(real(z_final),imag(z_final),'b')
          pause(0.001)
       end
    end
end
z_final = z_final(z_final~=0);

fprintf('\nRESULTS WITH COBRA');
fprintf('\nTotal steps : %d',k);
fprintf('\nTotal points on curve : %d',ic);
fprintf('\n');

if flag==1
%    plot(all,'b')
%    plot(real(z_final(1:end/2)),log(all(1:end/2)),'b');
   plot(log(all),'b');
end

if flag==2
   figure(1)
   plot3(real(z_final),imag(z_final),all,'b');
   grid on
   view(3)
end

if flag==3
   plot(real(spectrum),imag(spectrum),'r*'); % Plot eigenvalues of A
   hold on
%    plot(real(zeta),imag(zeta),'r*');
%    plot(real(z),imag(z),'g*');
%    plot(real(z(error_min)),imag(z(error_min)),'b*');
   plot(real(z_final),imag(z_final),'b');
%    plot(real(z_final),imag(z_final),'go');
   % Plot initial point of curve
   xlabel('R');
   ylabel('Im');
end
end

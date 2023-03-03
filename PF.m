%% Curve Tracing algorithm to compute pseudospectra

% lamda0  : Initial eigenvalue
% epsilon : Perturbation value
% A       : Matrix
% tol2    : Tolerance for needed points
% tk      : Step of PF
% flag    : 1 : Plot errors
%           2 : Plot boundary of pseudospectrum with errors in 3D
%           3 : Plot boundary of pseudospectrum

function PF(lamda0,epsilon,A,tol2,tk,flag)

% INPUTS
K = 6000; % The number of points on boundary of pseudpospectrum
d0 = 1; % Direction to look for the first point on boundary
tol = 1e-7; % Tolerance value to terminate while loop
show = 0; % Set equals 1 to see results while solving

spectrum = eig(A);

% STEP 0: Compute the first point z1 of the boundary
hold on
[n,m] = size(A);
I = eye(n);
theta0 = epsilon; % Set initial step 
z1_new = lamda0+theta0*d0; % Initialize firsti point z1
s_new = svd(z1_new*I-A); % Find singular values of matrix z1*I-A
s_new_min = min(s_new); % Find minimum singular value
s_old_min = s_new_min;
z1_all = [];
z1_all = [z1_all z1_new];

while (abs(s_old_min-epsilon)/epsilon)>tol
	z1_old = z1_new;
    s_old = svd(z1_old*I-A); % Find singular values of matrix z1*I-A
    s_old_min = min(s_old); % Find minimum singular value, that is g(z)
	[U,~,V] = svd(z1_old*I-A); % Find minimum singular vectors
	u_min = U(:,n); % Minimum singular left vector
	v_min = V(:,m); % Minimum singular right vector
    der_h = real(v_min'*u_min); % h'(theta)
	z1_new = z1_old-(s_old_min-epsilon)/der_h; % Update z1
    z1_all = [z1_all z1_new];
end
z1 = z1_new;
% plot(z1_all,'r')

% STEP 1,2: Prediction and Correction
r = (1i*(v_min')*u_min)/abs((v_min')*u_min); %i*grad(g(z))/|grad(g(z))|
z = zeros(1,K); % Here I'll store the calculated points of boundary
z(1) = z1;
k = 2;
error = (s_old_min-epsilon)/epsilon;
all = [];
all = [all error];

while k < K
	z_pr(k) = z(k-1)+tk*r; % Update prediction
	s_pr = svd(z_pr(k)*I-A); % Find singular values of matrix z_pr*I-A
	s_pr_min = min(s_pr); % Find minimum singular value
	[U_pr,~,V_pr] = svd(z_pr(k)*I-A); % Find minimum singular vectors
	u_min = U_pr(:,n); % Minimum singular left vector
	v_min = V_pr(:,m); % Minimum singular right vector
	z(k) = z_pr(k)-(s_pr_min-epsilon)/((u_min')*v_min); % Update z1
    s_z = svd(z(k)*I-A); % Find singular values of matrix z*I-A
	s_z_min = min(s_z); % Find minimum singular value
    error = abs((s_z_min-epsilon)/epsilon);
    all = [all error];
    if (abs(z(k)-z1)<tol2)
        break;
    end
	r = (1i*(v_min')*u_min)/abs((v_min')*u_min); % Update direction r
    k = k+1;
    if show==1
       if ceil(k/10)*10==k
          clf 
          hold on
          z = z(z~=0);
          plot(real(z),imag(z),'r')
          pause(0.001)
       end
    end
end
z = z(z~=0);

fprintf('\nRESULTS WITH PF');
fprintf('\nTotal steps : %d',k);
fprintf('\nTotal points on curve : %d',k);
fprintf('\n');

if flag==1
%    plot(all,'r')
%    plot(real(z(1:end/2)),log(all(1:end/2)),'r');
   plot(log(all),'r');
end

if flag==2
   figure(1)
   plot3(real(z),imag(z),all,'r');
   grid on
   view(3)
end

if flag==3
   plot(real(spectrum),imag(spectrum),'b*'); % Plot eigenvalues of A
   hold on
%    plot(real(z),imag(z),'bo')
   plot(real(z),imag(z),'r'); % Plot å-pseudospectrum
   xlabel('R');
   ylabel('Im');
end
end

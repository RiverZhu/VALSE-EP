function out = VALSE_EP( y_q, m, ha, x, Iter_max, B, yy_min, alpha, method_EP )
%VALSE algorithm for line spectral estimation
% INPUTS:
%   y_q  - measurement vector 
%   m  - is a vector containing the indices (in ascending order) of the M
%       measurements; subset of {0,1,...,m(end)}
%   ha - indicator determining which approximation of the
%       frequency posterior pdfs will be used:
%           ha=1 will use Heuristic #1
%           ha=2 will use Heuristic #2
%           ha=3 will use point estimation of the frequencies (VALSE-pt)
%   x  - the true signal - used for computing the MSE vs iterations
%
% OUTPUTS:
%   out - structure
%      .freqs      - vector of frequency estimates
%      .amps       - vector of amplitude estimates
%      .x_estimate - reconstructed signal
%      .noise_var  - estimate of the noise variance
%      .iterations - number of iterations until convergence
%      .mse        - evolution of the mse of x_estimate with iterations
%      .K          - evolution of the estimated number of components with iterations
%
M     = size(m,1);
N     = m(M)+1;     % size of full data
L     = N;          % assumed number of components
z_A_ext = zeros(M,1);
v_A_ext = 1e1*ones(M,1);
mse = zeros(Iter_max,1);
Kt = zeros(Iter_max,1);
yI = zeros(N,1);
if B<inf
    y_pre = yy_min + (y_q+0.5)* alpha;
    y_pre_c = y_pre(1:end/2)+1j*y_pre(end/2+1:end);
    yI(m+1) = y_pre_c;
else
    y=y_q;
    yI(m+1) = y;
end



R  = yI*yI';
sR = zeros(N-1,1);
for i=2:N
    for k=1:i-1
        sR(i-k) = sR(i-k) + R(i,k);
    end
end
Rh  = toeplitz([sum(diag(R));sR])/N; 
try
    evs = sort(real(eig(Rh)),'ascend');
catch 
     warning('Assigning a value of NaN.');
end
nu  = mean(evs(1:floor(N/4)));

if B<inf
    z_A_ext_real = [real(z_A_ext);imag(z_A_ext)];
    v_A_ext_real = [v_A_ext;v_A_ext]/2;
    [z_B_post_real, v_B_post_real] = GaussianMomentsComputation_MJH(y_q, z_A_ext_real, v_A_ext_real, yy_min, B, alpha, nu/2);
    v_B_post = v_B_post_real(1:end/2)+v_B_post_real(end/2+1:end);
    z_B_post = z_B_post_real(1:end/2)+1j*z_B_post_real(end/2+1:end);

    
    switch method_EP
        case 'scalar_EP'
            v_B_post = mean(v_B_post)*ones(M,1);
        case 'diag_EP'
%             v_B_post = v_B_post;
    end
    
    v_B_ext = v_B_post.*v_A_ext./(v_A_ext-v_B_post);
    z_B_ext = v_B_ext.*(z_B_post./v_B_post-z_A_ext./v_A_ext);
    
else
    v_B_ext = nu*ones(M,1);
    z_B_ext = y_q;
end


y = z_B_ext;
sigma = v_B_ext;


y2    = (y./sqrt(sigma))'* (y./sqrt(sigma));
A     = zeros(L,L);
J     = zeros(L,L);
h     = zeros(L,1);
w     = zeros(L,1);
C     = zeros(L);
t = 1;

% Initialization of the posterior pdfs of the frequencies
res   = y./sqrt(sigma);
for l=1:L
    % noncoherent estimation of the pdf
    yI = zeros(N,1);
    yI(m+1) = res;
    R  = yI*yI';
    sR = zeros(N-1,1);
    for i=2:N
        for k=1:i-1
            sR(i-k) = sR(i-k) + R(i,k);
        end
    end
    if l==1 % use the sample autocorrelation to initialize the model parameters
        K   = floor(L/2);
        rho = K/L;
%         tau = (y2/sum(1./sigma)-sum(sigma)/M)/(rho*L);
        tau = (y2-M)/(rho*L)/sum(1./sigma);
    end
    etaI   = 2*sR/(M+mean(sigma)/tau);
    ind    = find(abs(etaI)>0);
    if ha~=3
        [~,mu,kappa] = Heuristic2(etaI(ind), ind);
        A(m+1,l) = exp(1i*m * mu) .* ( besseli(m,kappa,1)/besseli(0,kappa,1) );
    else
        [~,mu] = pntFreqEst(etaI(ind), ind);
        A(m+1,l) = exp(1i*m * mu);
    end
    % compute weight estimates; rank one update
    w_temp = w(1:l-1); C_temp = C(1:l-1,1:l-1);
    J(1:l-1,l) = A(m+1,1:l-1)'*diag(1./sigma)*A(m+1,l); J(l,1:l-1) = J(1:l-1,l)'; J(l,l) = sum(1./sigma);
    J = (J+J')/2;
    h(l) = A(m+1,l)'*diag(1./sigma)*y;
    v = 1/(sum(1./sigma)+1/tau-real(J(1:l-1,l)'*C_temp*J(1:l-1,l)));  % diagonal case
    u = v .* (h(l) - J(1:l-1,l)'*w_temp);

    w(l) = u;
    ctemp = C_temp*J(1:l-1,l);
    w(1:l-1) = w_temp - ctemp*u;
    C(1:l-1,1:l-1) = C_temp + v*(ctemp*ctemp');
    C(1:l-1,l) = -v*ctemp;  C(l,1:l-1) = C(1:l-1,l)'; C(l,l) = v;

    % the residual signal
    res = (y - A(m+1,1:l)*w(1:l))./sqrt(sigma);

    if l==K % save mse and K at initialization
        xro    = A(:,1:l)*w(1:l);
        if B==1
            deb_c = xro'*x/(xro'*xro+eps); 
        else
           deb_c = 1;
        end
        mse(t) = 10*log10(norm(x-deb_c*xro)^2/norm(x)^2);
        Kt(t)  = K;
    end
end


%%% Start the VALSE algorithm
cont = 1;
while cont
    
        t = t + 1;
        % Update the support and weights
        [ K, s, w, C ] = maxZ_diag( J, h, M, sigma, rho, tau );
        if K==0
            sprintf('No signal detected');
            out = struct('freqs',th,'amps',w(s),'x_estimate',xr,'noise_var',nu,'iterations',t,'MSE',mse,'K',Kt);
            return;
        end
        % Update the noise variance, the variance of prior and the Bernoulli probability
        if K>0
%                 z_A_post  = A(m+1,s)*w(s);
%                 v_A_post = real(diag(A(m+1,s)*C(s,s)*A(m+1,s)'));
%                 TT_mat = (A(m+1,s).*conj(A(m+1,s)));
%                 add_var1 = w(s)'*w(s)*ones(M,1)-TT_mat*(w(s).*conj(w(s)));
%                 add_var2 = trace(C(s,s))*ones(M,1)-TT_mat*diag(C(s,s));
%                 v_A_post = v_A_post+real(add_var1)+real(add_var2);
%                 sigma = abs(z_B_ext-z_A_post).^2+v_A_post;
            tau = real( w(s)'*w(s)+trace(C(s,s)) )/K;
            if K<L
                rho = K/L;
            else
                rho = (L-1)/L; % just to avoid the potential issue of log(1-rho) when rho=1
            end
        else
            rho = 1/L; % just to avoid the potential issue of log(rho) when rho=0
        end
        inz = 1:L; inz = inz(s); % indices of the non-zero components
        th = zeros(K,1);
        for i = 1:K
            if K == 1
                r = y;
               eta = 2 * ( (r./sigma) * w(inz)' );
            else
                A_i = A(m+1,inz([1:i-1 i+1:end]));
                r = y - A_i*w(inz([1:i-1 i+1:end]));
               eta = 2*diag(1./sigma) * ( r * w(inz(i))' - A_i * C(inz([1:i-1 i+1:end]),i) );
            end
            if ha == 1
                [A(:,inz(i)), th(i)] = Heuristic1( eta, m, 1000 );
            elseif ha == 2
                [A(:,inz(i)), th(i)] = Heuristic2( eta, m );
            elseif ha == 3
                [A(:,inz(i)), th(i)] = pntFreqEst( eta, m );
            end
        end
        
        z_A_post  = A(m+1,s)*w(s);

        v_A_post = real(diag(A(m+1,s)*C(s,s)*A(m+1,s)'));
        TT_mat = (A(m+1,s).*conj(A(m+1,s)));
        add_var1 = w(s)'*w(s)*ones(M,1)-TT_mat*(w(s).*conj(w(s)));
        add_var2 = trace(C(s,s))*ones(M,1)-TT_mat*diag(C(s,s));
        v_A_post = v_A_post+real(add_var1)+real(add_var2);

        switch method_EP
            case 'scalar_EP'
                v_A_post = mean(v_A_post)*ones(M,1);
            case 'diag_EP'
%                 v_A_post = v_A_post;
        end
        
        
        v_A_ext = v_A_post.*v_B_ext./(v_B_ext-v_A_post);
        v_A_ext = v_A_ext.*(v_A_ext>0)+max(v_A_ext)*10*(v_A_ext<=0);
        z_A_ext = v_A_ext.*(z_A_post./v_A_post-y./v_B_ext); 

        if B<inf
                z_A_ext_real = [real(z_A_ext);imag(z_A_ext)];
                v_A_ext_real = [v_A_ext;v_A_ext]/2;
                [z_B_post_real, v_B_post_real] = GaussianMomentsComputation_MJH(y_q, z_A_ext_real, v_A_ext_real, yy_min, B, alpha, nu/2);
                v_B_post = v_B_post_real(1:end/2)+v_B_post_real(end/2+1:end);
                z_B_post = z_B_post_real(1:end/2)+1j*z_B_post_real(end/2+1:end);

                switch method_EP
                    case 'scalar_EP'
                        v_B_post = mean(v_B_post)*ones(M,1);
                    case 'diag_EP'
%                         v_B_post = v_B_post;
                end
                v_B_ext = v_B_post.*v_A_ext./(v_A_ext-v_B_post+eps);
                z_B_ext = v_B_ext.*(z_B_post./v_B_post-z_A_ext./v_A_ext);
                nu = mean(abs(z_B_ext-z_B_post).^2+v_B_post);
        else
               v_B_post = nu.*v_A_ext./(nu+v_A_ext);
               z_B_post = v_B_post.*(z_A_ext./v_A_ext+y_q./nu);
               v_B_ext = v_B_post.*v_A_ext./(v_A_ext-v_B_post);
               z_B_ext = v_B_ext.*(z_B_post./v_B_post-z_A_ext./v_A_ext);  

        
               nu = mean(abs(z_B_ext-z_B_post).^2+v_B_post);
                   

        end
        sigma = v_B_ext;
        y = z_B_ext;
        J = A(m+1,:)'*diag(1./sigma)*A(m+1,:);
        J = J - diag(diag(J)) +  sum(1./sigma)*eye(N);
        h   = A(m+1,:)'*diag(1./sigma)*y;
        
        
        
        xr     = A(:,s)*w(s);
        if B==1
           deb_c = xr'*x/(xr'*xr+eps); 
        else
           deb_c = 1;
        end
        mse(t) = 10*log10(norm(x-deb_c*xr)^2/norm(x)^2);
        Kt(t)  = K;
        
    % stopping criterion:
    % the relative change of the reconstructed signalis below threshold or
    % max number of iterations is reached
        if (norm(xr-xro)/norm(xro)<1e-6) || (norm(xro)==0&&norm(xr-xro)==0) || (t >= Iter_max)||nu>1e40
            cont = 0;
            mse(t+1:end) = mse(t);
            Kt(t+1:end)  = Kt(t);
        end
        xro = xr;
end
out = struct('freqs',th,'amps',w(s),'x_estimate',xr,'noise_var',nu,'iterations',t,'MSE',mse,'K',Kt);
end

function [a, theta, kappa, mu] = Heuristic1( eta, m, D )
%Heuristic1 Uses the mixture of von Mises approximation of frequency pdfs
%and Heuristic #1 to output a mixture of max D von Mises pdfs

M     = length(m);
tmp   = abs(eta);
A     = besseli(1,tmp,1)./besseli(0,tmp,1);
kmix  = Ainv( A.^(1./m.^2) );
[~,l] = sort(kmix,'descend');
eta_q = 0;
for k=1:M
    if m(l(k)) ~= 0
        if m(l(k)) > 1
            mu2   = ( angle(eta(l(k))) + 2*pi*(1:m(l(k))).' )/m(l(k));
            eta_f = kmix(l(k)) * exp( 1i*mu2 );
        else
            eta_f = eta(l(k));
        end
        eta_q = bsxfun(@plus,eta_q,eta_f.');
        eta_q = eta_q(:);
        kappa = abs(eta_q);
        
        % to speed up, use the following 4 lines to throw away components
        % that are very small compared to the dominant one
        kmax  = max(kappa);
        ind   = (kappa > (kmax - 30) ); % corresponds to keeping those components with amplitudes divided by the highest amplitude is larger than exp(-30) ~ 1e-13
        eta_q = eta_q(ind);
        kappa = kappa(ind);
        
        if length(eta_q) > D
            [~, in] = sort(kappa,'descend');
            eta_q   = eta_q(in(1:D));
        end
    end
end
kappa   = abs(eta_q);
mu      = angle(eta_q);
kmax    = max(kappa);
I0reg   = besseli(0,kappa,1) .* exp(kappa-kmax);
Zreg    = sum(I0reg);
n       = 0:1:m(end);
[n1,k1] = meshgrid(n, kappa);
a       = sum( (diag(exp(kappa-kmax))* besseli(n1,k1,1) /Zreg ).*exp(1i*mu*n),1).';
theta   = angle(sum( (diag(exp(kappa-kmax))* besseli(1,kappa,1) /Zreg ).*exp(1i*mu*1),1));
end

function [a, theta, kappa] = Heuristic2( eta, m )
%Heuristic2 Uses the mixture of von Mises approximation of frequency pdfs
%and Heuristic #2 to output one von Mises pdf

N     = length(m);
ka    = abs(eta);
A     = besseli(1,ka,1)./besseli(0,ka,1);
kmix  = Ainv( A.^(1./m.^2) );
k     = N;
eta_q = kmix(k) * exp( 1i * ( angle(eta(k)) + 2*pi*(1:m(k)).' )/m(k) );
for k = N-1:-1:1
    if m(k) ~= 0
        phi   = angle(eta(k));
        eta_q = eta_q + kmix(k) * exp( 1i*( phi + 2*pi*round( (m(k)*angle(eta_q) - phi)/2/pi ) )/m(k) );
    end
end
[~,in] = max(abs(eta_q));
mu     = angle(eta_q(in));
d1     = -imag( eta' * ( m    .* exp(1i*m*mu) ) );
d2     = -real( eta' * ( m.^2 .* exp(1i*m*mu) ) );
if d2<0 % if the function is locally concave (usually the case)
    theta  = mu - d1/d2;
    kappa  = Ainv( exp(0.5/d2) );
else    % if the function is not locally concave (not sure if ever the case)
    theta  = mu;
    kappa  = abs(eta_q(in));
end
n      = (0:1:m(end))';
a      = exp(1i*n * theta).*( besseli(n,kappa,1)/besseli(0,kappa,1) );
end

function [a, theta] = pntFreqEst( eta, m )
%pntFreqEst - point estimation of the frequency

th     = -pi:2*pi/(100*max(m)):pi;

[~,i]  = max(real( eta'*exp(1i*m*th) ));
mu     = th(i);
d1     = -imag( eta' * ( m    .* exp(1i*m*mu) ) );
d2     = -real( eta' * ( m.^2 .* exp(1i*m*mu) ) );
if d2<0 % if the function is locally concave (usually the case)
    theta  = mu - d1/d2;
else    % if the function is not locally concave (not sure if ever the case)
    theta  = mu;
end
a      = exp(1i*(0:1:m(end))' * theta);
end

function [ K, s, w, C ] = maxZ_diag( J, h, M, sigma, rho, tau )
%maxZ maximizes the function Z of the binary vector s, see Appendix A of
%the paper

L = size(h,1);
cnst = log(rho/(1-rho)/tau);

K = 0; % number of components
s = false(L,1); % Initialize s
w = zeros(L,1);
C = zeros(L);
u = zeros(L,1);
v = zeros(L,1);
Delta = zeros(L,1);
if L > 1
    cont = 1;
    while cont
        if K<M-1
            v(~s) = 1 ./ ( sum(1./sigma) + 1/tau - real(sum(J(s,~s).*conj(C(s,s)*J(s,~s)),1)) ); % diagonal
            u(~s) = v(~s) .* ( h(~s) - J(s,~s)'*w(s));
            Delta(~s) = log(v(~s)) + u(~s).*conj(u(~s))./v(~s) + cnst;
        else
            Delta(~s) = -1; % dummy negative assignment to avoid any activation
        end
        if ~isempty(h(s))
            Delta(s) = -log(diag(C(s,s))) - w(s).*conj(w(s))./diag(C(s,s)) - cnst;
        end
        [~, k] = max(Delta);
        if Delta(k)>0
            if s(k)==0 % activate
                w(k) = u(k);
                ctemp = C(s,s)*J(s,k);
                w(s) = w(s) - ctemp*u(k);
                C(s,s) = C(s,s) + v(k)*(ctemp*ctemp');
                C(s,k) = -v(k)*ctemp;
                C(k,s) = C(s,k)';
                C(k,k) = v(k);
                s(k) = ~s(k); K = K+1;
            else % deactivate
                s(k) = ~s(k); K = K-1;
                w(s) = w(s) - C(s,k)*w(k)/C(k,k);
                C(s,s) = C(s,s) - C(s,k)*C(k,s)/C(k,k);
            end
            C = (C+C')/2; % ensure the diagonal is real
        else
            break
        end
    end
elseif L == 1
    if s == 0
        v = 1 ./ ( sum(1./sigma) + 1/tau );
        u = v * h;
        Delta = log(v) + u*conj(u)/v + cnst;
        if Delta>0
            w = u; C = v; s = 1; K = 1;
        end
    else
        Delta = -log(C) - w*conj(w)/C - cnst;
        if Delta>0
            w = 0; C = 0; s = 0; K = 0;
        end
    end
end
end

function [ k ] = Ainv( R )
% Returns the approximate solution of the equation R = A(k),
% where A(k) = I_1(k)/I_0(k) is the ration of modified Bessel functions of
% the first kind of first and zero order
% Uses the approximation from
%       Mardia & Jupp - Directional Statistics, Wiley 2000, pp. 85-86.
%
% When input R is a vector, the output is a vector containing the
% corresponding entries

k   = R; % define A with same dimensions
in1 = (R<.53); % indices of the entries < .53
in3 = (R>=.85);% indices of the entries >= .85
in2 = logical(1-in1-in3); % indices of the entries >=.53 and <.85
R1  = R(in1); % entries < .53
R2  = R(in2); % entries >=.53 and <.85
R3  = R(in3); % entries >= .85

% compute for the entries which are < .53
if ~isempty(R1)
    t      = R1.*R1;
    k(in1) = R1 .* ( 2 + t + 5/6*t.*t );
end
% compute for the entries which are >=.53 and <.85
if ~isempty(R2)
    k(in2) = -.4 + 1.39*R2 + 0.43./(1-R2);
end
% compute for the entries which are >= .85
if ~isempty(R3)
    k(in3) = 1./( R3.*(R3-1).*(R3-3) );
end

end

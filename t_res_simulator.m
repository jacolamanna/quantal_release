function [minis, Texp, T0] = t_res_simulator(Nmax, ll, alpha)
% code for simulating a release process based on fBm (vesicle motion) + exponential delay 1 (endocytosis) + exponential delay 2 (exocytosis)
% inputs -> Nmax: number of events; ll: [lambda1 lambda2]; alpha: Hurst index - 2
H = alpha + 2;
beta = 1 - H;
dT = 10^-3; % 0.01 ms - fbm sampling time
Nsim = 100000; 
no_fbm = 0;
 
% Variables
T0 = []; %tempi di primo ritorno
Texp = []; %componente esponenziale
Ttot = 0;
Nfus = 0;
Nout = 0;
 
while Nfus<Nmax
     
    % Vesicle flight
    n = 2^ceil(log2(Nsim));
    T = 3600;
     
    r = nan(n+1,1); 
    r(1) = 1; 
    idx = 1:n;
    r(idx+1) = 0.5*((idx+1).^(2*H) - 2*idx.^(2*H) + (idx-1).^(2*H));
    r =[r; r(end-1:-1:2)]; % first rwo of circulant matrix
    lambda = real(fft(r))/(2*n); % eigenvalues
    W = fft(sqrt(lambda).*complex(randn(2*n,1),randn(2*n,1)));
    W = n^(-H)*cumsum(real(W(1:n+1))); % rescale
    W = T^H*W; % t =(0:n)/n; % t = t*T;

    bm = [0;W];
     
    % Extraction of first passage time
    if bm(2)<0
        bm = -bm;
    end
    tt = find(bm <= 0, 2, 'first');
    if length(tt)>1
        Nfus = Nfus+1;
        t_fbm = tt(2) * dT;
        if no_fbm == 1
            t_fbm = 0;
        end
        t_exp = exprnd(1/ll(1))+exprnd(1/ll(2)); % endocytosis + exocytosis
        T0 = [T0 t_fbm];
        Texp = [Texp t_exp];
    else
        Nout = Nout + 1;
    end
end

Tfus = (Texp + T0);
minis = cumsum([0 Tfus]);







function [Ti, Ni_norm, alpha, ll, ll_full, C_full, alpha_full] = hyper_model2024(minis, fit, xmin, plotP, pars_in)
% Code for displaying interval distributions of minis (quantal release events) using log-binned histograms.
% The code includes histogram fitting using different models, from basic power-law decay to the diffusion-based model pdf (requires external functions)

% INPUTS
% minis: sequence of minis times | fit: type of fitting (see code) | xmin: minimum interval for decay fitting | plotP: = 1 produces a plot | pars_in: fixed models parameters (fit = 4)

% OUTPUTS
% Ti: histogram bins | Ni_norm: histogram normalized counts | alpha: fitted from decay | ll: fitted from biexponential | *_full: parameter fitted from the complete model

if nargin < 2
    fit = 0;
    alpha = NaN;
    ll = NaN(2,1);
    plotP = 1;
end
if nargin < 4
    plotP = 1;
end

alpha_full = NaN;


isi = diff(minis);

% remove outliers
isi = isi(isi>10^-4);

% Log Binning
dX = 0.1;
Xs = floor(min(log(isi)));
mlog = round(max(log(isi)))+1;
Xi = Xs:dX:mlog;
ii = 0:length(Xi)-1;
Ti = exp(Xs) .* exp(ii * dX);

Ni = zeros(size(Xi));
for i=1:length(Xi)
    Ni(i) = sum( log(isi) >= Xi(i) & log(isi) < Xi(i)+dX );
end

Ni_norm = Ni./[Ti(1) diff(Ti)] / sum(Ni);
Ni_norm_save = Ni_norm;
Ni_norm(Ni_norm<=10^-5) = 10^-5;

if plotP==1
    stairs(Ti, Ni_norm,'k');
    xlimi = get(gca,'XLim');
    ylimi = get(gca,'YLim');
end   

xlabel('Time [s]');
ylabel('Probability density [1/sec]');

if fit>0
    % Power-law decay fitting
    Ti_d = [0 diff(Ti)];
    Ti_c = Ti+Ti_d/2;
    
    xmax = find(Ni_norm > 10^-5,1,'last');
    Ti_fit_i = Ti_c>=xmin; % & Ti_c<=Ti_c(xmax);
    Ti_fit = Ti_c(Ti_fit_i);
    Ni_fit = Ni_norm(Ti_fit_i);

    Ni_fit_ok = Ni_fit > 10^-5;
    Ni_fit = Ni_fit(Ni_fit_ok);
    Ti_fit = Ti_fit(Ni_fit_ok);

    Ti_fitL = log(Ti_fit);
    Ni_fitL = log(Ni_fit);
    
    p1 = polyfit(Ti_fitL,Ni_fitL,1);
    alpha = p1(1);

    Ni_vL = polyval(p1,Ti_fitL);

    if plotP==1 && fit~=6 && fit~=7
        hold on, plot(exp(Ti_fitL),exp(Ni_vL),'b--','LineWidth',2);
    end

    % LME
    startpt = [1/mean(isi) 1/15];
    Fexp = @(x,l1,l2) (x>=0) .* ( (l1*l2/(l1-l2))*(exp(-l2*x)-exp(-l1*x)) ); 
    [ll,~] = mle(isi(isi<xmin),'pdf',Fexp,'Start',startpt, 'LowerBound',[0 0]);
    h_tot = 1;
    ll = sort(ll,2,"descend"); % So that lambda1 > lambda2 and u1 (release) < u2 (endocytosis = 15 s
    if plotP==1 && fit>=1.5 && fit~=4 && fit~=6 && fit~=7
        hold on, plot(Ti, hypo_exp(ll,Ti) * h_tot,'r--','LineWidth',2)
    end

    % Closed-form model (lama_model2024) with H = 0.5 -> alpha = -1.5 (standard BM)

    Dfix = 0.11; % um^2/s % D = 4 * 10^-11 * 10^-4; % m^2/s
    Cdiff = (2 * sqrt(Dfix));
    Cfix = Cdiff^2;
    Cfix2 = Cdiff;
    levy_l = Cfix;

    clear C l1 l2
    if fit == 2

        clear C l1 l2 x_val x
        initial_params = [1e-12, ll(1), ll(2), alpha];
       
        % MLE
        options = statset('MaxIter',10000,'MaxFunEvals',10000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');        
        pars = mle(isi, 'pdf', @(x,C,l1,l2) diffusion_pdf_conv2024(x,C,l1,l2,-1.5), 'Start', initial_params(1:end-1), 'LowerBound',[1e-12 1e-12 1e-12], 'UpperBound', [Inf Inf Inf], 'Options',options);
        pars(4) = -1.5;

        ll_full = pars(2:3);
        C_full = pars(1);
        D_full = (C_full / 2)^2;
        alpha_full = -1.5;

        if plotP==1
            hold on, plot(Ti_c, lama_model2024([ll_full C_full],Ti_c) * h_tot,'g-','LineWidth',2)
         end

    end

    if fit == 3 % Fit numerical convolution pdf (diffusion_pdf_conv2024.m) using fixed alpha at -1.5 (BM)

        initial_params = [1e-12, ll(1), ll(2), alpha];
        clear C l1 l2 x_val x 
       
        % MLE FOR STANDARD BM
        options = statset('MaxIter',10000,'MaxFunEvals',10000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');        
        pars = mle(isi, 'pdf', @(x,C,l1,l2) diffusion_pdf_conv2024(x,C,l1,l2,-1.5), 'Start', initial_params(1:end-1), 'LowerBound',[1e-12 1e-12 1e-12], 'UpperBound', [Inf Inf Inf], 'Options',options);
        pars(4) = -1.5;
        
        alpha_full = pars(4); % alpha here if fixed
        C_full = pars(1);
        ll_full = pars(2:3);
        % C = 2sqrt(D); -> C_full = 2sqrt(D) -> D = (Cfull / 2)^2
        D_full = (C_full / 2)^2;

        if plotP==1
            % cannot use conv() because not equally spaced Ti_c
            yy = diffusion_pdf2024(Ti_c, C_full,ll_full(1),ll_full(2), alpha_full) * h_tot;
            sumP = sum(Ni_norm);
            yy = sumP * yy / sum(yy);
            hold on, plot(Ti_c, yy,'g-','LineWidth',2);
        end

    end

    if fit == 5 % Fit pdf (diffusion_pdf2024.m) using alpha fitted from power-law (model approximation for generic alpha)

        initial_params = [0.1, ll(1), ll(2), alpha];
        clear C l1 l2 x_val x 

        % MLE FOR GENERAL ALPHA (H)
        options = statset('MaxIter',10000,'MaxFunEvals',10000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');        
        pars = mle(isi, 'pdf', @(x,C) diffusion_pdf2024(x,C,ll(1),ll(2),alpha), ...
            'Start', initial_params(1), 'LowerBound',[0], 'UpperBound', [Inf], 'Options',options);
        pars(4) = alpha;
        pars(2:3) = ll;

        alpha_full = pars(4); % alpha here if fixed
        C_full = pars(1);
        ll_full = pars(2:3);
        % C = 2sqrt(D); -> C_full = 2sqrt(D) -> D = (Cfull / 2)^2
        D_full = (C_full / 2)^2;

        if plotP==1
            % cannot use conv() because not equally spaced Ti_c
            yy = diffusion_pdf2024(Ti_c, C_full,ll_full(1),ll_full(2), alpha_full) * h_tot;
            sumP = sum(Ni_norm);
            yy = sumP * yy / sum(yy);
            hold on, plot(Ti_c, yy,'g-','LineWidth',2);
        end

    end

    if fit == 4
        yy = diffusion_pdf2024(Ti_c,pars_in(1),pars_in(2),pars_in(3),pars_in(4)) * h_tot;
        sumP = sum(Ni_norm);
        yy = sumP * yy / sum(yy);
        hold on, plot(Ti_c, yy,'g-','LineWidth',2);
    end

    if fit == 6 % Fit biexponential model only
               
        % MLE
        clear l1 l2 x
        pdf_biexp = @(x,l1,l2) exppdf(x,l1)+exppdf(x,l2);
        options = statset('MaxIter',10000,'MaxFunEvals',10000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');        
        pars = mle(isi, 'pdf', pdf_biexp, 'Start', [1 1], 'LowerBound',[0 0], 'UpperBound', [Inf Inf], 'Options',options);

      
        if plotP==1
            % cannot use conv() because not equally spaced Ti_c
            yy = pdf_biexp(Ti_c,pars(1),pars(2));
            sumP = sum(Ni_norm);
            yy = sumP * yy / sum(yy);
            hold on, plot(Ti_c, yy,'b-','LineWidth',2);
        end

    end

    if fit == 7 % Fit triexponential model only
               
        % MLE
        clear l1 l2 l3 x
        pdf_biexp = @(x,l1,l2,l3) exppdf(x,l1)+exppdf(x,l2)+exppdf(x,l3);
        options = statset('MaxIter',10000,'MaxFunEvals',10000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');        
        pars = mle(isi, 'pdf', pdf_biexp, 'Start', [2 2 2], 'LowerBound',[0 0 0], 'UpperBound', [Inf Inf Inf], 'Options',options)

      
        if plotP==1
            % cannot use conv() because not equally spaced Ti_c
            yy = pdf_biexp(Ti_c,pars(1),pars(2),pars(3));
            sumP = sum(Ni_norm);
            yy = sumP * yy / sum(yy);
            hold on, plot(Ti_c, yy,'g-','LineWidth',2);
        end

    end

    if fit == 1
        % Display estimations
        text(10, 10^-0.4, {['\alpha = {\color{blue}',num2str(alpha,3),'}']},'Interpreter','tex');
    elseif fit == 1.5
        % Display estimations
        text(10, 10^-0.4, {['\alpha = {\color{blue}',num2str(alpha,3),'}'], ['\lambda_1 = {\color{red}',num2str(ll(1),3),'} s^{-1}'], ...
            ['\lambda_2 = {\color{red}',num2str(ll(2),3),'} s^{-1}']},'Interpreter','tex');        
    elseif fit == 2 || fit == 3 || fit == 5
        text(10, 10^-0.4, {['\alpha = {\color{blue}',num2str(alpha,3),'} | {\color{green}',num2str(alpha_full,3),'}'], ['\lambda_1 = {\color{red}',num2str(ll(1),3),'} | {\color{green}',num2str(ll_full(1),3),'} s^{-1}'], ...
            ['\lambda_2 = {\color{red}',num2str(ll(2),3),'} | {\color{green}',num2str(ll_full(2),3),'} s^{-1}'], ...
            ['C = {\color{green}',num2str(C_full(1),2),'}'], ['D = {\color{green}',num2str(D_full(1),2),' m^{2}/s}']},'Interpreter','tex');
    elseif fit == 4
        D_full = (pars_in(1) / 2)^2;
        text(10, 10^-0.4, {['\alpha = {\color{blue}',num2str(alpha,3),'} | \alpha = {\color{green}',num2str(pars_in(4),3),'} '], ...
            ['\lambda_1 = {\color{green}',num2str(pars_in(2),3),'} s^{-1}'], ...
            ['\lambda_2 = {\color{green}',num2str(pars_in(3),3),'} s^{-1}'], ...
            ['C = {\color{green}',num2str(pars_in(1),2),'}'], ['D = {\color{green}',num2str(D_full(1),2),' m^{2}/s}']},'Interpreter','tex');
    end

    
end

if plotP==1
    set(gca,'XScale', 'Log', 'YScale', 'Log');
    set(gca,'XLim',[10^-2.5 10^3],'YLim',[1.1*10^-5 max(Ni_norm)*2])
end









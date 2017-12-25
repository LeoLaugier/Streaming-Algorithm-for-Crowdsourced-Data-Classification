% Plot average estimation error with respect to n 

clear all
close all

alpha = 1; % 0 <= alpha <= 1
m = 50;      % numbe of labellers, natural number >= 3
s = 50;      % number of tasks, natural number
% q= ones(3,1);
q = [0.25 , 0.3 , 0.4];
p_n = []; % represents p(s,i)_(1<=i<=n), for 3<=n<=m
p_hat_n = []; % represents p_hat(s,i)_(1<=i<=n), for 3<=n<=m
average_error=zeros(m/2 -1,1);

q_select = q(1) %%%%%%%%%%% choice of q
%for j = 1:3
    for p = 2:(m/2)
        n=2*p;
        % Ground truth (p_i)1<=i<=n
        %% Methode 1
        % while q(j) >= 1/2 -1/n
        %     p = rand(n,1);
        %     q(j) = 1/n * sum(p);
        % end
        
        %% Methode 2
        % p = 0.3*ones(1,n)';
        
        %% Methode 3
        p_n = [(2*q_select-1/2) * ones(1,n/2)  1/2 * ones(1,n/2) ];
        
        % Ground truth (G_t)1<=t<=s
        G = randi([0 1], s,1);
        G(G==0) = -1;
        
        % Generation of (x_(t,i))1<=t<=s,1<=i<=n
        x = zeros(s,n);
        for t = 1:s
            for i = 1:n
                r1 = rand;
                if r1  < alpha * p_n(i)
                    x(t,i) = -G(t);
                else if (r1 < alpha )
                        x(t,i) = G(t);
                    else
                        x(t,i) = 0;
                    end
                end
            end
        end
        
        %%% Estimators (p_hat_i)1<=i<=n and (G_hat_t)1<=t<=s
        step = 2; % pas intervenant dans la recherche de v(a)
        epsilon = 0.0001; % precision in the search of v(a)
        
        %% Calculus of (a_hat_(t,i))1<=t<=s,1<=i<=n
        a_hat = zeros(s,n);
        
        % Initialization
        a_hat(1,:)=( sum(x(1,:)) * x(1,:) + ( sum(abs(x(1,:))) - 2 ) * abs(x(1,:)) ) / (2*(n-1)*alpha^2);
        for t = 2:s
            a_hat(t,:) = (t-1)/t * a_hat(t-1,:) + ( sum(x(t,:)) * x(t,:) + ( sum(abs(x(t,:))) - 2 ) * abs(x(t,:)) ) / (2*t*(n-1)*alpha^2);
        end
        
        %% Calculus of (v(a_hat_(t,i)))1<=t<=s,1<=i<=n
        v_0 = zeros(s,1); % v_0 represents the vector v_0(a_hat)
        v = zeros(s,1); % v represents the vector v(a_hat) : the solution of f(a_hat,v) = v
        p_hat = zeros(s,n);
        w_hat = zeros(s,n);
        G_hat = zeros(s,1);
        for t = 1:s
            v_0(t) = max(4 * (n-1)/n^2 * max(2*a_hat(t,:)-1) , 0 );
            y = v_0(t);
            if f(a_hat(t,:),v_0(t),n) <= v_0(t)
                
                %affichage f
                %         abscisse = [v_0(t):0.01:v_0(t)+1];
                %         for abs_i=1:length(abscisse)
                %             image(abs_i) = f(a_hat(t,:),abscisse(abs_i),n) - abscisse(abs_i);
                %         end
                %         plot(abscisse, image);
                
                while f(a_hat(t,:),y,n) < y
                    f(a_hat(t,:),y,n);
                    y = y*step;
                end
                
                vg = v_0(t);
                vd = y;
                y = (vg + vd)/2;
                while abs(vd - vg) > epsilon
                    if f(a_hat(t,:),y,n) < y
                        vg = y;
                        y = (vg + vd)/2;
                    else if f(a_hat(t,:),y,n) > y
                            vd = y;
                            y = (vg + vd)/2;
                        else
                            vg = y;
                            vd = y;
                        end
                    end
                end
                
                v(t) = y;
                
                % Calculus of (p_hat_(t,i))1<=t<=s,1<=i<=n and (w_hat_(t,i))1<=t<=s,1<=i<=n
                p_hat(t,:) = g(a_hat(t,:),v(t),n);
            else
                p_hat(t,:) = 1/2 * ones(1,n)';
            end
            w_hat(t,:) = log( 1./p_hat(t,:) - 1);
            
            % Calculus of (G_hat(t))1<=t<=s
            z = w_hat(t,:) * x(t,:)';
            if z>0
                G_hat(t) = 1;
            else if z<0
                    G_hat(t) = -1;
                else
                    r2 = rand;
                    if r2 < 0.5
                        G_hat(t) = 1;
                    else
                        G_hat(t) = -1;
                    end
                end
            end
        end
        
        p_hat_n =  p_hat(s,:);
        average_error(p-1) = 1/n *sum(abs(p_hat_n - p_n));
        % abscisse = [1:s];
        %
        % for t=1:s
        %     average_error(t) = 1/n * sum(abs(p_hat(t,:)' - p));
        % end
        % figure (1)
        % plot(abscisse,  average_error, 'LineWidth', 2);
        % legend(['q = '  num2str(q(1))],['q = '  num2str(q(2))],['q = '  num2str(q(3))]);
        % hold on
        %
        % figure (2)
        % for t=1:s
        %     maximum_error(t) = max(abs(p_hat(t,:)' - p));
        % end
        %
        % plot(abscisse,  maximum_error, 'LineWidth', 2);
        % legend(['q = '  num2str(q(1))],['q = '  num2str(q(2))],['q = '  num2str(q(3))]);
        % hold on
        
    end
    
    figure(1)
    abscisse = [4:2:m];
    
    figure (1)
    plot(abscisse,  average_error, 'LineWidth', 2);
    legend(['q = '  num2str(q_select)]);
    hold on
%end

figure (1)
title('Average estimation error $\frac{1}{n} * ||\hat{p}(t) - p||_{1}$', 'interpreter','latex')
xlabel('Number of labellers n')
ylabel('Average estimation error')
 
% figure (2)
% title('Maximum estimation error $||\hat{p}(t) - p||_{\infty}$', 'interpreter','latex')
% xlabel('Number of tasks t')
% ylabel('Maximum estimation error')
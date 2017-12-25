% Plot estimation errors with respect to alpha 
clear all
close all

alpha = [0.28:0.08:1]; % 0 <= alpha <= 1
n = 50;      % numbe of labellers, natural number >= 3
s = 50;      % number of tasks, natural number
% q= ones(3,1);
q = [0.25 , 0.3 , 0.4];
p_alpha = ones(10,n); 
p_hat_alpha = ones(10,n); 
%for j = 1:3
    for j = 1:10

        % Ground truth (p_i)1<=i<=n
        %% Methode 1
        % while q(j) >= 1/2 -1/n
        %     p = rand(n,1);
        %     q(j) = 1/n * sum(p);
        % end
        
        %% Methode 2
        % p = 0.3*ones(1,n)';
        
        %% Methode 3

        p_alpha(j,:) = [(2*q(2)-1/2) * ones(1,n/2)  1/2 * ones(1,n/2) ];
        
        % Ground truth (G_t)1<=t<=s
        G = randi([0 1], s,1);
        G(G==0) = -1;
        
        % Generation of (x_(t,i))1<=t<=s,1<=i<=n
        x = zeros(s,n);
        for t = 1:s
            for i = 1:n
                r1 = rand;
                if r1  < alpha(j) * p_alpha(j,i)
                    x(t,i) = -G(t);
                else if (r1 < alpha(j) )
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
        a_hat(1,:)=( sum(x(1,:)) * x(1,:) + ( sum(abs(x(1,:))) - 2 ) * abs(x(1,:)) ) / (2*(n-1)*alpha(j)^2);
        for t = 2:s
            a_hat(t,:) = (t-1)/t * a_hat(t-1,:) + ( sum(x(t,:)) * x(t,:) + ( sum(abs(x(t,:))) - 2 ) * abs(x(t,:)) ) / (2*t*(n-1)*alpha(j)^2);
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
        
        p_hat_alpha(j,:) = p_hat(s,:);
        
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
    abscisse = alpha;
    average_error = zeros(j,1);
    for j=1:10
        average_error(j) = max((abs(p_hat_alpha(j,:) - p_alpha(j,:))));
    end
    figure (1)
    plot(abscisse,  average_error, 'LineWidth', 2);
    legend(['t = ' num2str(s)  ]);
    hold on
%end

figure (1)
title('Maximum estimation error $||\hat{p}(t) - p||_{\infty}$', 'interpreter','latex')
xlabel('Answer probability \alpha ')
ylabel('Average estimation error')
 
% figure (2)
% title('Maximum estimation error $||\hat{p}(t) - p||_{\infty}$', 'interpreter','latex')
% xlabel('Number of tasks t')
% ylabel('Maximum estimation error')
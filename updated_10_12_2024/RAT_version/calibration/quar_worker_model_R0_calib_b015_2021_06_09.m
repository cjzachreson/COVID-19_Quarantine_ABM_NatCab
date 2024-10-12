


clear all
close all

%rng(1)
%rng(2) % used for 2021 10 05 - first set
rng(3)
output_dirname = [pwd, '\output\Delta_RAT_calib_b015_2021_10_14'];
if ~isdir(output_dirname)
    mkdir(output_dirname)
end

Hellewell_data = readtable([pwd '\references\Hellewell_et_al_test_sensitivity_vs_time.csv']);

print_sensitivity_timeline_flag = true;
print_TOST_flag = true;

plot_sensitivity_timeline_flag = true;
plot_single_instances_flag = false;
print_single_instances_flag = false;

n_replicates = 1000;

N = 500;
%N = 5
%N = 4
%betas = 0.1:0.05:4

%betas = 1.5:0.1:2.5
%betas = 0.1:0.1:3

betas = 1.52

R0_beta = zeros(size(betas));

TOST_vs_beta = [];

% messing around with viral load

Vmax = 7;

% this is acutally reciprocal dispersion, I think...
dispersion = 0.15;


for b_i = 1:numel(betas)
    
   % clearvars -except n_replicates betas b_i R0_beta output_dirname N TOST_vs_beta Vmax t_plat dispersion
    
    TOST_b_i = NaN(n_replicates * 2, 1);
    
    %beta = 0.115; % R0 = 1.1; % force of infection per day, from infected individual
    
    beta = betas(b_i); % corresponding to the mean
    
    
    inc_mu = 1.62;
    
    % for mean of 4 days, use
    %inc_mu = 1.39;
    
    inc_sig = 0.418;
    
    symp_min = 5; %lower-bound of viral load wind-down period (i.e., symptomatic, contagious period)
    symp_max = 10; % upper-bound of viral load wind-down period
    
    
    inc_plat_fac = 0.1;
    % fraciton of incubation period for which infectiousness plateaus
    
    rec_plat_fac = 0;
    % fraction of recovery period for which infectiousness platau continues.
   
    
    asymp_fraction = 0.33; %fraction of agents that do not express symptoms
    
    %beta = 0.115; % R0 = 1.1
    
    %0.0023 * 50; % force of infection per day, for R0 =1.1
    % needs to be divided by workplace population***
    
    % calibrated for N = 50
    
    % to check R0, turn test sensitivity to 0, and asymp_fraction to 1.0
    % this will prevent termination upon detection or symptom expression.
    
    %tf = 100;
    dt = 0.1; % in days
    
    %%
    
    %n_replicates = 1000;
    
    R0_n = zeros(n_replicates, 1);
    
    TOSD = NaN(n_replicates, 1);
    
    infections_n = zeros(n_replicates, 1);
    
    
    timeline = -30:dt:60;
    t0_i = find(timeline == 0);
    
    sensitivity_timeseries_rel = NaN(n_replicates, numel(timeline));
    sensitivity_timeseries_absolute = NaN(n_replicates, numel(timeline));
    
    for r = 1:n_replicates
        
        %%
        r
        
        
        %output structures
        sz = 100;
        secondary_infections = zeros(N, ceil(sz/dt)+1);
        
        %% intialize agent properties
        infection_status = zeros(N, 1);
        symptomatic_status = zeros(N, 1);
        latent_status = zeros(N, 1);
        incubating_status = zeros(N, 1);
        recovering_status = zeros(N, 1);
        recovery_status = zeros(N, 1);
        infection_time = zeros(N, 1);
        
        k_inc = zeros(N, 1);
        k_r = zeros(N, 1);
        
        beta_i = zeros(N, 1);
        
        
        infection_clocks = zeros(N, 1);
        
        incubation_periods = zeros(N, 1);
        
        recovery_periods = zeros(N, 1);
        
        % initialize incubation periods, recovery periods and latent
        % periods
        
        % chop off the top of the viral load curve, call it a plateau
        
        %incubation period distribution
        inc_pd = makedist('Lognormal', 'mu', inc_mu, 'sigma', inc_sig);
        
        p_vals = [0:0.001:1];
        
        inc_inverse_cdist = inc_pd.icdf(p_vals);
        
        
        % parameters for piecewise sigmoid of test sensitivity:
        
        
        %         b1 = 1.5;
        %         b2 = 2.19;
        %         b3 = -1.1;
        
        % position of breakpoint relative to incubation period...
%         C_min = 2.01;
%         C_max = 5.11;

        C_min = 1;%2.01;
        C_max = 5.11;


        
        % maximum (value at breakpoint)
        % this is the part I need to change 2021-10-14 to reduce RAT
        % sensitivity
%         b1_min = 0.8;
%         b1_max = 2.31;
        b1_min = -0.282; % s = 0.43
        b1_max = 1.1; %s = 0.75
        
        % growth rate before breakpoint
        b2_min = 1.26;
        b2_max = 3.47;%3.47;
        %b2_min = 0.5; % allowing slower growth...
        % Note - this diverges from Hellewell's parameter estimates. NOTE:
        % produced strange results...
        
        
        % decay rate after breakpoing
%         b3_min = 1.05;
%         b3_max = 1.2;
        b3_min = 1.05;
        b3_max = 1.14;
        
        
        c_i = zeros(N, 1);
        b1_i = zeros(N, 1);
        b2_i = zeros(N, 1);
        b3_i = zeros(N, 1);
        
        index_case = ceil(rand() * N);
        
        
        
        for i = index_case
            
            
            while incubation_periods(i) < 1
                incubation_periods(i) = lognrnd(inc_mu, inc_sig);
            end
            
            % correlate our recovery period to incubation period via
            % quantile:
            
            [~, q_index] =  min(abs(inc_inverse_cdist - incubation_periods(i)));
            
            q = p_vals(q_index);
            
            % uncorrelated
            recovery_periods(i) = symp_min + rand() * (symp_max - symp_min);
            
            %correlated
            %recovery_periods(i) = symp_min + q * (symp_max - symp_min);
            
            
            %c_i(i) = C_min + q * (C_max - C_min); % longer incubation, longer growth period
            %c_i(i) = C_min + rand() * (C_max - C_min); % incubation, uncorrelated to growth period
            
            %c_i(i) = incubation_periods(i);% - min(q *(C_max - C_min), 2.2); % map C directly to incubation period - no cap.
            
            c_i(i) = incubation_periods(i) - q *(C_max - C_min); % map C to incubation period, with delay mapped to incubation
            
            %b1_i(i) = b1_min + (1 - q) * (b1_max - b1_min); % shorter incubation higher maximum detection probability
            
            b1_i(i) = b1_min + rand() * (b1_max - b1_min); % maximum detection probability not correlated to  incubation period
            
            %b1_i(i) = b1_min + (q) * (b1_max - b1_min); % longer incubation higher maximum detection probability
            
            b2_i(i) = b2_min + (1 - q) * (b2_max - b2_min); % shorter incubation period, faster growth
            
            %b3_i(i) = -1.1; % try holding this constant - there's already a correlation built in (faster growth, faster decay)
            
            %b3_i(i) = -0.25;
            
            %b3_i(i) = -1 *(b3_min + (1 - q) * (b3_max - b3_min)); % shorter incubation period, faster decay post-incubation
            
            %b3_i(i) = -1 *(b3_min + q * (b3_max - b3_min)); % shorter incubation period, longer decay post-incubation
            
            b3_i(i) = -1 *(b3_min + rand() * (b3_max - b3_min)); % incubation period, uncorrelated to decay post-incubation
            
            % note - correlation is built in to hellewell's framework -
            % the decay coefficient is a function of b2 and b3. if b3 is
            % fixed, decay will always occur as a fixed factor of b2.
            
            
            % these will be gamma distributed with mean beta, and the
            % specified dispersion parameter (shape parameter k )
            
            
            % set up timecourse of transmission potential:
            beta_dist_scale = beta/dispersion;
            
            beta_i(i) = gamrnd(dispersion, beta_dist_scale);
            
            t_plat_inc = ceil(incubation_periods(i) * inc_plat_fac / dt) * dt;
            t_plat_postinc = ceil(recovery_periods(i) * rec_plat_fac / dt) * dt;
            
            k_inc(i) = log(Vmax) / (incubation_periods(i) - t_plat_inc);
            k_r(i) = log(1/Vmax) / (recovery_periods(i) - t_plat_postinc);
            
            
            
            
            % for plotting infectiousness over time
            x1 = 0:dt:(incubation_periods(i) - t_plat_inc);
            x2 = x1(end) + dt : dt : x1(end) + t_plat_inc + t_plat_postinc;
            x3 = x2(end) + dt : dt : x2(end) + recovery_periods(i) - t_plat_postinc;
            
            y0 = beta_i(i) .* (1/Vmax);
            
            y1 = beta_i(i) .* (exp( k_inc(i) * x1 ) / Vmax);
            y2 = beta_i(i) .* ones(size(x2));
            y3 = beta_i(i) * exp( k_r(i) * (x3 - x2(end)) );
            
            x = [x1, x2, x3];
            y = [y1, y2, y3] - y0;
            
            
            %integral test
%             A1 = sum((y1 - y0)*dt)
%             
%             A1_test = (beta_i(i) / Vmax) * 1/k_inc(i) * (exp(k_inc(i) * x1(end)) - 1) - (y0 * x1(end))
%             
%             A2 = sum((y2 - y0) * dt)
%             
%             A2_test = beta_i(i) * (x2(end) - x2(1)) - (y0 * (x2(end)-x2(1)))
%             
%             A3 = sum((y3 - y0) * dt)
%             
%             A3_test = (beta_i(i)) * 1/k_r(i) * (exp(k_r(i) * (x3(end) - x3(1))) - 1) - (y0 * (x3(end) - x3(1)))
%             
%             A3_2 = sum((y3(1:50) - y0) * dt);
%             A3_2_test = (beta_i(i)) * 1/k_r(i) * (exp(k_r(i) * (x3(50) - x3(1))) - 1) - (y0 * (x3(50) - x3(1)))
%             
%             A3_diff = A3 - A3_2
%             A3_diff_Test = A3_test - A3_2_test
            
            % % set up timecourse of test sensitivity
            
            
            C_i = ceil(c_i(i)/dt)*dt;
            %
            t1 = 0:dt:C_i;
            %t2 = C_i + dt : dt : C_i + ceil(recovery_periods(i)/dt)*dt ;
            t2 = C_i + dt : dt : 50;%C_i + ceil(recovery_periods(i)/dt)*dt ;
            
            z1 = t1 - C_i;
            z2 = t2 - C_i;
            
            p1 = 1 ./ (1 + exp(-1.*(b1_i(i) + b2_i(i).*z1)));
            p2 = 1 ./ (1 + exp(-1.*(b1_i(i) + b2_i(i).*z2 + b2_i(i).*b3_i(i).*z2)));
            
            if plot_single_instances_flag
            
                Test_Sensitivity_i = [p1, p2];
                
                figure(3)
                yyaxis left
                plot(x, y./max(y) )
                ylabel('infectiousness')
                
                yyaxis right
                
                plot([t1, t2], Test_Sensitivity_i, '-', 'Color', [1, 0, 0, 0.5]);
                
                ylabel('test sensitivity')
                
                t_FoI = x';
                FoI_norm_i = [y./max(y)]';
                
                t_Sens = [t1, t2]';
                Test_sens_i = Test_Sensitivity_i';
                
                
                if print_single_instances_flag
                    
                    instance_table_FoI = table(t_FoI, FoI_norm_i);
                    instance_table_TSens = table(t_Sens, Test_sens_i);
                    
                    if incubation_periods(i) < 3.5
                        writetable(instance_table_FoI, [output_dirname, '\FoI_single_instance_example_short.csv']);
                        writetable(instance_table_TSens, [output_dirname, '\TSens_single_instance_example_short.csv'])
                        
                    else
                        if incubation_periods(i) > 5 && incubation_periods(i) < 6
                            
                            writetable(instance_table_FoI, [output_dirname, '\FoI_single_instance_example_avg.csv']);
                            writetable(instance_table_TSens, [output_dirname, '\TSens_single_instance_example_avg.csv'])
                        else
                            if incubation_periods(i) > 9 && incubation_periods(i) < 11
                                writetable(instance_table_FoI, [output_dirname, '\FoI_single_instance_example_long.csv']);
                                writetable(instance_table_TSens, [output_dirname, '\TSens_single_instance_example_long.csv'])
                            end
                        end
                    end
                end
                %pause
                
            end
            
            
            delay = ceil((incubation_periods(i) - c_i(i)) / dt);
            %
            if print_sensitivity_timeline_flag || plot_sensitivity_timeline_flag
                
                sensitivity_timeseries_rel(r, t0_i-delay-numel(p1)+1 : t0_i-delay) = p1;
                sensitivity_timeseries_rel(r, t0_i-delay +1 : t0_i-delay+numel(p2)) = p2;
                sensitivity_timeseries_absolute(r, 1:numel([p1, p2])) = [p1, p2];
                %
                %
                med_ts_rel = nanmedian(sensitivity_timeseries_rel);
                q90_top_rel = quantile(sensitivity_timeseries_rel, 0.95);
                q90_bottom_rel = quantile(sensitivity_timeseries_rel, 1-0.95);
                
                med_ts_ab = nanmedian(sensitivity_timeseries_absolute);
                q90_top_ab = quantile(sensitivity_timeseries_absolute, 0.95);
                q90_bottom_ab = quantile(sensitivity_timeseries_absolute, 1-0.95);
            end
            
            if plot_sensitivity_timeline_flag
                
                figure(6)
                plot([0:dt:numel(med_ts_ab)*dt-dt], med_ts_ab, 'k-')
                hold on
                plot([0:dt:numel(med_ts_ab)*dt-dt], q90_top_ab, 'r-')
                plot([0:dt:numel(med_ts_ab)*dt-dt], q90_bottom_ab, 'r-')
                %plot(Hellewell_data.days_since_infection, Hellewell_data.median, 'k-')
                plot(Hellewell_data.days_since_infection, Hellewell_data.lower_95, 'k--')
                plot(Hellewell_data.days_since_infection, Hellewell_data.upper_95, 'k--')
                hold off
                
                figure(7)
                plot(timeline, med_ts_rel, 'k-')
                hold on
                plot(timeline, q90_top_rel, 'r-')
                plot(timeline, q90_bottom_rel, 'r-')
                hold off
                
            end
            
            
            
        end
        
        latent_periods = dt;%incubation_periods * latent_period_median;
        
        %%
        %infect a random agent
        %         index_case = ceil(rand() * N);
        
        infection_status(index_case) = 1;
        latent_status(index_case) = 1;
        incubating_status(index_case) = 1;
        symptomatic_status(index_case) = (rand() > asymp_fraction);
        
        %%
        %iterate through time
        iteration = 0;
        
        TOSD_r = NaN;
        
        tf = incubation_periods(index_case) + recovery_periods(index_case) + dt;
        
        for t = 0:dt:tf
            iteration = iteration + 1;
            
            % iterate infection clocks
            infection_clocks = infection_clocks + dt * infection_status;
            
            %anyone over the sum of their incubation and recovery periods recovers
            recovery_status(infection_clocks > (incubation_periods + recovery_periods) & recovery_status == 0) = t;
            infection_status(infection_clocks > (incubation_periods + recovery_periods)) = 0;
            
            
            % transition from latent
            latent_status(infection_clocks > latent_periods) = 0;
            
            % transition from incubating to recovering
            incubating_status(infection_clocks > (incubation_periods)) = 0;
            recovering_status(infection_clocks > (incubation_periods) &...
                infection_clocks < (incubation_periods + recovery_periods)) = 1;
            recovering_status(infection_clocks > ( incubation_periods + recovery_periods)) = 0;
            
            
            % iterate through agents, and infect people.
            
            
            for i = index_case
                
                
                if infection_status(i) && ~latent_status(i)
                    
                    % compute detection
                    if mod(t, 1) == 0
                        p_detect = 0;
                        
                        %                         tau = infection_clocks(i) - incubation_periods(i);
                        
                        tau = infection_clocks(i) - c_i(i);
                        
                        %tau = infection_clocks(i) - (incubation_periods(i) - min(incubation_periods(i)*inc_peak_sens_fac, (5.5 - 3.18)));
                        
                        %tau = infection_clocks(i) - 3.19;
                        
                        if tau < 0%incubating_status(i)
                            
                            p_detect = 1 ./ (1 + exp(-1.*(b1_i(i) + b2_i(i).*tau)));
                            
                        else
                            
                            p_detect = 1 ./ (1 + exp(-1.*(b1_i(i) + b2_i(i).*tau + b2_i(i).*b3_i(i).*tau)));
                            
                        end
                        
                        if isnan(TOSD_r)
                            
                            if rand() < p_detect
                                TOSD_r = t - incubation_periods(i);
                            end
                        end
                        
                    end
                    
                    % compute pairwise infections
                    for j = 1:N
                        
                        if  ~infection_status(j) && ~recovery_status(j)
                            
                            beta0_i = beta_i(i) .* (1/Vmax);
                            
                            if incubating_status(i)
                                
                                t_inc =  infection_clocks(i);
                                
                                beta_it = beta_i(i) * (exp( k_inc(i) * t_inc ) / Vmax);
                                
                                
                            else
                                %beta_i = beta * (1 - (infection_clocks(i) - incubation_periods(i))/recovery_periods(i));
                                
                                tr = infection_clocks(i) - incubation_periods(i) ;
                                
                                beta_it = beta_i(i) *  exp( k_r(i) * tr );
                                
                            end
                            
                            if beta_it > beta_i(i) %plateau
                                beta_it = beta_i(i);
                            end
                            
                            beta_it = beta_it - beta0_i;
                            
                            pI = 1 - exp(-dt * beta_it / (N - 1));
                            
                            if pI > 1
                                'error'
                            end
                            
                            infect = rand() < pI;%((beta_it / (N-1)) * dt); the commented code is from when beta was constrained.
                            
                            
                            if infect
                                infection_status(j) = 1;
                                infection_time(j) = t;
                                latent_status(j) = 1;
                                incubating_status(j) = 1;
                                symptomatic_status(j) = rand() > asymp_fraction;
                                secondary_infections(i, iteration) = secondary_infections(i, iteration) + 1;
                                TOST_b_i(sum(~isnan(TOST_b_i)) + 1) = infection_clocks(i) - incubation_periods(i);
                            end
                        end
                        
                    end
                end
                
            end
            
            % for R0 calibration
            if infection_status(index_case) == 0 %|| (recovering_status(index_case) == 1 && symptomatic_status(index_case) == 1)
                break
            end
            
            
        end
        
        R_0 = sum(secondary_infections(index_case, :));
        
        R0_n(r) = R_0;
        
        TOSD(r) = TOSD_r ;
        
        
    end
    
    
    
    h = histogram(TOST_b_i(~isnan(TOST_b_i)), 'normalization', 'probability');
    t_TOST = [h.BinEdges(1:end-1) - h.BinWidth/2]';
    p_dens_TOST = [h.Values ./ h.BinWidth]';
    figure(1);
    plot(t_TOST, p_dens_TOST)
    title(['TOST, beta = ' num2str(betas(b_i))])
    
    %TOST fit from Feretti's preprint:
    
    mu = -0.0747;
    sigma = 1.8567
    nu = 3.3454
    
    TOST_Feretti = makedist('tLocationScale', 'mu', mu, 'sigma', sigma, 'nu', nu);
    x = -15:0.01:15;
    y = TOST_Feretti.pdf(x);
    
    hold on
    
    plot(x, y, 'r')
    hold off
    
    
    R0_beta(b_i) = nanmean(R0_n);
    
    figure(2)
    h2=histogram(R0_n, 'normalization', 'probability');
    xhist = h2.BinEdges(1:end-1) - h2.BinWidth/2;
    yhis = h2.Values;
    
    %fit to negative binomial to estimate dispersion
    nbin_pd = fitdist(R0_n, 'NegativeBinomial');
    x = 0:max(R0_n);
    y = pdf(nbin_pd, x);
    hold on
    plot(x, y)
    
    title(['Secondary case distribution, dispersion = ' num2str(nbin_pd.R)])
    
    hold off
    
    figure(4);
    h = histogram(TOSD, 'normalization', 'probability', 'BinWidth', 0.1);
    t_TOSD = [h.BinEdges(1:end-1) - h.BinWidth/2]';
    p_dens_TOSD = [h.Values ./ h.BinWidth]';
    p_cum_TOSD = [cumsum(h.Values)]';
    yyaxis left
    plot(t_TOSD, p_dens_TOSD, 'o')
    title(['TOSD'])
    yyaxis right
    plot(t_TOSD, p_cum_TOSD)
    
    
    
    disp(['*****'])
    disp(['beta = ' num2str(beta) ' ; R_0 = ' num2str(nanmean(R0_n)), ' ; std. = ' num2str(nanstd(R0_n))])
    disp(['*****'])
    
    if print_TOST_flag
        
        output_filename_1 = ['\secondary_cases_beta_' strrep(num2str(beta), '.', '-') '.csv'];
        dlmwrite([output_dirname, output_filename_1], R0_n);
        
        output_filename_2 = ['\TOST_beta_' strrep(num2str(beta), '.', '-') '.csv'];
        writetable(table(t_TOST, p_dens_TOST), [output_dirname, output_filename_2]);
        
        output_filename_3 = ['\TOSD_beta_' strrep(num2str(beta), '.', '-') '.csv'];
        writetable(table(t_TOSD, p_dens_TOSD, p_cum_TOSD), [output_dirname, output_filename_3]);
        
    end
    
    if print_sensitivity_timeline_flag
        
        t_rel =  [timeline]';
        
        t_abs = [0:dt:numel(med_ts_ab)*dt-dt]';
        
        med_ts_rel = med_ts_rel';
        q90_top_rel = q90_top_rel';
        q90_bottom_rel = q90_bottom_rel';
        
        med_ts_ab = med_ts_ab';
        q90_top_ab = q90_top_ab';
        q90_bottom_ab = q90_bottom_ab';
        
        sens_rel_incubation = table(t_rel, med_ts_rel, q90_top_rel, q90_bottom_rel);
        sens_rel_infection = table(t_abs, med_ts_ab, q90_top_ab, q90_bottom_ab);
        
        writetable(sens_rel_incubation, [output_dirname, '\Test_sensitivity_rel_symptoms_N_' num2str(N) '_reps_' num2str(n_replicates) '.csv']);
        writetable(sens_rel_infection, [output_dirname, '\Test_sensitivity_rel_infection_N_' num2str(N) '_reps_' num2str(n_replicates) '.csv']);
        
    end
    
    
end


beta = betas';
R0 = R0_beta';

SAR = R0 ./ (N - 1);

output_table = table(beta, R0, SAR);

output_filename = [output_dirname, '\R0_vs_beta_N_' num2str(N) '_reps_' num2str(n_replicates) '.csv'];

writetable(output_table, output_filename);

figure()
plot(betas, R0)
figure()
plot(betas, SAR.*100)






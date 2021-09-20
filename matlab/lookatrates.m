function lookatrates(varargin)
% Calculate best estimate and uncertainty of the element-specific (k) and
% cell-specific (r) rate of element assimilation by cells based on their
% isotopic composition measured by nanoSIMS.
%
% Theory behind the calculation method is explained in Polerecky et al. (submitted to Frontiers in Microbiology).
% Calculation approach is explained in https://github.com/lpolerecky/LARS.
%
% USAGE:
% lookatrates(Input_file, Nsimul, pause_for_each_cell, export_graphs_as_png);
%
% INPUT PARAMETERS:
% Input_file = name of the xlsx input file (assumed to be in the data subfolder)
% Nsimul = number of Monte-Carlo simulations for each cell (2000 or more)
% pause_for_each_cell = pause after calculating rate for each cell (0/1)
% export_graphs_as_png = export calculation results as PNG (0/1)
% Use [] if you want to use a default value for the argument.
%
% EXAMPLES:
% lookatrates; % all input parameters will have default values
% lookatrates('DataCells1.xlsx', 5000, 1, 1);
% lookatrates('DataCells1.xlsx', 2000, [], 0);
% lookatrates('DataCells1.xlsx');
%
% Written by Lubos Polerecky, 2020-06-07, Utrecht University

% default function parameters (see below for explanation)
default_parameters={'DataCells1.xlsx', 2000, 1, 1};

% fill in parameters specified by the user
for i=1:4
    fmt='WARNING: %s not specified. Using a default value of %d.\n';
    switch i
        case 1, varname='Input_file'; fmt='WARNING: %s not specified. Using a default value of %s.\n';
        case 2, varname='Nsimul';
        case 3, varname='pause_for_each_cell';
        case 4, varname='export_graphs_as_png';
    end
    if length(varargin)>=i
        if ~isempty(varargin{i})
            default_parameters{i} = varargin{i};
        else
            fprintf(1,fmt, varname, default_parameters{i});
        end
    else
        fprintf(1,fmt, varname, default_parameters{i});
    end
end

%% update parameters and flags that determine the actions below

% name of the input file (assumed to be in the data subfolder)
input_xlsx_filename=default_parameters{1};
input_data_folder = 'data'; % sub-folder of the main matlab file

% number of simulations per cell to determine kC and its error
Nsimul = default_parameters{2};

% make a pause to be able to check the results for each cell (useful when
% you want to ponder about the result for each individual cell, but not
% when you just want to get the results) 
pause_for_each_cell = default_parameters{3};

% export results for each cell as a png file
export_graphs_as_png = default_parameters{4};

% number of the figure where the output will be displayed
fign = 1;

% if you really want to see *every individual model prediction*, set this
% value to 1 (useful for *detailed* debugging, but not when you just want
% to get the results)
display_individual_results = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's go! Please do not modify below unless you know what you are doing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load experimental data
% It is up to the user to ensure that the formatting of the xlsx file is
% the same as in the default input file.
in_file = [input_data_folder filesep input_xlsx_filename];
fprintf(1,'Loading data from %s\n',in_file);
Tin = readtable(in_file)
fprintf(1,'Done\n');

% generate output filename
[~, b, c] = fileparts(input_xlsx_filename);
output_xlsx_filename = [input_data_folder filesep b '-' datestr(now,'HH-MM-SS') c];

% change values in avgVcell and Vcell that are not numbers (because the
% values cannot be constrained by experimental data) to NaN
if 0
t7=table2array(Tin(:,7));
t7num = NaN(size(t7));
if iscell(t7)
    for ii=1:length(t7)
        tmp = str2num(t7{ii});
        if ~isempty(tmp)
            t7num(ii) = tmp;
        end
    end
else
    t7num = t7;
end
% convert data to an array
exp_data = [table2array(Tin(:,1:6)) t7num table2array(Tin(:,8:9))];
comments = table2cell(Tin(:,10));
end

exp_data = table2array(Tin(:,1:9));
comments = table2cell(Tin(:,10));


% allocate matrix for the output
T = zeros(size(exp_data,1), 19);

if ~isfolder('png'), mkdir('png'); fprintf(1,'Subfolder png created.\n'); end

%% process input line by line
for ind=1:size(exp_data,1)
    
    %% extract data from the input table
    % note default formatting of columns
    t_incubation = exp_data(ind,1);
    xt = exp_data(ind,2);
    dxt = exp_data(ind,3);
    xS = exp_data(ind,4);
    xini = exp_data(ind,5);
    rhoC = exp_data(ind,6);
    avgVcell = exp_data(ind,7);
    Vcell = exp_data(ind,8);
    dVcell = exp_data(ind,9);
    
    % C content of the average and measured cell
    Ccell = Vcell * rhoC;
    dCcell = dVcell * rhoC;
    
    avgCcell = avgVcell * rhoC;
    
    if ~isnan(avgCcell) && ~isnan(Ccell)
        
        % max C content, i.e., when the cell is dividing 
        Cmax = avgCcell*2*log(2); % zero-order kinetics model
        %Cmax1 = avgCcell/log(2); % first-order kinetics model

        % ensure that the cell-cycle of the measured cell is between 0 and 1,
        % or set Ccell to avgCcell if this is not the case
        s_end = Ccell/(Cmax/2)-1;
        ds_end = dCcell/(Cmax/2);
        if (s_end>=1 || s_end<0) && abs(ds_end)>3*eps
            % values outside this range are "suspicious" 
            % set Ccell to avgCcell by default
            Ccell=avgCcell;
            s_end = Ccell/(Cmax/2)-1;
            fprintf(1,'WARNING: Cell biovolume, Vcell, for cell %d is not between Vmax/2 and Vmax.\nThis is not allowed, and so the value was set to avgVcell.\n',ind);
        end

        %% generate randomly distributed values of x(13C) and s_end

        % x(13C) corresponds to the 13C atom fraction of the cell
        % x(13C) is assumed to be *normally* distributed (mean=xt, SD=dxt)
        xt_rnd = mvnrnd(xt, dxt^2, Nsimul); 

        % determine Crand for a dividing cell
        if abs(ds_end)<3*eps
            Crand = Ccell*ones(Nsimul,1);
        else
            Crand = Cdistrib_partsync(Ccell,dCcell,Cmax,Nsimul);
        end
        
        % also determine Crand for a non-dividing cell
        Crand_nondiv = Cdistrib_nondiv(Ccell,dCcell,Nsimul);

        
        %% test values used during debugging
        %s_end_rnd=[0.5 0.5 0.5]';
        %xt_rnd=[0.08 0.05 0.02]';

        %% STEP 1: convert x(13C) to x(13C)SE
        % xSEt = substrate-normalized excess 13C atom fraction of the cell
        xSEt_rnd = (xt_rnd-xini)/(xS-xini);
        xSEt     = (xt-xini)/(xS-xini);
        dxSEt    = dxt/(xS-xini);
                
        if display_individual_results
            fprintf(1,'Finding best estimates of kC and s_ini for cell %d\n',ind)
        end

        %% find estimates of rC, kC and s_ini for each combination of xSEt, Crand and s_end

        % first allocate output matrices
        rC_A = zeros(size(xSEt_rnd));
        rC_B = zeros(size(xSEt_rnd));
        rC_C = zeros(size(xSEt_rnd));
        rC_nondiv = zeros(size(xSEt_rnd));
        kC_avg = zeros(size(xSEt_rnd));
        s0_ini_rnd = zeros(size(xSEt_rnd));
        s1_ini_rnd = zeros(size(xSEt_rnd));
        ndiv0 = zeros(size(xSEt_rnd));
        ndiv1 = zeros(size(xSEt_rnd));
        Ca2Ci = zeros(size(xSEt_rnd));
        Ca2Cf = zeros(size(xSEt_rnd));
        fac   = zeros(size(xSEt_rnd));
        Nsimul2 = length(Crand);
        x13Ct100_0 = nan(100, Nsimul2);
        St100_0 = nan(100, Nsimul2);    
        x13Ct100_1 = nan(100, Nsimul2);
        St100_1 = nan(100, Nsimul2);    

        for j=1:Nsimul2

            if mod(j,1000)==0
                fprintf(1,'j=%d/%d\n',j,Nsimul2);
            end

            % choose specific values
            xSEtj = xSEt_rnd(j);
            C_endj = Crand(j);

            % prepare function for finding zero (based on help to fzero)
            myfun = @(rC, x, s, t) dx(rC, x, s, t, avgCcell);
            x = xSEtj;
            s = C_endj;
            t = t_incubation;    
            fun = @(rC) myfun(rC, x, s, t);

            if xSEtj>0 && xSEtj<1

                %% STEP 2: convert xSE to k
                kC_rnd = -1/t_incubation * log(1-xSEtj);
                kC_avg(j) = kC_rnd;
                
                %% STEP 3: estimate r
                
                %% approach A
                rC_A(j) = kC_rnd * avgCcell; % use avg C content over cell cycle
                
                %% approach B
                rC_B(j) = kC_rnd * Crand_nondiv(j); % use C content of the measured cell
                
                %% approach C
                rC_C(j) = fzero(fun, rC_A(j));
                
                %% calculate also rC for a non-dividing cell (approach often used in lit.)
                rC_nondiv(j) = xSEtj/t_incubation * Crand_nondiv(j);

                %% reconstruct cell cycle (St) and number of cell divisions during SIP incubation
                
                % zero-order kinetics
                [~, x13Ct0, t0, ~, St0] = x13C_time0(rC_C(j), C_endj, avgCcell, t_incubation, xSEtj);
                s0_ini_rnd(j) = St0(1);
                ndiv0(j) = size(St0,2)-1;
                
                % first-order kinetics
                [~, x13Ct1, t1, ~, St1] = x13C_time1(rC_A(j), C_endj, avgCcell, t_incubation);
                s1_ini_rnd(j) = St1(1);
                ndiv1(j) = size(St1,2)-1;
                
                % relate the amount of assimilated C to the initial and
                % final C content of the cell
                Ca2Ci(j) = rC_C(j) * t_incubation / ( (St0(1)+1)*Cmax/2 );
                Ca2Cf(j) = rC_C(j) * t_incubation / ( (St0(end)+1)*Cmax/2 );
                
                % relate rC_C to rC_nondiv
                fac(j) = rC_C(j) / rC_nondiv(j);
                                   
                % store x13Ct and St interpolated in 100 data points
                % these values will be used later for plotting histograms
                [t100, x13Ct100, St100] = interpol_CtSt(t0, x13Ct0, St0);
                x13Ct100_0(:,j) = x13Ct100;
                St100_0(:,j) = St100;
                [~, x13Ct100, St100] = interpol_CtSt(t1, x13Ct1, St1);
                x13Ct100_1(:,j) = x13Ct100;
                St100_1(:,j) = St100;
                                
                %% display the results of the model prediction, if requested
                if display_individual_results
                    plot_individual_result(fign, t0, x13Ct0, St0, t1, x13Ct1, St1, t_incubation, xSEtj, dxSEt, j, Nsimul);
                end

            else
                fprintf(1,'WARNING: xt cannot be negative or above xS. Data-point skipped.\n');
                rC_A(j)         = NaN;
                rC_B(j)         = NaN;
                rC_C(j)         = NaN;                
                rC_nondiv(j)    = NaN;
                kC_avg(j)       = NaN;
                ndiv0(j)        = NaN;
                ndiv1(j)        = NaN;
                s0_ini_rnd(j)   = NaN;
                s1_ini_rnd(j)   = NaN;
                Ca2Ci(j)        = NaN;
                Ca2Cf(j)        = NaN;
                fac(j)          = NaN;
            end

        end    

        %% display results for the current cell
        if export_graphs_as_png || pause_for_each_cell
            plot_results_for_current_cell(fign, Nsimul2, ...
                t100, St100_0, x13Ct100_0, St100_1, x13Ct100_1, ...
                xt_rnd, s0_ini_rnd, s1_ini_rnd, ...
                rC_A, rC_B, rC_C, rC_nondiv, ...
                avgCcell, Crand, Crand_nondiv, ...
                t_incubation, xSEt, dxSEt);
        end

        %% store results for the current cell in the output array
        %% average over the cell cycle
        T(ind,1) = mean(kC_avg,'omitnan');
        T(ind,2) =  std(kC_avg,'omitnan');
        T(ind,3) = mean(rC_A,'omitnan');
        T(ind,4) =  std(rC_A,'omitnan');
        T(ind,5) = T(ind,4)/T(ind,3);
        %% instantaneous value at the end of SIP incubation, dividing cell
        T(ind,6) = mean(rC_B,'omitnan');
        T(ind,7) =  std(rC_B,'omitnan');
        T(ind,8) = T(ind,7)/T(ind,6);
        %% average values during the SIP incubation, dividing cell
        T(ind,9)  = mean(rC_C,'omitnan');
        T(ind,10) =  std(rC_C,'omitnan');
        T(ind,11) = T(ind,10)/T(ind,9);
        %% average values over the SIP incubation, non-dividing cell
        T(ind,12) = mean(rC_nondiv,'omitnan');
        T(ind,13) =  std(rC_nondiv,'omitnan');
        T(ind,14) = T(ind,13)/T(ind,12);
        %% avg number of cell divisions during SIP incubation
        T(ind,15) = mean(ndiv0,'omitnan');
        T(ind,16) = mean(ndiv1,'omitnan');
        %% important ratios
        T(ind,17) = mean(fac,'omitnan');        
        T(ind,18) = mean(Ca2Ci,'omitnan');
        T(ind,19) = mean(Ca2Cf,'omitnan');
        
        %% make a pause, if requested
        if pause_for_each_cell
            fprintf(1,'Summary of results for cell %d displayed in figure %d.\nComment: %s\n',ind,fign,comments{ind});
            input('Press enter to continue (or Ctrl+c to break).');
        else
            pause(0.1);
        end

        %% export the graphs for the current cell as PNG
        if export_graphs_as_png
           [~,b,~] = fileparts(input_xlsx_filename);
           outpng = sprintf('png%c%s_%03d.png', filesep, b, ind);
           print(fign,outpng, '-dpng','-r150');
           fprintf(1,'Graphs exported in %s\n\n',outpng);
        end
    
    else
        
        % estimate kC by step 2;
        xSEt        = (xt-xini)/(xS-xini);
        dxSEt       = dxt/(xS-xini);
        kC          = -1/t_incubation * log(1-xSEt);
        dkC         = 1/t_incubation * dxSEt / (1-xSEt);
        
        % if avgCcell is available, estimate r_A
        if ~isnan(avgCcell)
            rC_A      = kC * avgCcell;
            drC_A     = dxSEt/xSEt * rC_A; % error propagation formula
            rC_B      = NaN;
            drC_B     = NaN;
            rC_nondiv  = NaN;
            drC_nondiv = NaN;
        end
        
        % if Ccell is available, estimate r_B and r_nondiv
        if ~isnan(Ccell)
            rC_A    = NaN;
            drC_A   = NaN;
            rC_B    = kC * Ccell;
            drC_B   = sqrt((dxSEt/xSEt)^2 + (dCcell/Ccell)^2) * rC_B; % error propagation formula
            rC_nondiv  = xSEt/t_incubation * Ccell;
            drC_nondiv = sqrt((dxSEt/xSEt)^2 + (dCcell/Ccell)^2) * rC_nondiv; % error propagation formula            
        end
           
        % if Ccell is available, estimate r_B and r_nondiv
        if isnan(Ccell) && isnan(avgCcell)
            rC_A    = NaN;
            drC_A   = NaN;
            rC_B    = NaN;
            drC_B   = NaN;
            rC_nondiv  = NaN;
            drC_nondiv = NaN;
        end
        
        %% average across the cell cycle
        T(ind,1) = kC;
        T(ind,2) = dkC;
        T(ind,3:4) = [rC_A, drC_A];
        T(ind,5) = T(ind,4)/T(ind,3);        
        %% instantaneous value at the end of SIP incubation, dividing cell, zero-order kinetics
        T(ind,6) = rC_B;
        T(ind,7) = drC_B;
        T(ind,8) = T(ind,7)/T(ind,6);
        %% average values during the SIP incubation, dividing cell, first-order kinetics
        T(ind,9:11)  = NaN;
        %% average values over the SIP incubation, non-dividing cell, zero-order kinetics
        T(ind,12) = rC_nondiv;
        T(ind,13) = drC_nondiv;
        T(ind,14) = T(ind,13)/T(ind,12);
        %% avg number of cell divisions during SIP incubation
        T(ind,15:16) = 0; % default
        %% important ratios
        T(ind,17:19) = NaN;
        
    end

end

%% convert arrays to nice tables
T = array2table(T, 'VariableNames', ...
    {'k', 'dk', 'r_A', 'dr_A', 'CV_A', ... % average over cell-cycle
     'r_B', 'dr_B', 'CV_B', ... % instantaneous, end of SIP experiment
     'r_C', 'dr_C', 'CV_C', ...% average over SIP, dividing cell
     'r_nondiv', 'dr_nondiv', 'CV_nondiv', ... % average over SIP, non-dividing cell
     'Ndiv0', 'Ndiv1', ... % number of cell divisions during SIP
     'ratio', 'Ea/Ei', 'Ea/Ef'}); % important ratios

%% display results in the console window
fprintf(1,'EXPERIMENTAL VALUES:\n');
Tin
fprintf(1,'PREDICTED VALUES:\n')
T
pause(0.1);

%% export results in an xls file
% first, copy the content of the input file to the output file, then add
% the rates as a new sheet
copyfile(in_file,output_xlsx_filename,'f');
writetable(T,output_xlsx_filename, 'Sheet', 'rates');
fprintf(1,'Rates exported in %s (sheet rates).\n',output_xlsx_filename);

end % end of lookatrates function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% auxillary functions used above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [t100, x13Ct100, St100] = interpol_CtSt(t, x13Ct, St)
% interpolate x13Ct and St in 100 time points
N = 100;
dt = diff(t(:));
dtmin = min(dt(dt>0))/10000;
% add a small dt at the beginning of each cell cycle so that
% interpolation can be done in unique time-points
if size(t,2)>1
    t(1,2:size(t,2)) = t(1,2:size(t,2)) + dtmin;
end            
t100     = linspace(min(t(:)), max(t(:)), N)';
x13Ct100 = interp1(t(:),x13Ct(:), t100);
St100    = interp1(t(:),St(:),    t100);
end


function dx_val = dx(rC, x_end, C_end, t_end, avgCcell)
[~, x13Ct0, ~, ~, ~] = x13C_time0(rC, C_end, avgCcell, t_end, x_end);
dx_val = x13Ct0(end)-x_end;
end

function [Ct, x13Ct, t, tau, St] = x13C_time0(rate, C_end, avgCcell, t_incubation, xSEt)
%% model for predicting 13C atom fraction in a *linearly* growing cell
% Lubos Polerecky, 2020-06-07, Utrecht University

% note that rate corresponds to the cell-specific rate, rC

%% maximal C content of cells, fg C cell-1
% cell content varies from Cmax/2 to Cmax 
% the stable C content distribution is described by Koch 1966
% (J.gen.Microbiol.), which implies value of Cmax based on the average
Cmax = avgCcell*2*log(2);    

% end cell-cycle stage
s_end = C_end/(Cmax/2)-1;
if s_end>=1 || s_end<0
    s_end = mod(s_end,1);
    %s_end=0;
end

%% reconstruct initial cell-cycle stage (s_ini) from the one at the end (s_end)
rC = rate;
kC = rC/avgCcell;
tau = log(2)/kC;
s1 = s_end - rC*t_incubation/(Cmax/2);
% ensure that s_ini is between 0 and 1!
if s1<0
    % debug
    a=0;
end
s_ini = mod(s1,1);    

%% initial C content of the cell
Cini = Cmax/2 * (1+s_ini);

% number of time-points between division time-points
Nt = 8;

%% division time-points, linearly growing cell
% division cycles
n_div = [1:10];
% 13C fractions after division cycles, if starting from s_ini=0
x13C = 1-1./2.^(n_div-1);
% minimum number of cell cycles the cell must have gone through
n_min = find(xSEt>x13C, 1, 'last' );
n_div = n_div(1:n_min);
% division time points
% note that 0 is added artificially at the beginning to make it work
% later on
t_div = [0 (n_div*Cmax/2-(Cini-Cmax/2))/rC];

%% time interval for each (potentially partial) cell cycle
t = zeros(Nt, n_min);
for j=1:n_min
    t(:,j) = linspace(t_div(j),min([t_incubation, t_div(j+1)]),Nt)';
end
if t_incubation>t_div(n_min+1)
    t = [t linspace(t_div(n_min+1), t_incubation, Nt)'];
end

%% calculate C(t) and x(t) for each time interval
Ct=zeros(size(t));
x13Ct=zeros(size(t)); 
for j=1:size(t,2)
    if j==1
        C_ini = Cini;
        x_ini = 0;
    else
        %% this is where biomass resetting due to cell division is applied!
        C_ini = Cmax/2;
        x_ini = x13Ct(Nt,j-1);
    end
    dt = t(:,j)-t(1,j);
    %% this is where the linear growth model is applied!
    Ct(:,j) = C_ini + rC*dt;
    x13Ct(:,j) = (x_ini + rC/C_ini*dt)./(1 + rC/C_ini*dt);
end

%% determine stage in the cell-cycle (value between 0 and 1)
St = Ct/(Cmax/2)-1;
end


function [Ct, x13Ct, t, tau, St] = x13C_time1(rate, C_end, avgCcell, t_incubation)
%% model for predicting 13C atom fraction in an *exponentially* growing cell
% Lubos Polerecky, 2020-06-07, Utrecht University

% note that the input parameter *rate* corresponds to the average
% cell-specific rate, <rC> = k*avgCcell;

%% maximal C content of cells, fg C cell-1
% cell content varies from Cmax/2 to Cmax 
% the stable C content distribution is described by Koch 1966
% (J.gen.Microbiol.), which implies value of Cmax based on the average
Cmax = avgCcell/log(2);

% end cell-cycle stage
s_end = C_end/(Cmax/2)-1;
if s_end>=1 || s_end<0
    s_end = mod(s_end,1);
end

%% determine initial cell-cycle stage (s_ini) from the final one (s_end)
rC = rate;              % average cell-specific rate
kC = rC/avgCcell;       % carbon-specific rate
tau = log(2)/kC;        % doubling time
% reconstruct initial C content (first without accounting for cell
% division)
Ci = C_end*exp(-kC*t_incubation);
% determine how many times the cell divided during t_incubation
Ndiv = find(Ci<Cmax./(2.^[0:10]), 1, 'last' ) - 1;
% now s_ini can be calculated from Ci and the "diminished" critical Cmax
s_ini = Ci/(Cmax/2^(Ndiv+1)) - 1;

if s_ini<0
    % debug, but this condition should never be true
    a=0;
end

%% true initial C content of the cell determined from s_ini
Cini = Cmax/2 * (1+s_ini);

% number of time-points between division time-points
Nt = 8;

%% determine division time-points, exponentially growing cell
% division cycles
n_div = [0:Ndiv];
% division time points
t_div = t_incubation - (n_div + log(C_end/Cmax)/log(2))*tau;
% note that 0 is added artificially at the beginning to make it work
% later on
t_div = [0 sort( t_div(t_div<t_incubation) )];
n_min = length(t_div)-1;

%% time interval for each (potentially partial) cell cycle
t = zeros(Nt, n_min);
for j=1:n_min
    t(:,j) = linspace(t_div(j),min([t_incubation, t_div(j+1)]),Nt)';
end
if t_incubation>t_div(n_min+1)
    t = [t linspace(t_div(n_min+1), t_incubation, Nt)'];
end

%% calculate C(t) and x(t) for each time interval
Ct=zeros(size(t));
x13Ct=zeros(size(t)); 
for j=1:size(t,2)
    if j==1
        C_ini = Cini;
        x_ini = 0;
    else
        %% this is where biomass resetting due to cell division is applied!
        C_ini = Cmax/2;
        x_ini = x13Ct(Nt,j-1);
    end
    dt = t(:,j)-t(1,j);
    %% this is where the exponential growth model is applied!
    Ct(:,j) = C_ini*exp(kC*dt);
    x13Ct(:,j) = 1-exp(-kC*t(:,j));
end

%% determine stage in the cell-cycle (value between 0 and 1)
St = Ct/(Cmax/2)-1;
end


function [m, pdf, m_mid, m_edge] = Cdistrib_unsync(Cmax,N,order)
% Generate distribution of cell sizes (expressed as C content) in a stable
% population of cells with perfectly unsynchronized cell cycles. Cells grow
% either linearly (order=0) or exponentially (order=1). Their sizes vary in
% time between Cmax/2 and Cmax, but the probability density function (pdf)
% does not. See Koch 1966 (J.gen.Microbio.)
x = rand(N,1);
if order==1
    m = Cmax./(2-x);
    pdf = Cmax./(x.^2);
elseif order==0
    m = Cmax*( 1 - log(2-x)/log(4) );
    pdf = 8*log(2)/Cmax * exp(-2*x*log(2)/Cmax);
end
% return also mid- and edge-points used for displaying the histogram of m
N0 = 50;
m_edge = linspace(Cmax/2,Cmax,N0+1)';
dm = diff(m_edge);
m_mid = m_edge(1:N0)+dm(1)/2;
end


function [m, pdf, m_mid, m_edge] = Cdistrib_partsync(C0,dC,Cmax,N)
% Generate distribution of cell sizes (expressed as C content) in a
% population of linearly growing cells with partially synchronized cell
% cycles. Cell sizes vary in time between Cmax/2 and Cmax due to binary
% cell division, and the probability density function (pdf) is given by Eq.
% 18 (see also Fig. S1) in Polerecky et al. (submitted to FrontMic).
N0=2000;
m_edge=linspace(Cmax/2,Cmax,N0+1)';
dm=diff(m_edge);
m_mid=m_edge(1:N0)+dm(1)/2;
m_ini = m_mid;
% first term of the pdf
pdf = exp(-(m_ini-C0).^2/(2*dC^2));
for j=1:10
    % add more terms to the pdf
    tL = 1/2^j*exp(-(m_ini-j*Cmax/2-C0).^2/(2*dC^2));
    tR = 2^j*exp(-(m_ini+j*Cmax/2-C0).^2/(2*dC^2));
    pdf = pdf + tL + tR;
end
% normalize pdf
intpdf = sum( pdf.*dm );
pdf = pdf/intpdf;
%% generate random numbers that follow pdf calculated above
% first determine the cdf
pdf0 = [0; pdf];
dm0 = [dm(1); dm];
cdf = cumsum( pdf0.*dm0 );
% now generate random numbers based on the inverse cdf
x = rand(N,1);
[A, ia] = unique(cdf);
m = interp1(A,m_edge(ia), x, 'linear');
end


function m = Cdistrib_nondiv(C0,dC,N)
% Generate distribution of cell sizes (expressed as C content) in a
% population of linearly growing cells with partially synchronized cell
% cycles. No cell division is assumed here, so the distribution is simply
% based on the Gaussian distribution.
dCr = dC/C0/2;
m = exp( randn(N,1)*dCr*log(C0) + log(C0) );
end

           
function plot_results_for_current_cell(fign, Nsimul2, ...
                t100, St100_0, x13Ct100_0, St100_1, x13Ct100_1, ...
                xt_rnd, s0_ini_rnd, s1_ini_rnd, ...
                rC_A, rC_B, rC_C, rC_nondiv, ...
                avgCcell, Crand, Crand_nondiv, ...
                t_incubation, xSEt, dxSEt)
%% plot results for the current cell

% calculate 2D histograms of the predicted values for plotting
xval=t100*ones(1,Nsimul2);
xval=xval(:);
frand=figure;
h1=histogram2(xval,St100_0(:),'DisplayStyle','tile','ShowEmptyBins','on','XBinEdges',linspace(min(t100),max(t100),101),'YBinEdges',linspace(0,1,101),'Normalization','pdf');
d1=h1.Values;
h2=histogram2(xval,x13Ct100_0(:),'DisplayStyle','tile','ShowEmptyBins','on','XBinEdges',linspace(min(t100),max(t100),101),'YBinEdges',linspace(0,1,101),'Normalization','pdf');
d2=h2.Values;
% renormalize so that they can be included in one image
d1=flipud(rescale(d1'));
d2=flipud(rescale(d2'));
rgb=zeros(size(d1,1),size(d1,2),3);
thr=0.01;
r=rgb(:,:,1);
r(d2>thr)=rescale(log10(d2(d2>thr))); g=r; b=r; r(d1>thr)=rescale(log10(d1(d1>thr)));
rgb(:,:,1)=r; rgb(:,:,2)=g; rgb(:,:,3)=b;
% calculate mean values for plotting
x13Ct_mean = mean(x13Ct100_0,2,'omitnan');
St_mean = mean(St100_0,2,'omitnan');
close(frand);

% plot data
fig1=figure(fign);
set(fig1,'units','normalized','outerposition',[0 0 1 1]);        

s1=subplot(2,3,1);
hold off
% histograms of x13C(t) and s(t)
image(t100,linspace(1,0,100),1-rgb);
hold on
% average prediction, linear growth model
plot(t100, x13Ct_mean, 'r-');
% average prediction, exponential growth model
plot(t100, 1-exp(-mean(rC_A/avgCcell,'omitnan')*t100), 'b-')
% average cell cycle, linear growth model
plot(t100, St_mean, ':', 'LineWidth',2,'Color',[0 1 1]*0.8);
% experimental value
plot(t_incubation, xSEt, 'rx');
plot(t_incubation*[1;1], xSEt+dxSEt*[1; -1], 'r-','LineWidth',2);
% set axes properties
set(gca,'ydir','normal')
xlabel('t');
ylabel('x_S^E, s');
legend({'x_S^E, 0-order', 'x_S^E, 1-order', 's, 0-order', 'experimental value'},...
    'location','southeast');
ylim([0 1]);
tit=sprintf('Full prediction');
title(tit)

s2=subplot(2,3,2);
hold off;
plot(xt_rnd,Crand,'bo');
xlabel('x_j')
ylabel('C_j')
tit=sprintf('Experimental values (C-nondiv = %.2f +/- %.2f)', mean(Crand_nondiv), std(Crand_nondiv));
title(tit)

s3=subplot(2,3,3);
histogram2(xt_rnd,Crand,'DisplayStyle','tile','ShowEmptyBins','on','EdgeAlpha',0,'NumBins',[30 30],'Normalization','pdf');
xlabel('x_j')
ylabel('C_j')
title('Experimental values')

s4=subplot(2,3,4);
hold off;
histogram(s0_ini_rnd,30,'EdgeColor',[1 0 0],'FaceColor',0.9*[1 1 1]);%,'Normalization','pdf');
hold on;
Cmax = avgCcell*2*log(2);
histogram(Crand/(Cmax/2)-1,30,'EdgeColor',[0 0 1],'FaceColor',0.9*[1 1 1]);%,'Normalization','pdf');
xlabel('s');
ylabel('count');
title('Predicted values')
legend({'s_{ini}','s_{end}'})
view([90 -90])

s5=subplot(2,3,5);
histogram2(rC_C,s0_ini_rnd,'DisplayStyle','tile','ShowEmptyBins','on','EdgeAlpha',0,'NumBins',[30 30],'Normalization','pdf');
xlabel('r, avg-SIP-div');
ylabel('s_{ini}');
title('Predicted values')

s6=subplot(2,3,6);
hold off;
histogram(rC_A,'EdgeColor',[0 0 1],'FaceColor',0.9*[1 1 1],'Normalization','pdf');
hold on;
histogram(rC_B,'EdgeColor',[1 0 0],'FaceColor',0.9*[1 1 1],'Normalization','pdf');
histogram(rC_C,'EdgeColor',[1 0 1],'FaceColor',0.9*[1 1 1],'Normalization','pdf');
histogram(rC_nondiv,'EdgeColor',[0 1 0],'FaceColor',0.9*[1 1 1],'Normalization','pdf');
xlabel('r');
ylabel('count');
leg1 = sprintf('avg-cell-cycle: %.3f+/-%.3f',mean(rC_A,'omitnan'),std(rC_A,'omitnan'));        
leg2 = sprintf('end-SIP: %.3f+/-%.3f',mean(rC_B,'omitnan'),std(rC_B,'omitnan'));
leg3 = sprintf('avg-SIP-div: %.3f+/-%.3f',mean(rC_C,'omitnan'),std(rC_C,'omitnan'));
leg4 = sprintf('avg-SIP-nondiv: %.3f+/-%.3f',mean(rC_nondiv,'omitnan'),std(rC_nondiv,'omitnan'));
tit=sprintf('Predicted rate');
title(tit)
legend({leg1,leg2,leg3,leg4},'location','northwest','box','off','FontSize',10)
end

function plot_individual_result(fign, t,x13Ct,St,t_incubation,xSEtj,dxSEt,j,Nsimul)
figure(fign);
s1=subplot(2,3,1); 
hold off
plot(t, x13Ct, 'k-');
hold on
plot(t, St, 'k--');
plot(t_incubation, xSEtj, 'rx');
plot(t_incubation*[1;1], xSEtj+dxSEt*[1; -1], 'r-','LineWidth',2);
xlabel('t');
ylabel('- x_S^E(t), -- s(t)');
ylim([0 1]);
title('Model prediction')
fprintf(1,'Check figure %d to quality check the results of iteration %d/%d.\n',fign,j,Nsimul);
input('Press enter to continue (or Ctrl+c to break).');
end
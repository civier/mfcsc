function [fitresult, gof] = createFit_fcsc_LAR(sc_avg_vec_sort, fc_avg_vec_sort, is_figures)

%NOTE: Power2. Levenbrg algorithm. LAR

%CREATEFIT(SC_AVG_VEC_SORT,FC_AVG_VEC_SORT)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : sc_avg_vec_sort
%      Y Output: fc_avg_vec_sort
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 08-Dec-2022 01:36:11


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( sc_avg_vec_sort, fc_avg_vec_sort );

% Set up fittype and options.
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [0.0161523321250906 0.669803695723152 -1.05455610851554];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if is_figures
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'fc_avg_vec_sort vs. sc_avg_vec_sort', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'sc_avg_vec_sort', 'Interpreter', 'none' );
    ylabel( 'fc_avg_vec_sort', 'Interpreter', 'none' );
    grid on
end

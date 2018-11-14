classdef NOxModelsUnitTests < matlab.unittest.TestCase
    
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Use 20 NOx concentrations between 0.1 and 10 ppb (2e19 is the
        % number density of air)
        %nox_concs = logspace(-10,-8,10)*2e19;
        %nox_concs = logspace(-10,-8,5)*2e19;
        nox_concs = logspace(-9,-8,5)*2e19;
        
        % Use weekend/weekday up to a factor of 2 decrease on weekends, try
        % some cases where it increases on the weekends just to stress test
        % it.
        %nox_ratios = [0.5 0.75 1.0 1.5 2.0];
        %nox_ratios = [0.5 1.0 2.0];
        nox_ratios = 0.5;
        
        phox = 6.25e6;
        alpha = 0.04;
        %test_vocrs = 1:10;
        %test_vocrs = [1 5 10];
        test_vocrs = [1 5 10];
        
        rel_tol = 0.01;
    end
    
    methods(Test)
        function test_hox_wkday_wkend_constraint(testCase)
            wkday_nox = testCase.nox_concs;
            vocrs = testCase.test_vocrs;
            ratios = testCase.nox_ratios;
            tol = testCase.rel_tol;
            test_results = false(numel(wkday_nox), numel(vocrs), numel(ratios));
            for i_vocr = 1:numel(vocrs)
                for i_ratio = 1:numel(ratios)
                    fprintf('  Testing VOCR = %.2f, wkend/wkday NOx = %.2f\n', vocrs(i_vocr), ratios(i_ratio));
                    wkend_nox = ratios(i_ratio) .* wkday_nox;
                    [wkday.tau,~,~,~,wkday.species] = nox_lifetime(wkday_nox, 'vocr', vocrs(i_vocr), 'alpha', testCase.alpha, 'phox', testCase.phox);
                    [wkend.tau,~,~,~,wkend.species] = nox_lifetime(wkend_nox, 'vocr', vocrs(i_vocr), 'alpha', testCase.alpha, 'phox', testCase.phox);
                    for i_nox = 1:numel(wkday_nox)
                        [wkday_test, wkend_test] = hox_solve_wkday_wkend_constraint(ratios(i_ratio), wkday.tau(i_nox), wkend.tau(i_nox), testCase.phox, testCase.alpha);
                        % Verify the NOx, OH, and VOCR
                        nox_check = abs(reldiff(wkday_test.nox, wkday_nox(i_nox))) <= tol && abs(reldiff(wkend_test.nox, wkend_nox(i_nox)));
                        oh_check = abs(reldiff(wkday_test.oh, wkday.species.oh(i_nox))) <= tol && abs(reldiff(wkend_test.oh, wkend.species.oh(i_nox)));
                        vocr_check = abs(reldiff(wkday_test.vocr, vocrs(i_vocr))) <= tol && abs(reldiff(wkend_test.vocr, vocrs(i_vocr)));
                        fprintf('    [NOx] = %.2g, nox_check = %d, oh_check = %d, vocr_check = %d\n', wkday_nox(i_nox), nox_check, oh_check, vocr_check);
                        test_results(i_nox, i_vocr, i_ratio) = nox_check && oh_check && vocr_check;
                    end
                end
            end
            
            my_dir = fileparts(mfilename('fullpath'));
            save_file = fullfile(my_dir, 'hox_ratio_test.mat');
            save(save_file, 'test_results', 'wkday_nox', 'ratios', 'vocrs');
            testCase.verifyTrue(all(test_results(:)));
        end
        
        function test_hox_tau_constraint(testCase)
            nox = testCase.nox_concs;
            vocrs = testCase.test_vocrs;
            test_results = false(numel(nox), numel(vocrs));
            for i_vocr = 1:numel(vocrs)
                fprintf('  Testing VOCR = %.2f\n', vocrs(i_vocr));
                taus = nox_lifetime(nox, 'phox', testCase.phox, 'alpha', testCase.alpha, 'vocr', vocrs(i_vocr));
                for i_nox = 1:numel(testCase.nox_concs)
                    [~,~,~,calc_vocr] = hox_solve_tau_constraint(nox(i_nox), testCase.phox, taus(i_nox), testCase.alpha);
                    test_results(i_nox, i_vocr) = abs(reldiff(calc_vocr, vocrs(i_vocr))) <= 0.01;
                end
            end
            
            testCase.verifyTrue(all(test_results(:)));
        end
    end
end


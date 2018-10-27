classdef NOxModelsUnitTests < matlab.unittest.TestCase
    
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Use 20 NOx concentrations between 0.1 and 10 ppb (2e19 is the
        % number density of air)
        nox_concs = logspace(-10,-8,10)*2e19;
        phox = 6.25e6;
        alpha = 0.04;
    end
    
    methods(Test)
        function test_hox_tau_constraint(testCase)
            nox = testCase.nox_concs;
            vocrs = 1:10;
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


classdef Likelihood4Coordinates < handle
    properties
    end
    
    methods
        function obj = Likelihood4Coordinates(varargin)
        end

        function [X, Y, L] = calLikelihood4Area(obj, area_size, pos_rx, rot_abs, received_signal_cell, el_num, nPower)
            % This script computes the Maximum Likelihood function L(x_tx,y_tx) based on the received signals at two receivers.
            x_tx = linspace(0, area_size, area_size+1);
            y_tx = linspace(0, area_size, area_size+1);
            [X, Y] = meshgrid(x_tx, y_tx);
            
            L = obj.likelihoodFromCoor(X, Y, pos_rx, rot_abs, received_signal_cell, el_num, nPower);
        end

        function L = likelihoodFromCoorSet(obj, coor, pos_rx, rot_abs, received_signal_cell, el_num, nPower)
            % This script computes the Maximum Likelihood function L(x_tx,y_tx) based on the received signals at two receivers.
            x_tx = coor(1); y_tx = coor(2);
            L = obj.likelihoodFromCoor(x_tx, y_tx, pos_rx, rot_abs, received_signal_cell, el_num, nPower);
        end

        function L = likelihoodFromCoor(obj, x_tx, y_tx, pos_rx, rot_abs, received_signal_cell, el_num, nPower)
            % This script computes the Maximum Likelihood function L(x_tx,y_tx) based on the received signals at two receivers.
            % Calculate the angular components for each receiver using function f
            sin_theta1 = obj.coors2sin(x_tx, y_tx, pos_rx(1,1), pos_rx(1,2), rot_abs(1));
            sin_theta2 = obj.coors2sin(x_tx, y_tx, pos_rx(2,1), pos_rx(2,2), rot_abs(2));
            L_d0=100; d0=100; alpha=4;
            gamma1 = obj.computeGamma(x_tx, y_tx, pos_rx(1,1), pos_rx(1,2), L_d0, d0, alpha);
            gamma2 = obj.computeGamma(x_tx, y_tx, pos_rx(2,1), pos_rx(2,2), L_d0, d0, alpha);
            % Create gamma with same dimensions as sin_theta
            % gamma1 = obj.createOnes(sin_theta1)
            % gamma2 = obj.createOnes(sin_theta2)
            % Compute the likelihood function based on received signals and steering vectors
            L = obj.likelihoodFromAngles(sin_theta1, sin_theta2, gamma1, gamma2, received_signal_cell, el_num, nPower);
        end

        function cells = createOnes(~, input)
            % Creates a cell array of ones with the same size/structure as the input
            if iscell(input)
                cells = cellfun(@(x) 1, input, 'UniformOutput', false);
            else
                cells = num2cell(ones(size(input)));
            end
        end

        %% ==================================== Local Functions
        function P = likelihoodFromAngles(obj, sin_theta1, sin_theta2, gamma1, gamma2, received_signal_cell, el_num, nPower)
            % Handle multiple time instances by averaging log-likelihood
            P_total = zeros(size(sin_theta1));
            time_samples = size(received_signal_cell{1}, 2);
            
            for t = 1:time_samples
                % Get signal at current time sample
                z_t = [received_signal_cell{1}(:,t).' received_signal_cell{2}(:,t).'].';
                
                % Compute covariance matrices
                Sigma_z = obj.computeSigmaZ(sin_theta1, sin_theta2, gamma1, gamma2, nPower, el_num);
                
                % Evaluate likelihood for this time sample
                P_t = zeros(size(Sigma_z));
                for i = 1:numel(Sigma_z)
                    try
                        [L, flag] = chol(Sigma_z{i});
                        if flag == 0
                            % Cholesky decomposition succeeded
                            logdet = 2*sum(log(diag(L)));
                            % Compute quadratic form using backslash
                            quadratic_form = z_t' * (Sigma_z{i} \ z_t);
                        else
                            % Matrix not positive definite, use a fallback
                            [U, S, ~] = svd(Sigma_z{i});
                            s = diag(S);
                            % Filter out very small eigenvalues for pseudo-inverse
                            tol = max(size(Sigma_z{i})) * eps(max(s));
                            r = sum(s > tol);
                            if r < length(s)
                                % Use pseudo-inverse for likelihood calculation
                                Sigma_inv = U(:,1:r) * diag(1./s(1:r)) * U(:,1:r)';
                                logdet = sum(log(s(1:r)));
                                % Compute quadratic form using explicit multiplication
                                quadratic_form = z_t' * Sigma_inv * z_t;
                            else
                                P_t(i) = -inf; % Mark as an impossible location
                                continue;
                            end
                        end
                        
                        % Apply likelihood formula
                        P_t(i) = real(-logdet - log(pi^(2*el_num)) - quadratic_form);
                        
                    catch
                        P_t(i) = -inf; % For any calculation errors
                    end
                end
                
                % Add to total log-likelihood
                P_total = P_total + P_t;
            end
            
            % Return average log-likelihood
            P = P_total / time_samples;
        end

        function sin_theta = coors2sin(~,x_tx, y_tx, x_rx, y_rx, phi)
            % Computes sin(theta) from grid points (x_tx,y_tx) relative to Rx position (x_rx,y_rx)
            thetaphi = atan2d(y_tx - y_rx, x_tx - x_rx);
            theta = thetaphi - phi;
            sin_theta = sind(theta);
        end

        function gamma = computeGamma(~, x_tx, y_tx, x_rx, y_rx, L_d0, d0, alpha)
            % Computes the steering vector gamma at each grid point.
            gamma = arrayfun(@(x_tx, y_tx) L_d0.^(-1/2) .* d0.^(alpha/2) .* ((x_rx-x_tx)^2+(y_rx-y_tx)^2).^(-alpha/4), x_tx, y_tx, 'UniformOutput', false);
        end
        function Sigma_z = computeSigmaZ(obj, sin_theta1, sin_theta2, gamma1, gamma2, nPower, el_num)
            % Computes the covariance matrix Sigma_z at each grid point.
            a1 = obj.steerVect_ULA(sin_theta1, el_num);
            a2 = obj.steerVect_ULA(sin_theta2, el_num);
            
            block1 = cellfun(@(a1, gamma1) abs(gamma1)^2 * (a1 * a1') + nPower * eye(el_num), a1, gamma1, 'UniformOutput', false);
            block2 = cellfun(@(a1, a2, gamma1, gamma2) gamma1 * conj(gamma2) * (a1 * a2'), a1, a2, gamma1, gamma2, 'UniformOutput', false);
            block3 = cellfun(@(a1, a2, gamma1, gamma2) gamma2 * conj(gamma1) * (a1 * a2'), a2, a1, gamma1, gamma2, 'UniformOutput', false);
            block4 = cellfun(@(a2, gamma2) abs(gamma2)^2 * (a2 * a2') + nPower * eye(el_num), a2, gamma2, 'UniformOutput', false);
            Sigma_z = cellfun(@(b1, b2, b3, b4) [b1, b2; b3, b4], ...
                block1, block2, block3, block4, 'UniformOutput', false);
            % Add a small identity term to each for numerical stability.
            % Sigma_z = cellfun(@(S) S + 1e-9 * eye(size(S)), Sigma_z, 'UniformOutput', false);
        end

        function a = steerVect_ULA(~, sin_theta, el_num)
            % Generates the ULA steering vector for each grid point.
            % sin_theta is assumed to be a matrix.
            a = arrayfun(@(sin_theta) exp(-1i * pi * (0:(el_num-1))' * sin_theta), sin_theta, 'UniformOutput', false);
        end
    
    end

end

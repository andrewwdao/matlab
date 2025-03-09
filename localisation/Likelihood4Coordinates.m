classdef Likelihood4Coordinates < handle
    properties
    end
    
    methods
        function obj = Likelihood4Coordinates(varargin)
        end

        function [X, Y, L] = calLikelihood4Area(obj, area_size, pos_rx, rot_abs, received_signal_cell, el_num, nPower)
            % This script computes the Maximum Likelihood function L(x_tx,y_tx) based on the received signals at K receivers.
            x_tx = linspace(0, area_size, area_size+1);
            y_tx = linspace(0, area_size, area_size+1);
            [X, Y] = meshgrid(x_tx, y_tx);
            
            L = obj.likelihoodFromCoor(X, Y, pos_rx, rot_abs, received_signal_cell, el_num, nPower);
        end

        function L = likelihoodFromCoorSet(obj, coor, pos_rx, rot_abs, received_signal_cell, el_num, nPower)
            % This script computes the Maximum Likelihood function L(x_tx,y_tx) based on the received signals at K receivers.
            x_tx = coor(1); y_tx = coor(2);
            L = obj.likelihoodFromCoor(x_tx, y_tx, pos_rx, rot_abs, received_signal_cell, el_num, nPower);
        end

        % filepath: /d:/workspaces/matlab/localisation/Likelihood4Coordinates.m
        function L = likelihoodFromCoor(obj, x_tx, y_tx, pos_rx, rot_abs, received_signal_cell, el_num, nPower)
            % This script computes the Maximum Likelihood function L(x_tx,y_tx) based on the received signals at K receivers.
            num_rx = size(pos_rx, 1);
            L_d0=100; d0=100; alpha=4;  % Parameters for the gamma function
            % Dynamically calculate for each receiver
            sin_theta_cell = cell(1, num_rx);
            gamma_cell = cell(1, num_rx);
            for rx_idx = 1:num_rx
                % Calculate the angular components
                sin_theta_cell{rx_idx} = obj.coors2sin(x_tx, y_tx, pos_rx(rx_idx,1), pos_rx(rx_idx,2), rot_abs(rx_idx));
                % Calculate gamma values
                gamma_cell{rx_idx} = obj.computeGamma(x_tx, y_tx, pos_rx(rx_idx,1), pos_rx(rx_idx,2), L_d0, d0, alpha);
                % gamma_cell{rx_idx} = obj.createOnes(sin_theta_cell{rx_idx});  % Create gamma with same dimensions as sin_theta - for testing gamma=1
            end
            
            % Compute the likelihood function based on received signals and steering vectors
            L = obj.likelihoodFromAngles(sin_theta_cell, gamma_cell, received_signal_cell, el_num, nPower);
        end

        function cells = createOnes(~, input)
            % Creates a cell array of ones with the same size/structure as the input
            if iscell(input)
                cells = cellfun(@(x) 1, input, 'UniformOutput', false);
            else
                cells = num2cell(ones(size(input)));
            end
        end

        %% ==================================== Local Functions ====================================
        % function safetyDistanceCheck(~, x_tx, y_tx, pos_rx)
        %     TX_SAFETY_DISTANCE = 2; % Define safety distance if not already defined
            
        %     % Only perform check for scalar inputs (single point evaluation)
        %     if isscalar(x_tx) && isscalar(y_tx)
        %         % Check distance to each receiver
        %         for i = 1:size(pos_rx, 1)
        %             distance = sqrt((x_tx - pos_rx(i,1))^2 + (y_tx - pos_rx(i,2))^2);
                    
        %             % Issue warning if too close
        %             if distance < TX_SAFETY_DISTANCE
        %                 warning('RX too close to TX %d (distance: %.2f m) - might cause numerical instability', i, distance);
        %             end
        %         end
        %     end
        %     % Skip the check for meshgrid inputs (likelihood map generation)
        % end
        
        function P = likelihoodFromAngles(obj, sin_theta_cell, gamma_cell, received_signal_cell, el_num, nPower)
            % Handle multiple time instances by averaging log-likelihood
            P_total = zeros(size(sin_theta_cell{1}));
            time_samples = size(received_signal_cell{1}, 2);
            
            for t = 1:time_samples
                % Get signal at current time sample
                z_t = cell2mat(cellfun(@(x) x(:,t), received_signal_cell, 'UniformOutput', false)); % z_t = [w_1t' ... wKt']' is a column vector
                
                % Compute covariance matrices using cell arrays
                Sigma_z = obj.computeSigmaZ(sin_theta_cell, gamma_cell, nPower, el_num);
                
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
                        P_t(i) = real(-logdet - quadratic_form);
                        
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

        function Sigma_z = computeSigmaZ(obj, sin_theta_cell, gamma_cell, nPower, el_num)
            % Computes the covariance matrix Sigma_z at each grid point using cell arrays
            num_rx = length(sin_theta_cell);
            % Get steering vectors for all receivers
            a_cell = cellfun(@(sin_theta) obj.steerVect_ULA(sin_theta, el_num), sin_theta_cell, 'UniformOutput', false);
            % Check if we're dealing with a single point or a grid
            isSinglePoint = ~iscell(a_cell{1});
            
            % create a single covariance matrix for a single point
            if isSinglePoint
                Sigma_z = {obj.buildCovarianceMatrix(a_cell, gamma_cell, nPower, el_num, num_rx)};
            % Grid case - create a cell array of covariance matrices for map visualisation
            else 
                % Determine the number of grid points
                num_el = numel(a_cell{1});
                % Initialize the cell array for full covariance matrices
                Sigma_z = cell(size(a_cell{1}));
                
                % Process each grid point
                for p = 1:num_el
                    % Extract steering vectors and gammas for this point
                    a_point = cell(1, num_rx);
                    gamma_point = cell(1, num_rx);
                    for i = 1:num_rx
                        a_point{i} = a_cell{i}{p};
                        gamma_point{i} = gamma_cell{i}{p};
                    end
                    % Build covariance matrix for this point
                    Sigma_z{p} = obj.buildCovarianceMatrix(a_point, gamma_point, nPower, el_num, num_rx);
                end
            end
        end

        function cov_matrix = buildCovarianceMatrix(~, a_vectors, gamma_values, nPower, el_num, num_rx)
            % Helper function to build a covariance matrix for a single point
            
            % Initialize the block matrix
            cov_matrix = zeros(num_rx * el_num, num_rx * el_num); % size: (num_rx * el_num) x (num_rx * el_num)
            for i = 1:num_rx  % Fill in each block of the covariance matrix
                for j = 1:num_rx
                    % Extract steering vectors and gamma values
                    a_i = a_vectors{i};
                    a_j = a_vectors{j};
                    gamma_i = gamma_values{i};
                    gamma_j = gamma_values{j};
                    
                    % Calculate block (i,j)
                    if i == j
                        % Diagonal block: signal covariance + noise
                        block_ij = abs(gamma_i)^2 * (a_i * a_i') + nPower * eye(el_num);
                    else
                        % Off-diagonal block: cross-correlation between receivers
                        block_ij = gamma_i * conj(gamma_j) * (a_i * a_j');
                    end
                    
                    % Insert this block into the full covariance matrix
                    row_indices = ((i-1)*el_num + 1):(i*el_num);
                    col_indices = ((j-1)*el_num + 1):(j*el_num);
                    cov_matrix(row_indices, col_indices) = block_ij;
                end
            end
            
            % Add small identity term for numerical stability
            cov_matrix = cov_matrix + 1e-6 * eye(size(cov_matrix));
        end

        function a = steerVect_ULA(~, sin_theta, el_num)
            % Generates the ULA steering vector for each grid point.
            % sin_theta is assumed to be a matrix.
            a = arrayfun(@(sin_theta) exp(-1i * pi * (0:(el_num-1))' * sin_theta), sin_theta, 'UniformOutput', false);
        end
    
    end

end

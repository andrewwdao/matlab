classdef Likelihood4Coordinates < handle
    properties
    end
    
    methods
        function obj = Likelihood4Coordinates(varargin)
        end

        function [X, Y, L] = calLikelihood4Area(obj, area_size, pos_rx, rot_abs, received_signal_cell, el_num, nPower)
            % This script computes the Maximum Likelihood function L(x,y) based on the received signals at two receivers.
            %% Compute ML Function L(x,y)
            x = linspace(0, area_size, area_size+1);
            y = linspace(0, area_size, area_size+1);
            [X, Y] = meshgrid(x, y);
            
            L = obj.calLikelihoodFromCoors(X, Y, pos_rx, rot_abs, received_signal_cell, el_num, nPower);
        end

        function L = calLikelihood4fmincon(obj, coor, pos_rx, rot_abs, received_signal_cell, el_num, nPower)
            % This script computes the Maximum Likelihood function L(x,y) based on the received signals at two receivers.
            %% Compute ML Function L(x,y)
            x = coor(1);
            y = coor(2);
            L = obj.calLikelihoodFromCoors(x, y, pos_rx, rot_abs, received_signal_cell, el_num, nPower);
        end

        function L = calLikelihoodFromCoors(obj, x, y, pos_rx, rot_abs, received_signal_cell, el_num, nPower)
            % This script computes the Maximum Likelihood function L(x,y) based on the received signals at two receivers.
            %% Compute ML Function L(x,y)
            % Calculate the angular components for each receiver using function f
            sin_theta1 = obj.coors2sin(x, y, pos_rx(1,1), pos_rx(1,2), rot_abs(1));
            sin_theta2 = obj.coors2sin(x, y, pos_rx(2,1), pos_rx(2,2), rot_abs(2));

            % Compute the likelihood function based on received signals and steering vectors
            L = obj.calLikelihoodFromAngles(sin_theta1, sin_theta2, received_signal_cell, el_num, nPower);
        end

        %% ==================================== Local Functions
        function sin_theta = coors2sin(~,x, y, x_rx, y_rx, phi)
            % Computes sin(theta) from grid points (x,y) relative to Rx position (x_rx,y_rx)
            thetaphi = atan2d(y - y_rx, x - x_rx);
            theta = thetaphi - phi;
            sin_theta = sind(theta);
        end

        function P = calLikelihoodFromAngles(obj, sin_theta1, sin_theta2, received_signal_cell, el_num, nPower)
            % Computes the likelihood function over the grid.
            % Step 1: Combine received signals into a vector.
            z = [received_signal_cell{1}.' received_signal_cell{2}.'].';
            % Step 2: Compute covariance matrices for each grid point.
            Sigma_z = obj.computeSigmaZ(sin_theta1, sin_theta2, nPower, el_num);
            % Step 3: Evaluate the log Likelihood function at each point.
            % P = cell2mat(cellfun(@(Sig) 1/(pi^(2*el_num)*det(Sig)) * exp(-z.' / Sig * z), ...
            %     Sigma_z, 'UniformOutput', false));
            P = cell2mat(cellfun(@(Sig) real(-log(pi^(2*el_num)*det(Sig)) - (z.' / Sig * z)), ...
                Sigma_z, 'UniformOutput', false));
        end

        function Sigma_z = computeSigmaZ(obj, sin_theta1, sin_theta2, nPower, el_num)
            % Computes the covariance matrix Sigma_z at each grid point.
            gamma1 = 1; gamma2 = 1;
            a1 = obj.steerVect_ULA(sin_theta1, el_num);
            a2 = obj.steerVect_ULA(sin_theta2, el_num);
            
            block1 = cellfun(@(vec) abs(gamma1)^2 * (vec * vec') + nPower * eye(el_num), a1, 'UniformOutput', false);
            block2 = cellfun(@(vec1, vec2) gamma1 * conj(gamma2) * (vec1 * vec2'), a1, a2, 'UniformOutput', false);
            block3 = cellfun(@(vec1, vec2) gamma2 * conj(gamma1) * (vec1 * vec2'), a2, a1, 'UniformOutput', false);
            block4 = cellfun(@(vec) abs(gamma2)^2 * (vec * vec') + nPower * eye(el_num), a2, 'UniformOutput', false);
            Sigma_z = cellfun(@(b1, b2, b3, b4) [b1, b2; b3, b4], ...
                block1, block2, block3, block4, 'UniformOutput', false);
            % Add a small identity term to each for numerical stability.
            Sigma_z = cellfun(@(S) S + 1e-9 * eye(size(S)), Sigma_z, 'UniformOutput', false);
        end

        function a = steerVect_ULA(~, sin_theta, el_num)
            % Generates the ULA steering vector for each grid point.
            % sin_theta is assumed to be a matrix.
            a = arrayfun(@(x) exp(-1i * pi * (0:(el_num-1))' * x), sin_theta, 'UniformOutput', false);
        end
    
    end

end

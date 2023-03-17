function [R, h] = intcor(u,y)

% Computes the periodic cross-correlation Ruy(h) between two signals u and
%
% Inputs:
%   u & y-  periodic signals
% Outputs:
%   R - the periodic cross-correlation between u and y
%   h - the lag values at which R is computed

    % Determine the length of the input signals
    M = length(u);
    N= length(y);

    % Compute the cross-correlation between u and y for each possible
    % value of the lag h
    R = zeros(1, M);
    for h = 1:M
        % Shift y by h samples to the right
        y_shifted = circshift(y, h-1);

        % Compute the dot product between u and the shifted y
        R(h) = sum(u .* y_shifted);
    end

    % Compute the range of lag values
    h = 1:M;

    % Normalize R by the length of the input signals
    R = R / M;
end

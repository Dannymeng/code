function [P,Q,W,B] = nipals(X,Y,a)
%% NIPALS implements an algorithm for PLSR
% Applies the NIPALS algorithm for PLSR (partial least-squares regression),
% (mostly) as described in Geladi (1986).
% Note that the routine does not do any centering or scaling.
%
% |[P,Q,W,B] = nipals(X,Y,a)|
%
% Inputs:
%
% * _X_ is a matrix of inputs
% * _Y_ is a matrix of outputs
% * _a_ is the number of principal components we wish to find
%
% Outputs:
% * _P_ is a matrix of 

%% Make a copy of the original X & Y
X0 = X;
Y0 = Y;

%% Determine the size of the matrices
n = size(X,1); % The number of observations
m = size(X,2); % The number of inputs
p = size(Y,2); % The number of outputs
% assert(a < n, ['Not enough observations to obtain the '...
%     'requested number of principal components'] );

tol = 1e-6; % The tolerance for the convergence test

%% Initialize the outputs
P = zeros(m,a); % The principal components of _X_
Q = zeros(p,a); % The principal components of _Y_
W = zeros(m,a); % The _X_ weights
B = zeros(a,a); % A vector of regression coefficients
T = zeros(n,a); % The loadings of _X_

Yres = Y;

%% Determine each of the principal components
for h=1:a
    tOld = zeros(n,1);
    u = Yres(:,1);
    while 1
        wOld = X'*u;
        w = wOld/norm(wOld);
        t = X*w;%/(w'*w);
        qOld = Yres'*t;
        q = qOld/norm(qOld);
        u = Yres*q;
        % Test for convergence
        err = norm(tOld - t)/norm(t);
        if err < tol
            % Calculate the final values for this component
            pOld = X'*t;
            p = pOld/norm(pOld);
            % Store component vectors into the arrays
            P(:,h) = p;
            Q(:,h) = q;
            W(:,h) = w;
            T(:,h) = t;
            % One difference from Geladi's algorithm is that
            % we never calculate the X-residual, so we need
            % to calculate the regression coefficients a
            % bit differently. (See Eq. 14 in Geladi.)
            B(1:h,h) = (T(:,1:h)'*T(:,1:h)) \ T(:,1:h)'*u;
            Ypred = T*B*Q';
            Yres = Y - Ypred;
            break
        else
            tOld = t;
        end
    end
end
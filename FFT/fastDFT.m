function X = fastDFT(inVec)

% fastDFT
% Implement the radix-2, decimation-in-time FFT algorithm
%
% input:  inVec - discrete time vector
% output: X     - fourier coefficient vector for the input
%
% author: Zichao Wang
% July 31, 2015


% get length of input vector
N = length( inVec );

% initialize DFT vector
X = zeros( 1, N );

% twiddle factor
W = exp( -1j * 2 * pi / N);


if N == 2
    % base case
    X(1) = inVec(1) + W^0 * inVec(2);
    X(2) = inVec(1) + W^1 * inVec(2);

else
    % recursive procedure
    
    % even samples (in matlab is odd because matlab start index is 1)
    E = fastDFT( inVec( 1:2:end ) );
    
    % odd samples (in matlab is even because matlab start index is 1)
    O = fastDFT( inVec( 2:2:end ) );

    % calculate fourier coefficients
    for k = 1:N/2
        X(k)       = E(k) + W^( k - 1 ) * O(k);
        X(k + N/2) = E(k) + W^( k - 1 + N/2 ) * O(k);
    end

end

end
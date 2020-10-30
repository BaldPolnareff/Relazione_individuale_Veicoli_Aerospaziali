% Implementation of the Jacobi iteration method

function [x_appr] = Jacobi (A, b, tol, max_iter)

	% Check wether A is a square matrix
	[N M] = size (A);
	if N ~= M
		error ('A should be a square Matrix!')
	end

	% Check whether A is diagonally dominant: this should guarantee convergence
	% For each row, the sum of the absolute value of all elements (minus the diagonal element) should be smaller than the diagonal element itself

	for i = 1 : N
		row = abs (A (i, :));
		comp = sum (row) - row (i);
		if row (i) <= comp 
			error ('A should be a diagonally dominant matrix')
		end
	end

	% Check whether b is a column vector of a compatible size
	[n m] = size (b);
	if ~(n == N && m == 1)
		error ('b should be a column vector compatible with A')
	
	else
	D = diag (diag (A));	% Matrix containing A's main diagonal and zeros
	E = tril (A, -1); 	% Matrix containing A's lower triangle (second argument exludes main diagonal)
	F = triu (A, +1); 	% Same as before, but with A's upper triangle
	LU = E + F;
		
	x_0 = zeros (N, 1); % First step of iteration, arbitrarily initialized with zeros
			
	% Begin iteration
	k = 0; x_appr = x_0; difference = ones (N, 1); x_prev = x_0;
	err = norm (difference);
	while abs (err) >= tol
		x_prev = x_appr;
		x_appr = D \ (b - LU * x_prev);
		difference = x_appr - x_prev;
		err = norm (difference);
		if k > max_iter
			error ('Iterations limit exceeded, you can change it manually')
		end
		k = k + 1;
	end
		disp ('The number of iterations is'), k
	end
end


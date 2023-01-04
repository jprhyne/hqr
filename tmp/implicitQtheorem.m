
   clear

   n = 10;

   A = randn(n);

   u1 = randn(n,1); u1 = u1/norm(u1,2);
   [U,~]=qr(u1);
   if ( norm( u1 - U(:,1), 'fro' ) > 1.0 ), U(:,1) = - U(:,1); end
   [Q,G]=hess(U'*A*U);
   U = U * Q;
   fprintf('1. repres = %6.2e\n', norm( A - U*G*U', 'fro' ) / norm( A, 'fro' ));
   fprintf('2. repres = %6.2e\n', norm( U'*A*U - G, 'fro' ) / norm( A, 'fro' ));
   fprintf('3. orthog = %6.2e\n', norm( eye(n) - U'*U, 'fro' ));
   fprintf('4. orthog = %6.2e\n', norm( eye(n) - U*U', 'fro' ));
   fprintf('5. u1?    = %6.2e\n',  norm( u1 - U(:,1), 'fro' ));


   [V,H]=hess(A);
   fprintf('1. repres = %6.2e\n', norm( A - V*H*V', 'fro' ) / norm( A, 'fro' ));
   fprintf('2. repres = %6.2e\n', norm( V'*A*V - H, 'fro' ) / norm( A, 'fro' ));
   fprintf('3. orthog = %6.2e\n', norm( eye(n) - V'*V, 'fro' ));
   fprintf('4. orthog = %6.2e\n', norm( eye(n) - V*V', 'fro' ));



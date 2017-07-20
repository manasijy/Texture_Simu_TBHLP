function [T_Solution] = calcLPSolution(b,c0,A)


%%
c = c0;
k = 0;
tol=1e-6; 
maxit=50; 

%%
[m,n]=size(A); b=b(:); c=c(:); it=0;  %b and c are made column matrix
if (length(c)~=n || length(b)~=m),error('wrong dimensions'); end
D=sign(sign(b)+.5); % 0.5 could be reduced if addition causes change in the sign of elements
D = diag(D);                    % initial (inverse) basis matrix
A = [A D];                      % incorporate slack/artificial variables.Utilizing 2 phase method

%%
B = n+1:n+m;    % initial basis                
N = 1:n;        % intitial nonbasis                
phase=1; xb=abs(b); s=[zeros(n,1);ones(m,1)];   % supercost 

%%
while phase<3
   df=-1; t=inf;
   yb= D'*s(B);  
   while (it < maxit)
      if isempty(N), break, end    
      r = s(N) - [A(:,N)]'*yb;      
      [rmin,q] = min(r);            
                                   
      if rmin>=-tol*(norm(s(N),inf)+1), break, end % optimal!            
      it=it+1;
      if df>=0                      % apply Bland's rule to avoid cycling
         if maxit==inf
            disp(['LINPROG(',int2str(it),'): warning! degenerate vertex']);
         end
         J=find(r<0); Nq=min(N(J)); q=find(N==Nq); % This line gets the firstmost column in NB which gives a -ve reduced cost. 
                                    % q is the location of that column in N                                                                                                     
      end 
      d = D*A(:,N(q)); % Descent Direction: d = -inv(B)*Nq.If dq >= 0, then the linear program is unbounded. Here d<=0 =>unbounded
      I=find(d>tol); %find out the entry no of those d values which are -ve
      if isempty(I), disp('Solution is unbounded'); it=-it; break; end
      xbd=xb(I)./d(I); [r,p]=min(xbd); p=I(p); % r is the step length alpha, p is position of that alpha, p will be going out
      if df>=0                     % apply Bland's rule to avoid cycling
         J=find(xbd==r); Bp=min(B(I(J))); p=find(B==Bp); 
      end 
      xb= xb - r*d; xb(p)=r;        % update x, bsic var xb(p) is replaced with alpha 
      df=r*rmin;                    % change in f 
      v = D(p,:)/d(p);              % row vector
      yb= yb + v'*( s(N(q)) - d'*s(B) );
      d(p)=d(p)-1;
      D = D - d*v;                  % update inverse basis matrix
      t=B(p); B(p)=N(q); % pth column of B is replaced with qth column of N. q is the one with minimum reduced cost. it has come in. 
      if t>n+k, N(q)=[]; else N(q)=t; end %The pth column of B is brought into N.
                                          %If t was artificial var then one col of N is reduced:single phase
   end                             % end of phase
   xb=xb+D*(b-A(:,B)*xb);          % iterative refinement
   I=find(xb<0);                   % must be due to rounding error
   if I, xb(I)=xb(I)-xb(I); end    % so correct
   if phase==2 || it<0, break; end % B, xb,n,m,res=A(:,B)*xb-b %| changed to ||
   if xb'*s(B)>tol,it=-it; disp('no feasible solution'); break;  end % if 
   phase=phase+1;                  % Go to Phase 2
   s=1e6*norm(c,'inf')*s; s(1:n)=c;% tol=tol*norm(s,inf);
end
x=sparse(n,1); x(B)=xb; x=x(1:n); if n<21, x=full(x); end
fval=c'*x;
if it>=maxit, disp('too many iterations'); it=-it; end

OP.B = B; OP.N = N; OP.D = D; OP.r = r; OP.xb = xb; OP.yb = yb; OP.var = var(xb); 
        
%% Getting multiple solutions and passing on the plastic spin as output

T_Solution = treeSol_function(OP,b,c0);


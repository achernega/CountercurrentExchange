% Multistage Chemical Extraction | A.Chernega J.Brown J.Anderson
clc;clear all;close all;format longG;
fprintf("-----------------------------------------------------\n");
fprintf("----------- Multistage Chemical Extraction ----------\n");
fprintf("-----------------------------------------------------\n\n");

% Input from user; prompts following values
%W = input("Enter water flow rate: ");
%S = input("Enter solvent stream flow rate: ");
%xin = input("Enter mass fraction of chemical in input water stream: ");
%yin = input("Enter mass fraction of chemical in input solvent stream: ");
%m = input("Enter ratio of mass fracation in solvent stream to water stream: ");
%n = input("Enter number of reactor stages: ");
%performw = input("For SOR: would you like to try many w values? 1 = YES; 2 = NO");
% Building of matrix A and vector b from input values

% <Hardcoded Values> Use for testing! For SPD, have W = S*m
W=200;S=50;xin=0.075;yin=0;m=4;n=15;performw = 1; % Set perform w to 1 if you want to run many values of w in SOR
%W=300;S=40;xin=0.075;yin=0;m=15;n=6;performw = 2;
% </Hardcoded Values>
              
A = zeros(n,n);
b = zeros(n,1);
diagDiag = -(W + S*m);
for i = 1 : n
  A(i,i) = diagDiag;
end
topDiag = S*m;
for i = 1 : n-1
  A(i,i+1) = topDiag;
end
bottomDiag = W;
for i = 1 : n-1
  A(i+1,i) = bottomDiag;
end
b(1) = (W*xin); % Note that the negatives of A and b
b(n) = (S*yin); % are are required because otherwise
A=-A           % A will never be SPD!
% Following code uses MatLab to find a solution
% which we will use as our 'true' solution
x_true = A\b;

% Tolerance declaration!
tol = 5*(10^(-6));
maxim = 1000;

% Symmmetric Positive Definite test - SPD --------------------------------------
symCheck = false;
eigCheck = false;
SPD = false;
minEig = min(eig(A));

if A == A'
    symCheck = true;
end
if minEig >= 0
    eigCheck = true;
end
if symCheck == true && eigCheck == true
    SPD = true;
end
fprintf("\nSPD Test: done\n");

% Ill-conditioning test - kA ---------------------------------------------------

[n_k, toss_k] = size(A);
invA_k = inv(A);
x1_k = zeros(n,1);
x2_k = zeros(n,1);
sum1_k = 0;
sum2_k = 0;
fprintf("\nIll-conditioning Test: done\n");

% Find sum of absolute value of the rows of A and invA
for row_k=1:n_k   
    sum1_k = 0;
    sum2_k = 0;
    for col_k=1:n_k
        sum1_k = sum1_k + abs(A(row_k, col_k));
        sum2_k = sum2_k + abs(invA_k(row_k, col_k));
    end
    x1(row_k) = sum1_k;
    x2(row_k) = sum2_k;
end

% Determine norms of A and invA, and calc condition number, k.
norm1_k = max(x1);
norm2_k = max(x2);
k_k = norm1_k * norm2_k;
fprintf("Condition number k(A) = %.2f\n", k_k)

% Jacobi Method ----------------------------------------------------------------

[n,toss_jac] = size(A);
xOld_jac = zeros(n,1);
xNew_jac = zeros(n,1);
errs_jac = [];
convs_jac = [];
conv_jac = 0;
absErrs_jac = [];
absErr_jac = 0;
iter_jac = 0;
sum1_jac = 0;
sum2_jac = 0;
err_jac = 1;
xArray_jac = [norm(xOld_jac , 2)];
xAct_jac = A\b;

%%  Implement Method

%Split A to build/test T
D = zeros(n,n);
L = zeros(n,n);
U = zeros(n,n);

%   Create diagonal matrix (D)
for Pass=1:n
    D(Pass, Pass) = A(Pass, Pass);
end

%   Create lower triangulare matrix (L)
for Pass=1:n-1
    for row=Pass+1:n
        L(row, Pass) = -A(row,Pass);
    end
end

%   Create upper traingular matrix (U)
for Pass=2:n
    for row=1:Pass-1
        U(row, Pass) = -A(row,Pass);
    end
end   

%Build T and test Spec Radius < 1
T_jac = inv(D) * (L+U);
    
%   Pull spec radius
e_jac = eig(T_jac);
specRad_jac = max(abs(e_jac));
if specRad_jac >= 1
    fprintf("The spectral radius for this matrix is: %.2f.  This is greater than or equal to 1 and convergence is not guaranteed for iterative methods using this matrix T.", specRad);
end

%Run Jacobi Method
while (err_jac > tol) && (iter_jac < maxim)
    %Previous xNew becomes xOld on each iteration.
    temp_jac = xNew_jac; 
    xOld_jac = temp_jac;    
    for Pass = 1:n %1 pass for each value in solution vector x
    sum1_jac = 0; %Reset sums for each pass
    sum2_jac = 0;
        for col = 1:n %Calculate first sum needed for method.
            %Build first sum iterating from 1 to Pass - 1
            if col < Pass
                sum1_jac = sum1_jac + A(Pass, col)*xOld_jac(col);
            end
            %Build second sum iterating freom Pass + 1 to n
            if col > Pass
                sum2_jac = sum2_jac + A(Pass, col)*xOld_jac(col);
            end
        end
        %Calc/Assign the pass'th value of solution vector x
        xNew_jac(Pass) = (1/A(Pass,Pass)) * ( b(Pass) - sum1_jac - sum2_jac);
    end
    %Update errors, error/convergence vectors for plots, etc., on each
    %iteration of while loop.
    err_jac = max(abs(xNew_jac - xOld_jac));
    errs_jac = [errs_jac ; err_jac];
    iter_jac = iter_jac + 1;
    conv_jac = max(abs(xNew_jac - xAct_jac)) / max(abs(xOld_jac - xAct_jac));
    convs_jac = [convs_jac; conv_jac];
    absErr_jac = max(abs(xNew_jac - xAct_jac));
    absErrs_jac = [absErrs_jac; absErr_jac];
    xArray_jac = [xArray_jac ; norm(xNew_jac,2)];
end
% Output to use from this code:
x_jac_final = xNew_jac;
count_jac = iter_jac;
errs_jac;
absErrs_jac;
convs_jac;
fprintf("\nJacobi Method: done\n");
fprintf("Spectral Radius = %.2f\n", specRad_jac);

% SOR Method -------------------------------------------------------------------

[m,n]=size(A);
xOld_sor = zeros(n,1);
xNew_sor = zeros(n,1);
errs_sor = [];
convs_sor = [];
conv_sor = 0;
absErrs_sor = [];
absErr_sor = 0;
iter_sor = 0;
iterArr_sor = [];
wArr_sor = [];
err_sor = 1;
xArray_sor = [norm(xOld_sor , 2)];
xAct_sor = A\b;

w_sor = 0.05;
%%checking for positive definite, if so uses optimal relaxation param
while w_sor <= 1.95
if performw == 2
    if ((specRad_jac^2 < 1) && all(eig(A) > eps))
        w_sor = 2 / (1+ sqrt(1- specRad_jac ^2))
    else
        w_sor = 1
    end
end
err_sor=1;
xOld_sor=zeros(n,1);
xNew_sor=zeros(n,1);
iter_sor=0;


%Build T and test Spec Radius < 1
T_SOR = inv(((1/w_sor)*D) - L) * ((((1/w_sor) - 1)*D) + U);   

%   Pull spec radius
e_SOR = eig(T_SOR);
specRad_SOR = max(abs(e_SOR));
if specRad_SOR >= 1
    fprintf("The spectral radius for this matrix is: %.2f.  This is greater than or equal to 1 and convergence is not guaranteed for iterative methods using this matrix T.", specRad);
end %end of specradius build/check

while err_sor > tol && iter_sor < maxim
    for i=1:n
        SL_sor=0;
        SU_sor=0;
        for j=1:i-1
            SL_sor=SL_sor+A(i,j)*xNew_sor(j);
        end
        for j=i+1:n
            SU_sor=SU_sor+A(i,j)*xOld_sor(j);
        end
        xNew_sor(i)=(1-w_sor)*xOld_sor(i)+w_sor*(b(i)-SL_sor-SU_sor)/A(i,i);    
    end
    err_sor=norm(xNew_sor-xOld_sor, inf);
    
    errs_sor = [errs_sor ; err_sor];
    conv_sor = max(abs(xNew_sor - xAct_sor)) / max(abs(xOld_sor - xAct_sor));
    convs_sor = [convs_sor; conv_sor];
    absErr_sor = max(abs(xNew_sor - xAct_sor));
    absErrs_sor = [absErrs_sor; absErr_sor];
    xArray_sor = [xArray_sor ; norm(xNew_sor,2)];

     
    xOld_sor=xNew_sor;
    iter_sor=iter_sor+1;
end   
    iterArr_sor = [iterArr_sor ; iter_sor];
    wArr_sor = [wArr_sor ; w_sor];
w_sor = w_sor + 0.05;
% Terminates if user does not want to try a bunch of values of w
if performw == 2
  w_sor = 2;
end
end
% Output to use from this code:
x_sor_final = xNew_sor;
count_sor = iter_sor;
errs_sor;
absErrs_sor;
convs_sor;
wArr_sor;
iterArr_sor;
xArray_sor;

fprintf("\nSOR Method: done\n");
fprintf("Spectral Radius = %.2f\n", specRad_SOR);
fprintf("w = %.2f\n", (2 / (1+ sqrt(1- specRad_jac ^2))));

% Gaussian w/ Scaled Partial Pivoting Method -----------------------------------

E_gau=[A b];
[n_gau,k_gau]=size(E_gau);
E1_gau=E_gau;
E_gau=E1_gau;
s_gau=zeros(n,1);
for i=1:n
    s_gau(i)=max(abs(E_gau(i,1:n)));
end
r_gau=[1:n]';
for pass=1:n-1
    [M_gau,I_gau]=max(abs(E_gau(r_gau(pass:n), pass))./s_gau(r_gau(pass:n)));
    I_gau=I_gau+pass-1;
    if I_gau>pass
        rr_gau=r_gau(pass);
        r_gau(pass)=r_gau(I);
        r_gau(I)=rr_gau;
    end
    r_gau = transpose(r_gau);

    for row=pass+1:n
        m_gau=-E_gau(r_gau(row), pass)/E_gau(r_gau(pass),pass);
        E_gau(r_gau(row),pass)=0;
        for col=pass+1:n+1
            E_gau(r_gau(row), col)=E_gau(r_gau(row),col)+m_gau*E_gau(r_gau(pass),col);
        end
    end
end
% Back substitution
x_gau=zeros(n,1);
x_gau(n)=E_gau(r_gau(n),n+1)/E_gau(r_gau(n),n);
for row=n-1:-1:1
    sum_gau=E_gau(r_gau(row), n+1);
    for col=row+1:n
        sum_gau=sum_gau-E_gau(r_gau(row),col)*x_gau(col);
    end
    x_gau(row)=sum_gau/E_gau(r_gau(row),row);
end
x_gau = x_gau;
fprintf("\nGaussian Elimination: done\n");

% Conjugate Gradient Method ----------------------------------------------------
% if SPD -> Perform ConjGrad

if SPD
  fprintf("\nMatrix is SPD; performing conjugate gradient method...\n");
  count_cg = 0;
  x_cg = zeros(n,1);
  r_cg = A*x_cg - b;
  d_cg = -r_cg;
  del_cg = r_cg'*r_cg;
  %<errorArrayInitialization>
  errs_cg = [];
  convs_cg = [];
  conv_cg = 0;
  absErrs_cg = [];
  absErr_cg = 0;
  xAct_cg = A\b;
  xArray_cg = [norm(x_cg , 2)];
  %</errorArrayInitialization>
  while(max(abs(r_cg)) > tol)
    u_cg = A*d_cg;
    lambda_cg = del_cg / (d_cg'*u_cg);
    prevx_cg = x_cg;
    x_cg = x_cg + lambda_cg*d_cg;
    %% <errorArrays>
    err_cg = max(abs(x_cg - prevx_cg));
    errs_cg = [errs_cg ; err_cg];
    conv_cg = max(abs(x_cg - xAct_cg)) / max(abs(prevx_cg - xAct_cg));
    convs_cg = [convs_cg; conv_cg];
    absErr_cg = max(abs(x_cg - xAct_cg));
    absErrs_cg = [absErrs_cg; absErr_cg];
    xArray_cg = [xArray_cg ; norm(x_cg,2)];
    %% </errorArrays>
    prevr_cg = r_cg;
    r_cg = r_cg + lambda_cg*u_cg;
    prevdel_cg = del_cg;
    del_cg = r_cg'*r_cg;
    if count_cg > 0
      preva_cg = a_cg;
    end
    a_cg = del_cg/prevdel_cg;
    prevd_cg = d_cg;
    d_cg = -r_cg + a_cg*prevd_cg;
    count_cg = count_cg + 1;
  end
  fprintf("\nConjugate Gradient Method: done\n");
end

% Output useful for actual chemist: --------------------------------------------
fprintf("\n----------------- Extraction Results -----------------\n");
fprintf("\nMass Fraction of Chemical in Water Stream at Stage n:\n\n");
if SPD == 1
  for i = 1 : length(x_cg)
    fprintf("Stage %02d | %.22f\n", i, x_cg(i));
  end
else
  for i = 1 : length(x_gau)
    fprintf("Stage %02d | %.22f\n", i, x_gau(i));
  end
end
fprintf("\n-------------- End of Extraction Results -------------\n");

% Analysis ---------------------------------------------------------------------

userInput = 1; % Just to enter while loop
while ~(userInput==7) % start menu
  userInput = menu('Choose method to analyze:', 'Jacobi', 'SOR', 'Conjugate Gradient', 'All Iterative Methods', 'Guassian Effectiveness', 'Visualization: Conjugate Gradient', 'Quit');
  close all;
  
  % Jacobi Analysis ------------------------------------------------------------
  if userInput == 1
    figure(1)
      set(figure(1), 'Position', [100, 100, 1000, 600]);
      plot(errs_jac);
      title('Jacobi: Error v Iteration');
      ylabel('Error');
      xlabel('Iteration');
    figure(2);
      set(figure(2), 'Position', [130, 130, 1000, 600]);
      plot(absErrs_jac);
      title('Jacobi: Absolute Error v Iteration');
      ylabel('Absolute Error');
      xlabel('Iteration');
    figure(3);
      set(figure(3), 'Position', [160, 160, 1000, 600]);
      plot(convs_jac);
      title('Jacobi: Convergence Constant v Iteration');
      ylabel('Convergence Constant');
      xlabel('Iteration');
    figure(4);
      set(figure(4), 'Position', [190, 190, 1000, 600]);
      plot(xArray_jac);
      title('Jacobi: Norm of X v Iteration');
      ylabel('Norm of X');
      xlabel('Iteration');
  end

  % SOR Analysis ---------------------------------------------------------------
  if userInput == 2
    if performw == 2
      figure(1);
        set(figure(1), 'Position', [100, 100, 1000, 600]);
        plot(errs_sor);
        SOR_proper_w = (2 / (1+ sqrt(1- specRad_jac ^2)));
        str = sprintf('SOR (w = %.2f): Error v Iteration', SOR_proper_w);
        title(str);
        ylabel('Error');
        xlabel('Iteration');
      figure(2);
        set(figure(2), 'Position', [130, 130, 1000, 600]);
        plot(absErrs_sor);
        str1 = sprintf('SOR (w = %.2f): Absolute Error v Iteration', SOR_proper_w);
        title(str1);
        ylabel('Absolute Error');
        xlabel('Iteration');
      figure(3);
        set(figure(3), 'Position', [160, 160, 1000, 600]);
        plot(convs_sor);
        str2 = sprintf('SOR (w = %.2f): Convergence Constant v Iteration', SOR_proper_w);
        title(str2);
        ylabel('Convergence Constant');
        xlabel('Iteration');
      figure(4);
        set(figure(4), 'Position', [190, 190, 1000, 600]); 
        plot(xArray_sor);
        str3 = sprintf('SOR (w = %.2f): Norm of X v Iteration', SOR_proper_w);
        title(str3);
        ylabel('Norm of X');
        xlabel('Iteration'); 
    end
    if performw == 1
      figure(5);
        set(figure(5), 'Position', [220, 220, 1000, 600]);
        plot(wArr_sor, iterArr_sor);
        title('SOR: Number Iterations to Converge v Value of w');
        ylabel('Number of Iterations');
        xlabel('w');
    end 
  end

  % Conjugate Gradient Analysis ------------------------------------------------
  if userInput == 3
    if SPD == 1
      figure(1)
        set(figure(1), 'Position', [100, 100, 1000, 600]);
        plot(errs_cg);
        title('Conjugate Gradient: Error v Iteration');
        ylabel('Error');
        xlabel('Iteration');
      figure(2)
        set(figure(2), 'Position', [130, 130, 1000, 600]);
        plot(absErrs_cg);
        title('Conjugate Gradient: Absolute Error v Iteration');
        ylabel('Absolute Error');
        xlabel('Iteration');
      figure(3)
        set(figure(3), 'Position', [160, 160, 1000, 600]);
        plot(convs_cg);
        title('Conjugate Gradient: Convergence Constant v Iteration');
        ylabel('Convergence Constant');
        xlabel('Iteration');
      figure(4)
        set(figure(4), 'Position', [190, 190, 1000, 600]);
        plot(xArray_cg);
        title('Conjugate Gradient: Norm of X v Iteration');
        ylabel('Norm of X');
        xlabel('Iteration');
    else
      fprintf("Matrix must be SPD to perform Conjugate Gradient analysis.");
    end
  end

  % All 3 Method Analysis ------------------------------------------------------
if SPD == 1
  % <Errs>
    arrLen_errs = [length(errs_jac), length(errs_sor), length(errs_cg)];
    maxLen_errs = max(arrLen_errs);
    arrIt_errs = [0:1:maxLen_errs-1];
    finErrs_jac = zeros(maxLen_errs,1);
    finErrs_sor = zeros(maxLen_errs,1);
    finErrs_cg = zeros(maxLen_errs,1);
    for i = 1 : length(errs_jac)
      finErrs_jac(i) = errs_jac(i);
    end
    for i = length(errs_jac) : length(finErrs_jac)
      finErrs_jac(i) = errs_jac(length(errs_jac));
    end  % --------------------
    for i = 1 : length(errs_sor)
      finErrs_sor(i) = errs_sor(i);
    end
    for i = length(errs_sor) : length(finErrs_sor)
      finErrs_sor(i) = errs_sor(length(errs_sor));
    end % --------------------
    for i = 1 : length(errs_cg)
      finErrs_cg(i) = errs_cg(i);
    end
    for i = length(errs_cg) : length(finErrs_cg)
      finErrs_cg(i) = errs_cg(length(errs_cg));
    end % --------------------

  % </Errs>
  
  % <AbsErrs>
    arrLen_AbsErrs = [length(absErrs_jac), length(absErrs_sor), length(absErrs_cg)];
    maxLen_AbsErrs = max(arrLen_AbsErrs);
    arrIt_AbsErrs = [0:1:maxLen_AbsErrs-1];
    finAbsErrs_jac = zeros(maxLen_AbsErrs,1);
    finAbsErrs_sor = zeros(maxLen_AbsErrs,1);
    finAbsErrs_cg = zeros(maxLen_AbsErrs,1);
    for i = 1 : length(absErrs_jac)
      finAbsErrs_jac(i) = absErrs_jac(i);
    end
    for i = length(absErrs_jac) : length(finAbsErrs_jac)
      finAbsErrs_jac(i) = absErrs_jac(length(absErrs_jac));
    end  % --------------------
    for i = 1 : length(absErrs_sor)
      finAbsErrs_sor(i) = absErrs_sor(i);
    end
    for i = length(absErrs_sor) : length(finAbsErrs_sor)
      finAbsErrs_sor(i) = absErrs_sor(length(absErrs_sor));
    end % --------------------
    for i = 1 : length(absErrs_cg)
      finAbsErrs_cg(i) = absErrs_cg(i);
    end
    for i = length(absErrs_cg) : length(finAbsErrs_cg)
      finAbsErrs_cg(i) = absErrs_cg(length(absErrs_cg));
    end % --------------------

  % </AbsErrs>
  
  % <convs>
    arrLen_convs = [length(convs_jac), length(convs_sor), length(convs_cg)];
    maxLen_convs = max(arrLen_convs);
    arrIt_convs = [0:1:maxLen_convs-1];
    finconvs_jac = zeros(maxLen_convs,1);
    finconvs_sor = zeros(maxLen_convs,1);
    finconvs_cg = zeros(maxLen_convs,1);
    for i = 1 : length(convs_jac)
      finconvs_jac(i) = convs_jac(i);
    end
    for i = length(convs_jac) : length(finconvs_jac)
      finconvs_jac(i) = convs_jac(length(convs_jac));
    end  % --------------------
    for i = 1 : length(convs_sor)
      finconvs_sor(i) = convs_sor(i);
    end
    for i = length(convs_sor) : length(finconvs_sor)
      finconvs_sor(i) = convs_sor(length(convs_sor));
    end % --------------------
    for i = 1 : length(convs_cg)
      finconvs_cg(i) = convs_cg(i);
    end
    for i = length(convs_cg) : length(finconvs_cg)
      finconvs_cg(i) = convs_cg(length(convs_cg));
    end % --------------------

  % </convs>
  
  % <xArray>
    arrLen_xArray = [length(xArray_jac), length(xArray_sor), length(xArray_cg)];
    maxLen_xArray = max(arrLen_xArray);
    arrIt_xArray = [0:1:maxLen_xArray-1];
    finxArray_jac = zeros(maxLen_xArray,1);
    finxArray_sor = zeros(maxLen_xArray,1);
    finxArray_cg = zeros(maxLen_xArray,1);
    for i = 1 : length(xArray_jac)
      finxArray_jac(i) = xArray_jac(i);
    end
    for i = length(xArray_jac) : length(finxArray_jac)
      finxArray_jac(i) = xArray_jac(length(xArray_jac));
    end  % --------------------
    for i = 1 : length(xArray_sor)
      finxArray_sor(i) = xArray_sor(i);
    end
    for i = length(xArray_sor) : length(finxArray_sor)
      finxArray_sor(i) = xArray_sor(length(xArray_sor));
    end % --------------------
    for i = 1 : length(xArray_cg)
      finxArray_cg(i) = xArray_cg(i);
    end
    for i = length(xArray_cg) : length(finxArray_cg)
      finxArray_cg(i) = xArray_cg(length(xArray_cg));
    end % --------------------
end
  % </xArray>
  
  if userInput == 4
    if (performw == 2) && (SPD == 1)
      figure(1)
        hold on
        set(figure(1), 'Position', [100, 100, 1000, 600]);
        plot(arrIt_errs, finErrs_jac, 'r', arrIt_errs, finErrs_sor, 'b', arrIt_errs, finErrs_cg, 'g')
        title('All Methods: Error v Iteration');
        ylabel('Error');
        xlabel('Iteration');
        legend('Jacobian', 'SOR', 'Conjugate Gradient')
        legend('Location','NorthWest')
        hold off;
      figure(2)
        hold on
        set(figure(2), 'Position', [130, 130, 1000, 600]);
        plot(arrIt_AbsErrs, finAbsErrs_jac, 'r', arrIt_AbsErrs, finAbsErrs_sor, 'b', arrIt_AbsErrs, finAbsErrs_cg, 'g')
        title('All Methods: Absolute Error v Iteration');
        ylabel('Absolute Error');
        xlabel('Iteration');
        legend('Jacobian', 'SOR', 'Conjugate Gradient')
        legend('Location','NorthWest')
        hold off;
      figure(3)
        hold on
        set(figure(3), 'Position', [160, 160, 1000, 600]);
        plot(arrIt_convs, finconvs_jac, 'r', arrIt_convs, finconvs_sor, 'b', arrIt_convs, finconvs_cg, 'g')
        title('All Methods: Convergence Constant v Iteration');
        ylabel('Convergence Constant');
        xlabel('Iteration');
        legend('Jacobian', 'SOR', 'Conjugate Gradient')
        legend('Location','NorthWest')
        hold off;
      figure(4)
        hold on
        set(figure(4), 'Position', [190, 190, 1000, 600]);
        plot(arrIt_xArray, finxArray_jac, 'r', arrIt_xArray, finxArray_sor, 'b', arrIt_xArray, finxArray_cg, 'g')
        title('All Methods: Norm of X v Iteration');
        ylabel('Norm of X');
        xlabel('Iteration');
        legend('Jacobian', 'SOR', 'Conjugate Gradient')
        legend('Location','NorthWest')
        hold off;
    else
      fprintf("\nMust choose to not perform multple w in beginning, or,  matrix must be SPD.\n");
    end
  end
  
  if userInput == 5
    fprintf("\nGaussian Elimiation w/ Scaled Partial Pivoting Output:\n");
    disp(x_gau);
    absErr_gau = max(abs((abs(x_true)-abs(x_gau))));
    fprintf("\nAbsolute error between Gaussian output and Matlab's solution: %.22f\n", absErr_gau);  
    fprintf("Note: k(A) = %.2f; error is accurate to 22 decimals places\n", k_k);
  end
  
  if userInput == 6
    [x_vis,y_vis] = meshgrid([-2:.2:2]);
    Z_vis = (x_vis.^2 + y_vis.^2)/2;
    figure(1)
    surf(x_vis,y_vis,Z_vis,gradient(Z_vis))
    title('Visualization of ConjGrad Functional | x^T = [x y] A = I(2x2)');
    colorbar
  end
  
  if userInput == 7
   fprintf("\n----------------- Program Terminated -----------------\n");
  end
  
end % end menu

% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program under the terms of 
% the GNU Lesser General Public License, either version 3 of the License, 
% or any later version.

% easy_runtest

% Set parameters if such parameters don't exist in work space
disp('Esay run test for Inverse power method program')
if ~exist('n','var');     n = 100;                            end          
if ~exist('m','var');        m= 500;                          end
if ~exist('mu','var');        mu = n/2;                       end
if ~exist('e','var');        e = 0.01;                        end
if ~exist('display','var');        display = 3;               end
if ~exist('seed','var');        seed = 3;                     end
more off

fprintf('Opening file "Report.txt"... ');
fid = fopen('Report.txt','w');
disp('done!')

report_ipm = "Report for easy_runtest for Inverse Power Method\nGenerated a random symmetric matrix by A = a'*a while a = randi(n,n), \nperform IPM on it, compare the result with in built linear solvers and \nsee the relationship between iterations and ||mu - dominant||/||mu - subdominant|| ";
% display the parameters
report_ipm = [report_ipm, sprintf('\nParameters:\n\n')];
report_ipm = [report_ipm, sprintf('n = %d\n',n)];
report_ipm = [report_ipm, sprintf('m = %d\n',m)];
report_ipm = [report_ipm, sprintf('mu = %d\n',mu)];
report_ipm = [report_ipm, sprintf('e = %f\n',e)];
report_ipm = [report_ipm, sprintf('seed = %f\n',seed)];
report_ipm = [report_ipm, sprintf('display = %d\n',display)];
report_ipm = [report_ipm, sprintf('You can check initial vector x0 in workspace.\n\n')];
report_ipm = [report_ipm, sprintf('n is the size of the matrix.\n')];
report_ipm = [report_ipm, sprintf('You can change the value of seed in the random number generator.\n')];
report_ipm = [report_ipm, sprintf('m is the maximum iteration rounds.\n')];
report_ipm = [report_ipm, sprintf('mu is y = (A - mu*I)^(-1) * x in each iteration.\n')];
report_ipm = [report_ipm, sprintf('e is the allowed error\nERR = || |x_(k+1) - x_(k)| ||_infty in each iteration\n\n')];
report_ipm = [report_ipm, sprintf('Parameters can be changed as variables in the workspace\n')];

rand('seed',seed)

fprintf('Generating a random symetric n ny n matrix...')
% A = rand(n,n);
a = rand(n,n);
A = a'*a;
x0 = (1:n)';
fprintf('done!\n')
    
% Compare the other eigenvalues around it
fprintf('Using eig function to find the eigenvalues...')
eig_list = eig(A);
distances = abs(eig_list - mu);
eig_disp = zeros(1,display);
fprintf('done!\n')


for j = 1 : display
    [~,index] = min(distances);
    distances(index) = [];
    eig_disp(j) = eig_list(index);
    eig_list(index) = [];
end

fprintf('Doing Inverse power method based on given mu and x0...')
if i < m
    [eig_value,~,i] = IPM(A,mu,x0,m,e);
    fprintf('done!\n')
    report_ipm = [report_ipm, sprintf('The eigenvalue nearest %.2f is:\n',mu)];
    report_ipm = [report_ipm, num2str(eig_value),"\n"];
    report_ipm = [report_ipm, sprintf('Using eig() function, the nearest %d eigenvalues to %.1f are:\n',display,mu)];
    for eigen = eig_disp
        report_ipm = [report_ipm, num2str(eigen),' '];
    end
    report_ipm = [report_ipm,"\n"];
    report_ipm = [report_ipm, sprintf('The iteration rounds for the nearest eigenvalue is:\n')];
    report_ipm = [report_ipm, num2str(i),"\n"];
    report_ipm = [report_ipm, sprintf('||mu - dominant||/||mu - subdominant|| = ')];
    report_ipm = [report_ipm, num2str(abs(mu-eig_disp(1))/abs(mu-eig_disp(2)))];
else
    report_ipm = [report_ipm, sprintf('Using eig() function, the nearest %d eigenvalues to %.1f are:\n',display,mu)];
    for eigen = eig_disp
        report_ipm = [report_ipm, num2str(eigen),' '];
    end
    report_ipm = [report_ipm,"\n"];
    report_ipm = [report_ipm, sprintf('||mu - dominant||/||mu - subdominant|| = ')];
    report_ipm = [report_ipm, num2str(abs(mu-eig_disp(1))/abs(mu-eig_disp(2))),"\n"];
    report_ipm = [report_ipm, sprintf('\nMaximum iteration reached, exiting the program.\n')];
end


for text = report_ipm
    fprintf(fid,text);
end
fprintf('Closing file "Report.txt"... ')
fclose(fid);
fprintf('done!\n')
disp('Bye!')
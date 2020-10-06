% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program under the terms of 
% the GNU Lesser General Public License, either version 3 of the License, 
% or any later version.

% runtest 

disp('This test finds the few largest eigenvalue of the symmetric square matrices.')
if ~exist('n','var');     n = 30;                             end          
if ~exist('m','var');        m= 800;                          end
if ~exist('e','var');        e = 0.005;                       end
if ~exist('e1','var');        e1 = 0.1;                       end
if ~exist('e2','var');        e2 = 0.07;                       end
if ~exist('e3','var');        e3 = 0.03;                      end
if ~exist('seed','var');        seed = 3;                     end
if ~exist('findn','var');     findn= 3;                       end
if ~exist('repeat','var');     repeat= 3;                     end
more off

fprintf('Opening file "Report_1.txt"... ');
fid = fopen('Report_1.txt','w');
disp('done!')

report = sprintf("Report for runtest for Inverse Power Method\nGenerated a random symmetric matrix by A = a'*a while a = randi(n,n), \nuse Inverse Power Method to find %d largest eigenvalues.",findn);
report = [report, sprintf('\nProcess:\n')];
report = [report, sprintf('1) Find the largest eigenvalue using Power Method.\n')];
report = [report, sprintf('2) Find the smallest eigenvalue using Inverse Power Method.\n')];
report = [report, sprintf('3) Choose mu between the smallest to largest, the step \nsizes are equally divided in log scale.\n')];

% display the parameters
report = [report, sprintf('\n=========================\n')];
report = [report, sprintf('\nParameters:\n\n')];
report = [report, sprintf('n = %d\n',n)];
report = [report, sprintf('m = %d\n',m)];
report = [report, sprintf('e = %f\n',e)];
report = [report, sprintf('e1 = %f\n',e1)];
report = [report, sprintf('e2 = %f\n',e2)];
report = [report, sprintf('e3 = %f\n',e3)];
report = [report, sprintf('seed = %f\n',seed)];
report = [report, sprintf('findn = %d\n',findn)];
report = [report, sprintf('repeat = %f\n',repeat)];
report = [report, sprintf('You can check initial vector x0 in workspace.\n\n')];
report = [report, sprintf('n is the size of the matrix.\n')];
report = [report, sprintf('m is the maximum iteration rounds.\n')];
report = [report, sprintf('repeat is the searching rounds executed in the second round.\n')];
report = [report, sprintf('mu is y = (A - mu*I)^(-1) * x in each iteration.\n')];
report = [report, sprintf('Change seed to 0 if you want the program to generate the random matrix.\n')];
report = [report, sprintf('e is the allowed error\nERR = || |x_(k+1) - x_(k)| ||_infty in each iteration.\n')];
report = [report, sprintf('The program considered |lambda1-lambda2|/|lambda1|>e1 as two distinct eigenvalues.\n')];
report = [report, sprintf('Parameters can be changed as variables in the workspace.\n\n')];


rand('twister',seed)
fprintf('Generating a random symmetric n by n matrix...')
% A = rand(n,n);
a = rand(n,n);
A = a'*a;
x0 = (1:n)';
fprintf('done!\n')

fprintf('Using Power Method to find the largest eigenvalue...')
largest_eig = PM(A,A(:,1),e,m);
fprintf('done!\n');

fprintf('Using Inverse Power Method to find the smallest eigenvalue...')
mu_1 = 0;
smallest_eig = IPM(A, mu_1,A(:,1),m,e);
fprintf('done!\n')

fprintf('Using in-built eig function to find the %d largest eigenvalues...', findn);
eig_list = flip(eig(A));
disp_list = eig_list(1:findn)';
report = [report, sprintf('Using in-built eig function, the largest %d eigenvalues are:\n',findn)];
for item = disp_list
    report = [report,num2str(item)," "];
end
report = [report, "\n"];
fprintf('done!\n')

fprintf('Finding %d largest eigenvalues using IPM ...', findn)
step_len = (2*largest_eig/n-smallest_eig)/(n);
eigenlist = [largest_eig];
j = 0;
mu_2 = 2*largest_eig/n;
eigen = largest_eig;

while length(eigenlist) < findn 
    mu_2 = mu_2 - step_len;
    [eigen,~,iterations] = IPM(A, mu_2, A(:,1), m, e);

    if eigen == 0 || iterations == 0 || isnan(eigen) 
        continue
    end 
    
    if eigen < smallest_eig
        continue
    end
    
    if iterations == m
        if j == 0
            fprintf('\nHere is a list of points which fail to converge with maximum iteration \nnumber %d in IPM:\n',m);
        end
        
        if mod(j,5) == 0 && j ~= 0
            fprintf("\n");
        end
        
        fprintf('%8d   ',mu_2)
        j = j+1;
        continue
    end
    
    if ~isnan(eigen) && ~eigen == 0 
        dist_list = abs((eigenlist - eigen)./eigen); 
    else
        continue
    end
    
    if all(dist_list>e1)
        eigenlist = [eigenlist,eigen];
    else
        continue
    end
    
end
report = [report, "\n"];
fprintf('\ndone!\n')
report = [report, sprintf('The %d largest eigenvalues using Inverse Power Method: \n',findn)];
fprintf('Checking if there are other eigenvalues between those %d eigenvalues...',findn);

for j = 1 : repeat
    eigenlist = flip(sort(eigenlist));
    for itval = 1: length(eigenlist)-1
        steplen = eigenlist(itval)-eigenlist(itval+1);
        
        if itval<3
            tol = e2 * steplen;
        else 
            tol = e3 * steplen;
        end
        
        miditval = 0.3*steplen;
        leftp = eigenlist(itval)-miditval;
        leftconv = IPM(A,leftp,A(:,1),m,e);
        rightp = eigenlist(itval+1) + miditval;
        rightconv = IPM(A,rightp,A(:,1),m,e);
        if rightconv-eigenlist(itval+1)>tol
            eigenlist = [eigenlist, rightconv];
        elseif eigenlist(itval)-leftconv > tol
            eigenlist = [eigenlist,leftconv];
        end
    end
    eigenlist = flip(sort(eigenlist));
    
    if length(eigenlist) < findn 
      continue
    else
      eigenlist = eigenlist(1:findn);
    end
end
fprintf('done!\n')

for display = eigenlist(1:length(eigenlist))
      report = [report, num2str(display),' '];
end

report = [report, "\n"];
report = [report,"The process ended!\n"];

fprintf('Writing into the file "Report_1.txt", closing the file...')
for line=report(1:end)
    fprintf(fid,line);
end
fclose(fid);
fprintf('done!\n')
fprintf('Bye!\n')
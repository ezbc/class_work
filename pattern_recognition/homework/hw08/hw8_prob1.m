clear;
close all;

a = [-2; -1; 1; 10];

feet = [5; 5; 6; 6;];
inches = [10; 11; 1; 10;];
a = feet + inches/12.0;

b = [-1; -1; 1; 1];
xhat_ls = (a'*a)\a'*b;
disp('Least squares xhat');
disp(xhat_ls)
disp(a*xhat_ls)

%svmStruct = svmtrain(heights,Basketball);

lambda = 1;
minVal = 1e+6;
grid_res = 0.1;
alpha_start = -1;
alpha_end = 1;
t0 = 1;

for ii = alpha_start:grid_res:alpha_end
    for jj = alpha_start:grid_res:alpha_end
        for kk = alpha_start:grid_res:alpha_end
            for ll = alpha_start:grid_res:alpha_end
                alpha = [ii; jj; kk; ll];
                totalSum = 0;
                for i=1:4
                    innerSum = 0;
                    innerSum2 = 0;
                    for j=1:4
                        innerSum = innerSum + alpha(j)*a(i)'*a(j);
                        innerSum2 = alpha(i)*alpha(j)*a(i)'*a(j);
                    end
                    totalSum = totalSum + (t0 -innerSum) + lambda*innerSum2;
                    if abs(totalSum) < minVal
                        alpha_min = alpha;
                        minVal = abs(totalSum);
                    end
                end
            end
        end
    end
end

for j = 1:4
    xhat = sum(alpha_min.*a);
end

disp('xhat');
disp(xhat);
disp('alpha');
disp(alpha_min);



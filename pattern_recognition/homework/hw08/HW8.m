a = [-2; -1; 1; 10];
b = [-1; -1; 1; 1];
xhat = (a'*a)\a'*b;

%svmStruct = svmtrain(heights,Basketball);

lambda = 1;
minVal = 1e+6;
grid_res = 0.01;
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
                        xhat = sum(alpha(j) * a(j));
                        alpha_min = alpha;
                        minVal = abs(totalSum);
                    end
                end
            end
        end
    end
end

disp('x');
disp(x);
disp('alpha');
disp(alpha);







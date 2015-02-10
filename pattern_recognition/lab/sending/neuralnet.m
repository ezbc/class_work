x=A';
t=zeros(2,70);
for i=1:length(b)
    if b(i,1)==1
        t(1,i)=1;
        t(2,i)=0;
    else
        t(1,i)=0;
        t(2,i)=1;
    end
end

setdemorandstream(391418381)

net = patternnet(30);
view(net)

[net,tr] = train(net,x,t);
nntraintool

plotperform(tr)

testX = x(:,tr.testInd);
testT = t(:,tr.testInd);

testY = net(testX);
testIndices = vec2ind(testY)

plotconfusion(testT,testY)

[c,cm] = confusion(testT,testY)

fprintf('Percentage Correct Classification   : %f%%\n', 100*(1-c));
fprintf('Percentage Incorrect Classification : %f%%\n', 100*c);

plotroc(testT,testY)






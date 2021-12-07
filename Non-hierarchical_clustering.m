fid = fopen('data2.txt', 'r');

X0 = fscanf(fid, '%g %g ', [2 inf]);  
X = X0';
fclose(fid);

figure;
plot(X(:,1), X(:,2), '*');


numClust = 2;
C = zeros(numClust,2); 

for i = 1 : numClust
    C(i, :) = X(i , :);
end

%???????????? ?????????? ?? ???????
U = zeros(length(X) , 2);
eps = 1e-6;
Q = 1000;
q = 0;

R = zeros(1,numClust);
for i = 1 : length(X) 
        for n = 1 : numClust
        R(n) = pdist([X(i , :) ; C(n , :)],'euclidean');
        end 
    for n = 1 : numClust
        if R(n) == min(R)
            U(i,1) = n;
            U(i,2) = R(n);
        end
    end
end

while abs(Q - q)>eps
    Q = q    
    C = zeros(numClust,2);
    for l = 1 : numClust 
        n = 0;
        for j = 1 : length(U)
            if (U(j,1)==l)
                n = n+1;
                C(l, :) = C(l, :) + X(j , :);
            end
        end
        C(l, :) = C(l, :)./n; 
    end
    
    P = zeros(1,numClust);
    for i = 1 : length(X) 
        for n = 1 : numClust
        P(n) = pdist([X(i , :) ; C(n , :)],'euclidean');
        end 

        for n = 1 : numClust
            if P(n) == min(P)
                U(i,1) = n;
                U(i,2) = P(n);
            end
        end
    end 
    
    u = zeros(1,numClust);
    for l = 1 : numClust  
        for j = 1 : length(U)
            if (U(j,1)==l)
                u(l) = u(l)+U(j,2);
            end
        end
    end
    q = sum(u); 
  
end

figure
gscatter(X(:,1), X(:,2), U(:,1))
hold on
scatter(C(:,1),C(:,2),'*')    
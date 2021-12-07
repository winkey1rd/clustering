
fid = fopen('data1.txt', 'r');

X0 = fscanf(fid, '%g %g', [2 inf]); %25x2 
fclose(fid);

X=X0';
figure
plot(X(:,1),X(:,2), '*');
title('Input data');

qualityFactor = zeros(3);

%??????? ??????????: 1 ? ?????????, 2 ? ??????????????????? ?????????, 3 ? ??????
%                       'euclidean'                 'seuclidean'            'cityblock'
%?????? ??????????: a ? ???????? ??????, b ? ???????? ??????, c ? ??????? ?????
%                       'single'                'complete'               'average'

%?????????
DEuclidean = pdist(X, 'euclidean'); 
EuclideanSingle = linkage(DEuclidean,'single');
qualityFactor(1,1) = cophenet(EuclideanSingle, DEuclidean);
EuclideanComplete = linkage(DEuclidean,'complete'); 
qualityFactor(1,2) = cophenet(EuclideanComplete, DEuclidean);
EuclideanAverage = linkage(DEuclidean,'average');
qualityFactor(1,3) = cophenet(EuclideanAverage, DEuclidean);

%??????????????????? ?????????
Dseuclidean = pdist(X, 'seuclidean');
SeuclideanSingle = linkage(Dseuclidean,'single');
qualityFactor(2,1) = cophenet(SeuclideanSingle,Dseuclidean);
SeuclideanComplete = linkage(Dseuclidean,'complete');
qualityFactor(2,2) = cophenet(SeuclideanComplete,Dseuclidean);
SeuclideanAverage = linkage(Dseuclidean,'average');
qualityFactor(2,3) = cophenet(SeuclideanAverage,Dseuclidean);


%??????
Dcityblock = pdist(X, 'cityblock');
CityblockSingle = linkage(Dcityblock,'single');
qualityFactor(3,1) = cophenet(CityblockSingle,Dcityblock);
CityblockComplete = linkage(Dcityblock,'complete');
qualityFactor(3,2) = cophenet(CityblockComplete,Dcityblock);
CityblockAverage = linkage(Dcityblock,'average');
qualityFactor(3,3) = cophenet(CityblockAverage,Dcityblock);

%?????????????
k = qualityFactor'

figure
dendrogram(EuclideanAverage);
title('MostEffective - EuclideanAverage');
T = cluster(EuclideanAverage, 'maxclust',4);

figure
gscatter(X(:,1),X(:,2),T);
hold on
k1=1; k2=1; k3=1; k4=1;

%???? ?? ??????? ????????
for i=1:24;
     if (T(i) == 1)
        clust1(k1,:) = X(i,:);
        k1 = k1+1;
     end
     if (T(i) == 2)
        clust2(k2,:) = X(i,:);
        k2 = k2+1;
     end
     if (T(i) == 3)
        clust3(k3,:) = X(i,:);
        k3 = k3+1;
     end
     if (T(i) == 4)
        clust4(k4,:) = X(i,:);
        k4 = k4+1;
     end
end

%??????
M = [mean(clust1);mean(clust2);mean(clust3);mean(clust4)];
VAR = [var(clust1);var(clust2);var(clust3);var(clust4)];

DCenters = pdist(M,'euclidean');

Dclust1 = squareform(pdist([clust1;mean(clust1)],'euclidean'));
Dclust1_ = Dclust1(end,:);

Dclust2 = squareform(pdist([clust2;mean(clust2)],'euclidean'));
Dclust2_ = Dclust2(end,:);

Dclust3 = squareform(pdist([clust3;mean(clust3)],'euclidean'));
Dclust3_ = Dclust3(end,:);

Dclust4 = squareform(pdist([clust4;mean(clust4)],'euclidean'));
Dclust4_ = Dclust4(end,:);

plot(M(:,1),M(:,2), '+black');
hold off
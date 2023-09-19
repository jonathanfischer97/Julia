%Read in a dataset with D datapoints and M features
%Our M features are the rate constants and concentrations that are
%optimally selected to produce oscillations, with some of them fixed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Currently it is hard-coded in feature space:
%Column 1: fitness, Col2: Period, Col3: Amplitude
%Col20: 'A'
%
%Use the log of the data, as we are sampling across several orders of 
%magnitude. In this way, we do not bias against small rates in favor of
%large rates. 
%When using the log(data), do not divide out means during normalization, 
% as the exponents vary around zero.
%
%Evaluate the principal components to identify different classes of
%solutions
 

%return the 
%Number of clusters
%total size of solutions based on summing over clusters
%average period and amplitude of each cluster
%The top variables that generate variance between clusters
%the top variables that generate variance between data (top ones are
%similar across both).
%the entire sorted variance vector along with the names of the variables,
%so we can see the top and bottom ones. Ignore the fixed variables. 
%Return the variance as well of the fitness, period and amplitude.

%numCluster, avgDist, sumD, avgAmp, avgPeriod, sortSdev, idsSdev,
%sortData, idData, labels

%distances between clusters? 
function[varM,  cluster, dataStats]=find_solution_classes_analyze(filename)

%hard coded column numbers. If they change, update these variables.
colPeriod=2; %column with period
colAmplitude=3; %column with amplitude
colA=20; %column with A concentration.
colStart=4; %1,2,3 are the fitness,period, amplitude, so parameters start in column 4
%Load in data that is in a numeric array with a header line.
M=importdata(filename);


%names of variables in order.
M.colheaders 

%data in a numeric array
M.data; %M is of size nxm.
[N, Nfeature]=size(M.data);%N is the number of data points.

labels=M.textdata(colStart:end);
%U is nxn, S is nxm and V is mxm the right singular vectors. 
%[U, S, V]=svd(M.data); 
%S are singular values. 
% M= U*S*V'
%sigma=diag(S);

%sigma=sqrt(lambda) where lambda are the eigenvales of M'M

% renormalize="FALSE"
% %divide out the standard deviation
% if(renormalize=="TRUE")
%     sdev=std(M.data);
%     tmp=find(sdev<1e-10);
%     sdev(tmp)=1; %replace 0s with ones. 
%     varM=M.data./sdev;%if std=0, keep at fixed value.e
% else
%     varM=M.data
%     display('NOT renormalizing by standard deviation of each variable!')
% end

display('Taking Log of data')
%don't divide it out!
varM=log10(M.data)

%Skip the first 3 columns, fitness, period and amplitude
varM=varM(:,colStart:end);

if(N==1)
    display('ONLY ONE SOLUTION')
    dataStats.sortSdev=zeros(1, Nfeature);
    cluster.N=1;
    cluster.sizes=1;
    cluster.kmeans(1,:)=M.data(1,colStart:end);
   
    cluster.avgAmp=M.data(colAmplitude)/M.data(colA)*100;
    cluster.avgPeriod=M.data(colPeriod);
    return
end
%Coefficients are the eigenvectors of the covariance matrix (M'*M). 
%score are the data projected onto the eigenvectors, so the first 2 columns
%are all data points projected onto the first 2 eigenvectors.
%latent are the eigenvalues.
[cof, score, latent]=pca(varM);


%calculate the standard deviation in each feature of the data. 
[sortData, idData]=sort(std(varM),'descend')

%Report statistics on the raw dataset after the log transform
dataStats.sortSdev=sortData
dataStats.ids=idData
dataStats.sortlabels=labels(idData) %requires that the labels are the same size as varM. 


%Broad estimate of the cluster sizes, just chose a cutoff in the relative eigenvalues 
fraction=latent./latent(1);
cutOff=0.1
numCluster=length(find(fraction>cutOff));

%L=2;
%firstL=cof(:,1:L);
%transformData=varM*firstL; %this is the same as the columns of score except they have the mean subtracted off.


%PLOT the data onto the first two principal components (PCs)
figure(2)
cmap=colormap('jet');
plot(score(:,1),score(:,2),'o','Color',cmap(1,:))


firstNC=cof(:,1:5);%Calculate for the first 5 because we need the means.
initPCs=varM*firstNC; %these are the data projected onto the firstL PCs, without means subtracted off.
mu=mean(initPCs); %calculate the means of the data onto the initial PCS so it can be subtracted off


%N is the number of data points found.
if(N<256)
    Ncolor=4;
    Nblock=floor(N/Ncolor);
else
    Ncolor=256;
    Nblock=floor(N/256);
end

%PLOT the data onto the first 3 PCs, color by appearance in the array.
f=figure(3);
ax=axes('Parent',f,'FontSize',20,'LineWidth',1);%,'XScale','linear','YScale','linear');
hold(ax);
title('PC1, 2, 3. Colored: first to last')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')

for n=1:1:Ncolor
    start=(n-1)*Nblock+1;
    fin=start+Nblock;
    if(n==Ncolor)
        fin=N;
    end
    plot3(initPCs(start:fin,1),initPCs(start:fin,2),initPCs(start:fin,3),'o','Color',cmap(n,:))
end


%Plot the data onto PCs 1,4,5
f=figure(4);
ax=axes('Parent',f,'FontSize',20,'LineWidth',1);%,'XScale','linear','YScale','linear');
hold(ax);
title('PC1, 4,5. Colored: first to last')
xlabel('PC1')
ylabel('PC4')
zlabel('PC5')
for n=1:1:Ncolor
    start=(n-1)*Nblock+1;
    fin=start+Nblock;
    if(n==Ncolor)
        fin=N;
    end
    plot3(score(start:fin,1),score(start:fin,4),score(start:fin,5),'*','Color',cmap(n,:))
end

%%%%%%%%%%%%%%%%%%%%
%Cluster the data
nClus=numCluster; %we got numCluster from the PCA using a rough cutoff on eigenvalues.
[idx, clusMeans, sumD]=kmeans(varM, nClus);

for i=1:nClus
    numC(i)=length(find(idx==i));
end
%fano=std(C)./mean(C);

%[sortFano,idsFano]=sort(fano,'descend')
[sortSdev,idsSdev]=sort(std(clusMeans),'descend')%When using log(data), do not divide by mean,

M.colheaders{idsSdev}

%Plot along the first 3 PCs, color by kmeans cluster. 
f=figure(5);
ax=axes('Parent',f,'FontSize',20,'LineWidth',1);%,'XScale','linear','YScale','linear');
hold(ax);
title('PCs 1,2,3. Color by kmeans Cluster')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
clusD=zeros(nClus-1,1);
for i=1:nClus
    ids=find(idx==i);
    color=(i-1)*floor(256/nClus)+1;
    %Plot the data
    plot3(score(ids,1),score(ids,2),score(ids,3),'s','Color',cmap(color,:))
    %Plot the cluster means
    plot3(clusMeans(i,:)*cof(:,1)-mu(1), clusMeans(i,:)*cof(:,2)-mu(2), clusMeans(i,:)*cof(:,3)-mu(3),'k+','MarkerSize',8)
        
end

%find the squared distance between means, use the largest separation as a metric
maxClusterSep=0;
for i=1:nClus
    for j=i+1:nClus
        df=clusMeans(i,:)-clusMeans(j,:);
        dist2=df*df';
        if(dist2>maxClusterSep)
            maxClusterSep=dist2;
        end
    end
end

%Plot along the first 3 PCs, color by kmeans cluster. 
f=figure(6);
ax=axes('Parent',f,'FontSize',20,'LineWidth',1);%,'XScale','linear','YScale','linear');
hold(ax);
title('PC 1,4,5: Color by kmeans clusters')
xlabel('PC1')
ylabel('PC4')
zlabel('PC5')
for i=1:nClus
    ids=find(idx==i);
    color=(i-1)*floor(256/nClus)+1;
    %plot the data
    plot3(score(ids,1),score(ids,4),score(ids,5),'s','Color',cmap(color,:))
    %plot the cluster means
    plot3(clusMeans(i,:)*cof(:,1)-mu(1), clusMeans(i,:)*cof(:,4)-mu(4), clusMeans(i,:)*cof(:,5)-mu(5),'k+','MarkerSize',8)
end

%what are the sizes of the clusters.
avgDist=zeros(nClus,1);
avgPeriod=zeros(nClus,1);
avgAmp=zeros(nClus,1);
stdPeriod=zeros(nClus,1);
stdAmp=zeros(nClus,1);
pointsPerClus=zeros(nClus,1);
relAmp=M.data(:,colAmplitude)./M.data(:,colA)*100; %in units of Percent, max 100%
period=M.data(:,colPeriod);
for i=1:nClus
    ids=find(idx==i);
    %matlab automatically calculates sum distance in the kmeans function call.
    %to get average distance we can just divide the sum by the number of
    %points.
    %dtoC=varM(ids,:)-clusMeans(i,:);
    len=length(ids);
    avgAmp(i)=mean(relAmp(ids));
    avgPeriod(i)=mean(period(ids));
    stdAmp(i)=std(relAmp(ids));
    stdPeriod(i)=std(period(ids));
    pointsPerClus(i)=len;
    %for j=1:len
    %    avgDist(i)=avgDist(i)+dtoC(j,:)*dtoC(j,:)';
    %end
    %avgDist(i)=avgDist(i)/len;
end

%Plot the data along the first 2 PCs, and color the points according to
%period and amplitude. Do they correlate with clusters? 
%period is colPeriod, amplitude is colAmplitude.
  %Rel amplitude is colAmplitude/colA 

delAmp=100/255; %use 256 colors. 
figure(10)
hold on
numError=0;
for i=1:1:N
    color=floor(relAmp(i)/delAmp)+1;
    if(color>256)color=256;
        display('WARNING, REL AMP> 100%')
        relAmp(i)
        numError=numError+1;
    end
    plot3(score(i,1),score(i,2), score(i,3),'o','Color',cmap(color,:))
    
end
for i=1:nClus
    plot3(clusMeans(i,:)*cof(:,1)-mu(1), clusMeans(i,:)*cof(:,2)-mu(2), clusMeans(i,:)*cof(:,3)-mu(3),'k+','MarkerSize',8)
end
title('PCs 1,2,3. Amplitude is color')

%Plot now coloring by period.
figure(11)
hold on
delPeriod=max(period)/255;
for i=1:1:N
    color=floor(period(i)/delPeriod)+1;
    plot3(score(i,1),score(i,2),score(i,3),'o','Color',cmap(color,:))
   
end
for i=1:nClus
    plot3(clusMeans(i,:)*cof(:,1)-mu(1), clusMeans(i,:)*cof(:,2)-mu(2), clusMeans(i,:)*cof(:,3)-mu(3),'k+','MarkerSize',8)
end
title('PCs 1,2,3. Period is color')


%Store properties of the clusters
cluster.N=numCluster
cluster.sizes=pointsPerClus
cluster.avgDist2=sumD./pointsPerClus;%this is the average distance from points to the cluster centroids, in distance^2 units. 
%cluster.sumDist=sumD
cluster.avgAmp=avgAmp
cluster.avgPeriod=avgPeriod
cluster.stdAmp=stdAmp
cluster.stdPeriod=stdPeriod
cluster.sortedStd=sortSdev
cluster.ids=idsSdev
cluster.sortLabels=labels(idsSdev)
cluster.latent=latent
%calculate the cluster means in the same sorted order
cmeans=clusMeans(:,idsSdev);
%Report clusters back in the same units as data, so not on the log scale
cluster.kmeans=10.^cmeans(:,1:end-4); %ignore the last 4 values, they are fixed.
cluster.numError=numError;
cluster.maxClusterSep=maxClusterSep;





function [sumhog,answer]=dp_hog_function(n,hogarray,m,newindices,exp_thick)

carraycopy=hogarray;
carraycopy(find(hogarray==0))=Inf;
minarray=min(carraycopy,[],2); 
maxarray=max(hogarray,[],2);
norm_array=(carraycopy-repmat(minarray,[1,size(hogarray,2)]))./repmat((maxarray-minarray),[1,size(hogarray,2)]);
score=norm_array;
scale=exp_thick/25;
threshold=0.2;
minscale=scale*(1-threshold);
maxscale=scale*(1+threshold);
%minimize hog difference
f=ones(m,n)*(inf);
%step=(endslice-startslice)/220;%2.35;
%f i j is the summation of the first i slices to the first j atlas slices
for i=1:m
	for j=1:n
		if(i==1)
			f(i,:)=score(i,:);
			sol{i,j}=[j];%best solution is just j for each  i j ...f[i,j]gives the distance matching slice 1 to j.
		else
			for k=1:j-1%if there are more than one experimental slice, for atlas slices 1 to j-1
				distance=newindices(1,i)-newindices(1,i-1);
				%score(i,j)+f(i-1,k)+Inf*((j-k)<minscale*distance || (j-k)>maxscale*distance)
				%sol{i-1,k}
				if  (j-k>=floor(minscale*distance)) && (j-k<=ceil(maxscale*distance)) && (score(i,j)+f(i-1,k)<f(i,j))

					f(i,j)=score(i,j)+f(i-1,k);
					%f(i,j)=score(i,j)+f(i-1,k)-abs((ori(i)-ori(i-1))*step-(j-k))*alpha;
					sol{i,j}=[sol{i-1,k} j];
				end
			end
		end
	end
end
%f is the distance from m slices in experimental slice to n slices in the atlas

[~,best]=min(f(m,:));
% m
% best
% size(sol)
answer=transpose(sol{m,best});
selectedhog=zeros(1,m);
%sum the original hog instead of the normalized hog
for i=1:m
    selectedhog(1,i)=hogarray(i,answer(i));
end
sumhog=sum(selectedhog);

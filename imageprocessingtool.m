function [image] = imageprocessingtool(image,n1,med1);


paddedA = padarray(image, [round(n1-1)/2 round(n1-1)/2],"replicate",'both');

% imtool(paddedA)

%Array Size
[M, N] = size(image);

NextStep = zeros(size(image));

for i=1:M
    for j=1:N
        A=paddedA(i:n1+i-1,j:n1+j-1);

        flatA=A(:);

        sortedA=sort(flatA,'descend');

        Myarray = sortedA(1:med1);

        Myvalue = median(Myarray);

        if image(i,j)>Myvalue
            NextStep(i,j) = image(i,j);
        else
            NextStep(i,j) = Myvalue;
        end

    end
end
NextStep=uint8(NextStep);
image=NextStep;
% % figure(2);imshow(NextStep)
% 
% PaddedNextStep= padarray(NextStep, [round(n1-1)/2 round(n1-1)/2],"replicate",'both');
% 
% for i=1:M
%     for j=1:N
%         A=PaddedNextStep(i:n1+i-1,j:n1+j-1);
%         NextStep(i,j) = mean(mean(A));
% 
%     end
% end
% 
% NextStep=uint8(NextStep);
% image=NextStep;
% figure(3);imshow(NextStep)

end
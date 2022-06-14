%% getCentroids_v1_1
%  Version 1.1
%  Author: Adeyinka Lesi
%  Date: 4/6/17

function [cents] = getCentroids_v1_1(img)
% Find the centroids of the objects in the image; individual objects are
% identified with a unique number (so the object doesn't have to be
% connected)

list = nonzeros(unique(img))';
if(isempty(list))
    cents = struct([]);
else
    cents(length(list)).Centroid = zeros(1,2);
    cents(length(list)).MedianX = zeros(1,2);
    cents(length(list)).MedianY = zeros(1,2);
end
for oid = list
    [y,x] = find(img==oid);
    cents(oid).Centroid = [mean(x),mean(y)];
    indlist = find(img==oid);
    sel = max(1,round(length(indlist)*0.5));
    medind = indlist(sel);
    [cents(oid).MedianX(2),cents(oid).MedianX(1)] = ind2sub(size(img),medind);
    indlist2 = find(img'==oid);
    sel2 = max(1,round(length(indlist2)*0.5));
    medind2 = indlist2(sel2);
    [cents(oid).MedianY(1),cents(oid).MedianY(2)] = ind2sub(size(img'),medind2);
end
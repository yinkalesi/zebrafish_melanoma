function  L_final  = return_sub_listL(L,listOfIndex)

L_final = zeros(size(L),class(L));

for i=1:length(listOfIndex)
    index = listOfIndex(i);    
    L_final(L==index) = index;
end

end


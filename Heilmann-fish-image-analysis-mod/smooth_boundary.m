function [x_s,y_s] = smooth_boundary(x,y,delta) 


x_s = x;
y_s = y;

max = length(x);
min = 1;


for i=1:length(x)
    
    count = 0;
    temp = 0;
    
    for k=-delta:delta 
        
        index = i+k;
        
        if (i + k)<min
            index = max + (i + k);    
        end    
        if (i + k)>max
            index = -(max - (i + k));    
        end    

        temp = temp + x(index);
   
        count = count + 1; 
    end
    
    x_s(i)= temp/count;
end

for i=1:length(y)
    
    count = 0;
    temp = 0;
    
    for k=-delta:delta 
        
        index = i+k;
        
        if (i + k)<min
            index = max + (i + k);    
        end    
        if (i + k)>max
            index = -(max - (i + k));    
        end    

        
        temp = temp + y(index);
   
        count = count + 1; 
    end
    
    y_s(i)= temp/count;
end

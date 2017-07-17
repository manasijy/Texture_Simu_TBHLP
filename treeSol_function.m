function [sol] = treeSol_function(stPt,b,c)

new(1)=stPt;
sol(1)=stPt;

while ~isempty(new)
    first = new(1);
    new(1) = [];
    temp = multipleSol_function(first,b,c);

    for i=1:1:numel(temp)
        repeat = 0;
       
            for j = 1:1:numel(sol)
            if isequaln(sort(sol(j).B),sort(temp(i).B))
                repeat=1;
                break;
            end
            end
        
        if ~repeat                 
           new = [new,temp(i)];
           sol = [sol,temp(i)];
        end
    end
end
end
    

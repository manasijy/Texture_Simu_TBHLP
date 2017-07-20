function Arc = get_Arc(DCMtrx,eR,CRelaxed)
%GETARC Summary of this function goes here
%   Detailed explanation goes here

    RC13 = DCMtrx'*eR(:,:,1)*DCMtrx;
    RC23 = DCMtrx'*eR(:,:,2)*DCMtrx;
    RC12=  DCMtrx'*eR(:,:,3)*DCMtrx;
    A13 = [RC13(1,1),-RC13(1,1);RC13(2,2),-RC13(2,2);RC13(2,3),...
        -RC13(2,3);RC13(1,3),-RC13(1,3);RC13(1,2),-RC13(1,2)];
    A23 = [RC23(1,1),-RC23(1,1);RC23(2,2),-RC23(2,2);RC23(2,3),...
        -RC23(2,3);RC23(1,3),-RC23(1,3);RC23(1,2),-RC23(1,2)];
    A12 = [RC12(1,1),-RC12(1,1);RC12(2,2),-RC12(2,2);RC12(2,3),...
        -RC12(2,3);RC12(1,3),-RC12(1,3);RC12(1,2),-RC12(1,2)];
    
    switch CRelaxed
        case 1 %only 13
           Arc = A13;
        case 2  %only 23
           Arc = A23;
        case 3 %only 12
           Arc = A12;
        case 4 %13 and 23
           Arc = [A13,A23]; 
        case 5 % 13,23 and 12
           Arc = [A13, A23, A12];
        otherwise
           Arc = [];
    end
end


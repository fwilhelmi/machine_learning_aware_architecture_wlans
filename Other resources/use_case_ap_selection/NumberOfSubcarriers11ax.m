% -----------------------------------------
% -----------------------------------------
% Number of Subcarriers
% -----------------------------------------
% -----------------------------------------

function Ysb=NumberOfSubcarriers11ax(B)

switch B
    case 2.5 
        Ysb = 24;
    case 5 
        Ysb = 48;
    case 10 
         Ysb = 102;    
    case 20 
         Ysb = 234;
    case 40
         Ysb = 468;    
    case 80
         Ysb = 980;
    case 160
         Ysb = 1960;

    otherwise
        Ysb = 234;
end

end

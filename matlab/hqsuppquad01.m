function [xs,ws]=hqsuppquad01(inode)
%
%       This subroutine returns to the user special-purpose
%       quadratures. The quadrature number k (k=1) integrates
%       functions of the form
%       
%       f(x)=phi(x)+psi(x)*log(|x-x_k|)+1/(x-x_k)+1/(x-x_k)^2,   (1)
%       
%       with phi and psi polynomials of order 2, and x_k one of the
%       Gaussian nodes on the interval [-1,1] (in a 1-point quadrature)
%
%       approximate accuracy of the quadratures:
%
%       smooth, smooth x log - 15 digits ( 34 in extended precision )
%       hilbert - 14 digits ( 33 in extended precision calculations )
%       hypersingular - 11 digits ( 30 in extended precision calculations )
%
%       Input parameters: 
%
%       inode - the node number of 1-point Legendre quadrature
%
%       Output parameters:
%
%       xs - the nodes of the quadrature
%       ws - the weights of the quadrature
%
%       In this routine, all auxiliary quadratures have 4-points
%       
%

xs1 = [
       -3.9021169606396825495153935800135850E-01,
        7.9384539051302290615351549765250394E-02,
        2.3175104705465428376862752542177513E-01,
        7.5144000798208693575457987189466582E-01];
ws1 = [
        9.5781216530792264883208189896560360E-01,
       -1.7323690400943302832555528369457516E-01,
        1.0117916168666548753268745403809136E+00,
        2.0363312183485550416659884434805797E-01];

if ( inode ==  1 ) xs = xs1(:); ws = ws1(:); end

end

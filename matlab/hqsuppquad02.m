function [xs,ws]=hqsuppquad02(inode)
%
%       This subroutine returns to the user special-purpose
%       quadratures. The quadrature number k (k=1,2) integrates
%       functions of the form
%       
%       f(x)=phi(x)+psi(x)*log(|x-x_k|)+1/(x-x_k)+1/(x-x_k)^2,   (1)
%       
%       with phi and psi polynomials of order 4, and x_k one of the
%       Gaussian nodes on the interval [0,1] (in a 2-point quadrature)
%
%       approximate accuracy of the quadratures:
%
%       smooth, smooth x log - 15 digits ( 34 in extended precision )
%       hilbert - 14 digits ( 33 in extended precision calculations )
%       hypersingular - 11 digits ( 30 in extended precision calculations )
%
%       Input parameters: 
%
%       inode - the node number of 2-point Legendre quadrature
%
%       Output parameters:
%
%       xs - the nodes of the quadrature
%       ws - the weights of the quadrature
%
%       In this routine, all auxiliary quadratures have 6-points
%       
%

xs1 = [
       -8.9018677908910660463189003294432592E-01,
       -6.2870693768375495808309609585526878E-01,
       -5.8464825915187138698880412714446677E-01,
       -4.7799194632563004934063903137707996E-01,
        9.3742857272706347732652393454913629E-02,
        7.8358672098522799016949945645681307E-01];
ws1 = [
        2.5055978107570433593066567783447397E-01,
        1.7212370388827705426846696190976309E-01,
       -5.5985185826108120360462096645293828E-03,
        3.1958648966723507430961642561912064E-01,
        7.3810977053556049965269099996670659E-01,
        5.2521877341583384787460614433446509E-01];

if ( inode ==  1 ) xs = xs1(:); ws = ws1(:); end
if ( inode ==  2 ) xs = -xs1(:); ws = ws1(:); end

end

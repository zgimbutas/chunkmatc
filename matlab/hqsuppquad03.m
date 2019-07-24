function [xs,ws]=hqsuppquad03(inode)
%
%       This subroutine returns to the user special-purpose
%       quadratures. The quadrature number k (k=1,2,3) integrates
%       functions of the form
%       
%       f(x)=phi(x)+psi(x)*log(|x-x_k|)+1/(x-x_k)+1/(x-x_k)^2,   (1)
%       
%       with phi and psi polynomials of order 6, and x_k one of the
%       Gaussian nodes on the interval [-1,1] (in a 3-point quadrature)
%
%       approximate accuracy of the quadratures:
%
%       smooth, smooth x log - 15 digits ( 34 in extended precision )
%       hilbert - 14 digits ( 33 in extended precision calculations )
%       hypersingular - 11 digits ( 30 in extended precision calculations )
%
%       Input parameters: 
%
%       inode - the node number of 3-point Legendre quadrature
%
%       Output parameters:
%
%       xs - the nodes of the quadrature
%       ws - the weights of the quadrature
%
%       In this routine, all auxiliary quadratures have 8-points
%       
%

xs1 = [
       -9.8834081101812372755101208145644071E-01,
       -9.1106342251358617124237981718155729E-01,
       -7.9017250398601508350230001724889062E-01,
       -7.8362207848298192609238369251251881E-01,
       -7.1053945493412104876216442778981374E-01,
       -3.0502282163895649524271654878598379E-01,
        3.2543499255879755015543810790203096E-01,
        8.5379729942367812704652547826708594E-01];
ws1 = [
        3.2100478936684566923596061585155329E-02,
        1.2107629210152451611464065988244450E-01,
        1.0638118893223525186237171513583022E-01,
       -4.1242952916840648188174373258451049E-02,
        2.1497585699874996561203619476592146E-01,
        5.6375135982132651916092941876263901E-01,
        6.3941310947053241847062680784138584E-01,
        3.6354466665578741004397351528507469E-01];

xs2 = [
       -8.7362688438437260785266579227391650E-01,
       -5.0261003671421670134197704997994573E-01,
       -3.6003083507176399168981937234820361E-01,
       -4.4744131307220227118113152966213738E-02,
       -1.5287691374923584286750777812130971E-02,
        7.3464236936785937707765770508643170E-02,
        4.4550540748659417818485445046136613E-01,
        8.7094653777748545666357317042825118E-01];
ws2 = [
        3.0721561788249431435563416956359458E-01,
        2.6671158722955189395056123352045251E-01,
        2.3955280134630714032691697065633129E-01,
        2.3257178004615352020254345681145783E-01,
       -3.8080972558452362012909129819067962E-02,
        2.0882328024381968167785730895125218E-01,
        4.6764948155804011806137414818046929E-01,
        3.1555642425208569343802184213551029E-01];

if ( inode ==  1 ) xs = xs1(:); ws = ws1(:); end

if ( inode ==  2 ) xs = xs2(:); ws = ws2(:); end

if ( inode ==  3 ) xs = -xs1(:); ws = ws1(:); end

end

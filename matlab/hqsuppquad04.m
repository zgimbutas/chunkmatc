function [xs,ws]=hqsuppquad04(inode)
%
%       This subroutine returns to the user special-purpose
%       quadratures. The quadrature number k (k=1,2,3,4) integrates
%       functions of the form
%       
%       f(x)=phi(x)+psi(x)*log(|x-x_k|)+1/(x-x_k)+1/(x-x_k)^2,   (1)
%       
%       with phi and psi polynomials of order 8, and x_k one of the
%       Gaussian nodes on the interval [-1,1] (in a 4-point quadrature)
%
%       approximate accuracy of the quadratures:
%
%       smooth, smooth x log - 15 digits ( 34 in extended precision )
%       hilbert - 14 digits ( 33 in extended precision calculations )
%       hypersingular - 11 digits ( 30 in extended precision calculations )
%
%       Input parameters: 
%
%       inode - the node number of 4-point Legendre quadrature
%
%       Output parameters:
%
%       xs - the nodes of the quadrature
%       ws - the weights of the quadrature
%
%       In this routine, all auxiliary quadratures have 10-points
%       
%

xs1 = [
       -9.7564792804921816851984386237126283E-01,
       -9.0839251703224047422604030528548469E-01,
       -8.6451345053869770280855070910529140E-01,
       -8.6345212242682579000558749809379786E-01,
       -8.3104967122996860916475041607198450E-01,
       -6.1252434621057351973844444209405678E-01,
       -2.1547703798089989650945783221927837E-01,
        2.6287761695885293507894792170078555E-01,
        6.8641318267066385754410993914952496E-01,
        9.4123907229389612725424383668475750E-01];
ws1 = [
        5.6698344629208295180266710646266243E-02,
        6.1845614697125125642708683680737619E-02,
        2.6227297722917397649779841556945409E-02,
       -1.3245449184929613197048151592620155E-02,
        1.1187238804792708852687584975705934E-01,
        3.1897341981359658318654371884577039E-01,
        4.5861788497536319939749836455774611E-01,
        4.7429812873256818193422431082465482E-01,
        3.5251711954289912968265527601524346E-01,
        1.5219525102332461199649539570819675E-01];

xs2 = [
       -9.4769626885843111764068779866700967E-01,
       -7.5617984109662602766760459978130501E-01,
       -5.2048330329671503846309074920816108E-01,
       -3.6188387349135839629331141203429169E-01,
       -3.4411518599277971011212276939456573E-01,
       -3.0019564029787006807571559118784166E-01,
       -6.1098583114108985890447068444782075E-02,
        3.2131455057535777046934980396420312E-01,
        6.9855867605490552756322955951099448E-01,
        9.4130805832488728319427943987734668E-01];
ws2 = [
        1.3037258625532394428292099442699821E-01,
        2.3458550376030678424736428376103906E-01,
        2.1500895419216868550595388256777536E-01,
        8.3932117608365382966050989039063640E-02,
       -4.6207233655864082710635802508014734E-03,
        1.2540848722027976579761228375229433E-01,
        3.3327379887776417268727931004833051E-01,
        4.0550680182285252757505529988109238E-01,
        3.2632348030237129286307539729288112E-01,
        1.5020899332615385234575113948132687E-01];

if ( inode ==  1 ) xs = xs1(:); ws = ws1(:); end
if ( inode ==  2 ) xs = xs2(:); ws = ws2(:); end

if ( inode ==  4 ) xs = -xs1(:); ws = ws1(:); end
if ( inode ==  3 ) xs = -xs2(:); ws = ws2(:); end

end

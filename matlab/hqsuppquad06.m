function [xs,ws]=hqsuppquad06(inode)
%
%       This subroutine returns to the user special-purpose
%       quadratures. The quadrature number k (k=1,2,...,6) integrates
%       functions of the form
%       
%       f(x)=phi(x)+psi(x)*log(|x-x_k|)+1/(x-x_k)+1/(x-x_k)^2,   (1)
%       
%       with phi and psi polynomials of order 12, and x_k one of the
%       Gaussian nodes on the interval [-1,1] (in a 6-point quadrature)
%
%       approximate accuracy of the quadratures:
%
%       smooth, smooth x log - 15 digits ( 34 in extended precision )
%       hilbert - 14 digits ( 33 in extended precision calculations )
%       hypersingular - 11 digits ( 30 in extended precision calculations )
%
%       Input parameters: 
%
%       inode - the node number of 6-point Legendre quadrature
%
%       Output parameters:
%
%       xs - the nodes of the quadrature
%       ws - the weights of the quadrature
%
%       In this routine, all auxiliary quadratures have 14-points
%       
%

xs1 = [
       -9.9263539674109127951208289875338411E-01,
       -9.6900133677332799985661815430730046E-01,
       -9.4577862695548045942422683687033824E-01,
       -9.3340655954049423889661983497951881E-01,
       -9.3326490226847310557212072591995062E-01,
       -9.2186769004835663003019862855539729E-01,
       -8.4111091008303046748816939182953589E-01,
       -6.7612293081805142244273288324096140E-01,
       -4.3032974761425662754273111618809442E-01,
       -1.2182084195222929481744625618195092E-01,
        2.1596490998337817144198586887026984E-01,
        5.3898268644243360871966740338903819E-01,
        8.0060515822462550078434071488098086E-01,
        9.6084710558935951044534373371995997E-01];
ws1 = [
        1.7759892583111480076865354108885711E-02,
        2.6103557802852222591270744823007779E-02,
        1.7459199447556483117143707023529640E-02,
        1.4109350897967189689648641640424017E-02,
       -1.0489919181418865856736954378426588E-02,
        3.9782371851763641387987849624108210E-02,
        1.2249421283740452104426069073019616E-01,
        2.0704761799926175785717096695879296E-01,
        2.8152906628649766656077681875262145E-01,
        3.2983170017261991283943797021116446E-01,
        3.3826026413396604548674772311118098E-01,
        2.9984553341296966566328386489743342E-01,
        2.1658947179814221767664950644809970E-01,
        9.9677679957306061865493116048982101E-02];

xs2 = [
       -9.8223096863339573716992651309042026E-01,
       -9.1351471374819150702682109684819832E-01,
       -8.1594910634778340291380982511344579E-01,
       -7.2368580166848394760881043052569361E-01,
       -6.6759108369554886262596169218187621E-01,
       -6.6341859726846234658162132414300899E-01,
       -6.4733144289991528139555114876133872E-01,
       -5.5635831514383147327164210186065958E-01,
       -3.7138013210542899129734004021205847E-01,
       -1.0405942581019993399495392947727451E-01,
        2.1271887706544872236337838829775301E-01,
        5.3029261352075848156380039343924971E-01,
        7.9488488841436519241643338498489054E-01,
        9.5951495152699091534388236786964563E-01];
ws2 = [
        4.4755648469260165388911644325300559E-02,
        8.8433922664748870715659200701819467E-02,
        1.0079002581368123410788858512913621E-01,
        7.8222599064518523270537908038239413E-02,
        3.1785310526608897740828010725238477E-02,
       -5.1437282859043003601897193127402695E-03,
        4.3674158624101622910663722716822483E-02,
        1.3820494449105709192826144861747348E-01,
        2.2990054921218380031001727049508450E-01,
        2.9915092906350724920635800855744597E-01,
        3.2610945589318822108651864621498570E-01,
        2.9987089583413638948016666079735131E-01,
        2.2130059004685918239107532243129188E-01,
        1.0294469858205305182330329056255082E-01];

xs3 = [
       -9.7213550472023728629029947919002685E-01,
       -8.6083382378829845916131553709765282E-01,
       -6.8977770562380483546972878491458713E-01,
       -5.0135377873065468409740735520140462E-01,
       -3.4215803649309629277104263522522194E-01,
       -2.5125850592619046006780451151893416E-01,
       -2.4011709817935887446434487820149153E-01,
       -2.2053662285757197566203528715108472E-01,
       -1.0822518673766023407236120209411290E-01,
        9.5206922526254982952429396279273504E-02,
        3.5136919290750034426963647458824428E-01,
        6.1260085070648728740456280044323097E-01,
        8.3092142117285983957320628617000864E-01,
        9.6665132020297154033813903963298121E-01];
ws3 = [
        7.0616961105411000930185129502804386E-02,
        1.4731601672748015247723375278516353E-01,
        1.8754490223138541044217440041339665E-01,
        1.8133946070268383503469109832476019E-01,
        1.3025883915370569655785920883484744E-01,
        4.5654387832923007138268411777560245E-02,
       -1.0984004579504122868799817138059798E-03,
        5.6364936344650078736090057506934296E-02,
        1.6305351358315586842637841573880544E-01,
        2.3725818684212538159107976459344908E-01,
        2.6697000745781424684328661859838070E-01,
        2.4735029716486736738540865036066602E-01,
        1.8255562994608317152250218273629091E-01,
        8.4815261365665195201722290540747089E-02];

if ( inode ==  1 ) xs = xs1(:); ws = ws1(:); end
if ( inode ==  2 ) xs = xs2(:); ws = ws2(:); end
if ( inode ==  3 ) xs = xs3(:); ws = ws3(:); end

if ( inode ==  6 ) xs = -xs1(:); ws = ws1(:); end
if ( inode ==  5 ) xs = -xs2(:); ws = ws2(:); end
if ( inode ==  4 ) xs = -xs3(:); ws = ws3(:); end

end
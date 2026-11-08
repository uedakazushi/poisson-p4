# 最後HQQ core全T19曲線：独立査読

2026-11-08、integrator `/root`。新Ga²曲線を独立手検証した。旧D102/M47の接23/軌道22診断は不変であり、局所成分性を推測しない。新point/operator/jetは無し。同一e/s曲線の別作者exactreadbackは別証拠。

真T19全completionのcoupledPac=2tu*bH/(r-v)を保持し、t=(r-v)/2,u=r-vとする。H-Casimir全六はproducer(3)、rawとの全Euler差はL=(0,-rb,(-3r+v)c/2,(-r-v)d/2,0)。q3=cd-b²,q4=ab+dH、actualD=-r(r-v)Hq3q4。

X=d∂b+2b∂c-a∂H,Y=-d∂a+b∂Hは[X,Y]=0、X/Y(q3)=X/Y(q4)=0。合成実linearflow a'=a-mu*d,b'=b+lambda*d,c'=c+2lambda*b+lambda²d,d'=d,H'=H-lambda*a+mu*b+lambda*mu*d はdet1、逆lambda/mu符号反転。lambda=mu=1/eでsourceHのinversepullbackはH+(a-b)/e+d/e²。forwardH'の逆符号を使わない。

dchart計算補助(U,V,I,J)=(b,-a,c-b²,H+ab)、逆(a,b,c,H)=(-V,U,I+U²,J+UV)。非線形補助をPGL同値には用いない。chart括弧はPij−xj*Pid+xi*Pjd。
例えば{a,b}_chart=-vab+(r-v)H/2だから{U,V}=(r+v)UV/2+(r-v)J/2。
{b,c}_chart=-(r-v)b(c-b²)だから{U,I}=-(r-v)UI。
-{a,c}_chart+2b{a,b}_chart=(r+v)a(c-b²)だから{V,I}=-(r+v)VI。
同様に全残三式は{U,J}=vUJ,{V,J}=-vVJ,{I,J}=2rIJ。Ga²はU,Vtranslation、I,J固定。
r=e,v=-e+e²sigma、sigma=4j*s,j²=-3を代入しU,V→U-1/e,V-1/eとすれば、全六の負e冪が消えproducer(10)の全e次数≤2式となる。sは一切反転しない。dchart restrictionはglobal quadratic projective spaceから固定injectivelinearなので、そのimageにgeneric所属する全polynomial curveは全eでimageに入り、canonicaldivzero GL-equivariantsectionが全150polynomial liftを与える。genericd-Casimirを仮定しない。全Jacobiはgenericparentからpolynomialidentityとして延長される。

中心のh=d-Casimir全六(a,b,c,H)は
Pab=sigma*h²/2,Pac=sigma*b*h,PaH=hH+ab-sigma*a*h/2,
Pbc=2(hc-b²),PbH=hH+ab+sigma*b*h/2,PcH=2ac+2bH+sigma*b²。
固定actualmapはx=(c-3H)/36,y=-(a+b)/4,z=-j(c+3H)/18,w=j(b-3a)/6,h=d、source順(a,b,c,d,H)でdet=-1/216、逆a=-y+jw/2,b=-3y-jw/2,c=18x+3jz,H=-6x+jz,d=h。

全六linear pushを手で展開した。R=hx-y²/2,S=hz-yw,U0=R-w²/24,W0=S-jw²/6に対し
Pxy=U0/2-shw/12,
Pxz=2xw/3+jzw/18-yz/3+s(y²-w²/12+jyw/3),
Pxw=W0/6+shy,Pyz=-W0/2+shy+js*hw/3,Pyw=sh²,Pzw=-2U0+shw。
詳細な独立符号チェックは次の通り。
Pxy=h(c-3H)/72-b(b+3a)/72-sigma*h(b-3a)/288。U0/2のbase部とsigma=4jsから一致。
Pxz=-j*(2ac+2bH+sigma*b²)/108、ac+bH=-6jyz+12jxw-3zw、b²=9y²+3jyw-3w²/4から全s項も一致。
Pyw=-j*sigma*h²/12=sh²。
Pxw=j*(-2h(c+3H)+2b(b-3a)+(9sigma*h/2)(a+b))/216、base=W0/6、mag=shy。
Pyz=j*(2h(c+3H)-2b(b-3a)+sigma*h*(5b-3a)/2)/72、base=-W0/2、mag=shy+js*hw/3。
Pzw=(-2h(c-3H)+2b(b+3a)+(3sigma*h/2)(b-3a))/36=-2U0+shw。
以上は元eta=1/12,xi=-j/36,kappa0=-1/6 sourceそのもの。canonicallift中心との差E∧N/5を明記し、literal同一を強要しない。

actualDはGa²inversepullbackとtotaldetから(2-e*sigma)*(d+e(a-b)+e²H)*q3*q4/216。
q3=18(R+jS/6+w²/24),q4=-6(R-jS/6-w²/8)なので中心D=-h*(R+jS/6+w²/24)*(R-jS/6-w²/8)=h*(-U0²-W0²/12)≠0。親period開-e²(2-e*sigma)≠0は全sで非空、有限除外根は明示polynomialliftの延長で扱う。

全sの中心は閉T19へ入る。xi他符号と元full20の全frame/磁気/lowerは既五族normalizer逆を戻す。最後coreの全支持T19帰属を受理、新成分増分零。HQQ全factor型の合成は旧SS/plane/empty全支持依存を別監査してから行う。全E/上界/主分類は未証明。

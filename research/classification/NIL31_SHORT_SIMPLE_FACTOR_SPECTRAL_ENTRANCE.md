# nil31短鎖 r=2/3：simple hyperplane を infinity に取る全源入口

2026-11-08、PRIMARY E。root が新因子flagを提案し、`/root/quartic_review` が全式を手で再導出した。
old infinity h を保つ marked-multipleH SS2 closure の否定はこの入口を扱わない。
actual simple factor Z を新 infinity に取るため、同じ fixed-old-H spectral profileの反復ではない。
本稿の時点で actualcurve の存在・全源E所属は未証明。

## 1. 全 projective tensor を保つ Casimir gauge と regular inverse

元全canonical源(15)の座標は(h,X,Y,Z,W)、t!=0,d任意。
actualD は非零定数倍 h²Z(hY−ZW)、Z は simple factor。

    L=(0,−th+2Y/3,2Z/3,0,2h/3), P'=P−E∧L

とすると P' の Z-row は正確に零。
divP'=N−5L=(5th−7Y/3)∂X−7Z∂Y/3−7h∂W/3。
scale−3/7 と actual invertible linear change

    H=Z,x=X+15tW/7,y=Y,z=W,w=h

を施す。inverse は h=w,X=x−15tz/7,Y=y,Z=H,W=z。
新divは N0=y∂x+H∂y+w∂z。
新 H-Casimir全六括弧は

    Pxy=−3dHw/7+27twy/49−27tHz/49−y²/7+2Hx/7,
    Pxz=−3dw²/7−3yz/7+2wx/7−9twz/49,
    Pxw=−9tw²/49−2wy/7,
    Pyz=wy/7−3Hz/7, Pyw=−2Hw/7, Pzw=−2w²/7。

Euler correction は同じ projective tensorを表す全 gauge変更、線形変換は actual PGLである。
全t,dを残し、後の有限jobは一つの稠密軌道 t=d1 のみを使う。
元六源→canonical inverse、td!=0一軌道、d0の元既約多項式源内閉性は旧支持ノートの証明を使う。
新actualDは H·w²(wy−Hz) の非零倍。
quartic は reducible なので、nil22longの既約quartic全源定理を直接適用しない。

## 2. 新 infinity に整合する spectral frame

    Nε=y∂x+(H+4ε²x)∂y+w∂z+ε²z∂w

は N0へ正則に戻る。generic linear eigenvalues は(-2,-1,1,2)ε、constant partはimage内。
zero hyperplane H は長い Jordan chain に残る。
全 spectralframeとinverseは

    X=H+4ε²x−2εy,Y=w−εz,Z=w+εz,W=H+4ε²x+2εy,
    x=(X+W−2H)/(8ε²),y=(W−X)/(4ε),
    z=(Z−Y)/(2ε),w=(Y+Z)/2。

Hは固定、determinant−32ε⁴。
frame=M·diag(1,ε²,ε,ε,1) なので、全 H-Casimir quadratic tensors の relativeweightは[-3,3]。
ε times inverse pullbackは pole高々二、regular oldsource を ε^-1 forwardした九係数のpoleは高々四。
parent tail ε³以上は負冪とendpointへ寄与しない。

## 3. 一回の有限 sufficient construction

公開 SS12九parameter全源と七真正Poisson relationsを保持し、
n_i の Laurent[-4,2]全63係数を許す。
最初に全negativepowers/150endpointを exact QQ線形消去する。
consistentなら元七関係を全epsilon次数で恒等的に課し、
linear又は明示linear-factor supportで閉じる場合だけ、boundedに一つのrational solutionを探す。
unsupported nonlinear residualでは保存停止し、GB/profile/jet列を追加しない。

この有限Laurent ansatzが nonlinear formal tails全てを尽くすとは仮定しない。
literal solutionが得られた場合だけ、別の sparse Fraction Laurent演算で
元六括弧/全350/全150endpoint/div/Pfを全展開し、generic parentの被覆hypothesisを個別に確認する。
source九係数をfree Poissonと呼ばず、真のproductparametersを保持する。
generic SS12 fieldだけで既知閉和所属を宣言しない。

## 4. 明示した未解決部分

新 infinity に移すこと自体はE所属を与えない。
成功したactualcurve＋generic既知parent認識が揃った場合だけ、
tdunit orbitと元六parameter源の閉性へ戻して全源被覆を結論する。
negative/unsupportedの場合にはこの全データを保存し、wholeEをOPENに残す。

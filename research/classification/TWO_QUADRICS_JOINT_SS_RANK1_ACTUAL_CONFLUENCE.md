# H＋相異二既約quadric・joint半単純rank一：全源のT16/T17実衝突

2026-11-08、PRIMARY E、rootの新しい純証明。独立読解前。
新CAS、operator、kernel、minor、profile、point又はjetは無い。

**定理。** actual退化因子がH＋相異なる二つの既約quadricで、二normal場のjoint半単純rankが一、
旧lambda零入口にある全rank四Poisson sourceは既知MPS成分の有限和に入る。
非planeのB三型は具体的にT16又はT17へ、A二型は旧空性又は既全plane被覆へ戻る。
全higher、全nil方向、全periodを保つactual曲線を全parameter上で構成する。

## 1. 元全sourceと因子条件

`TWO_QUADRICS_JOINT_SS_RANK1_SUPPORT.md` の全五型還元を使う。
二normal場自体がspan(S,N)に入るとは仮定せず、全追加nilpartを保持した必要supersourceである。
元全200×92方程式とJacobi/Pfの全枝は独立Fraction review済み。
入力は `CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_SUPPORT_20261107.json` の全括弧表と
`CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_SUPPORT_REVIEW_20261107.json`。
全H-Casimir chart h一の座標は(x,y,z,w)、hを戻すと以下は全global quadratic tensorである。

A1は全Pfにy因子。A2のsemisimple係数非零も全Pfにyw因子。
A2のsemisimple係数零は全nil22-image sourceへ実際に戻り、
後続の `TWO_QUADRICS_PLANE_WHOLE_COVERAGE.md` が全点を既知成分へ入れる。
旧note時点のplane未完文を現在の残余と誤読しない。
残るものはB1/B2/B3の全source表である。

## 2. 使用する真のMPS全parent

H-Casimir coordinates(H,A,B,C,D)を使う。

**T16全source**（旧named親同定で全periodとの正逆binding済み）：

    {A,B}=pAB+uH²=F,  {C,D}=qCD+vH²=G,
    cross = k(A∂A-B∂B) wedge(C∂C-D∂D)。

p,q,u,v,kは全五parameter、全Jacobi零、actualPf=HFG。

**T17全source**（未完成edge和でなくcompleted source）：

    {A,B}=pAB+uH²=F,  {C,D}=qCD+vF=G,
    cross=(A∂A-B∂B) wedge(aC∂C+bD∂D),  b=p-a。

全五parameter p,q,u,v,aを保ち、GのvF項を独立係数へ緩めない。
Jacobiはa+b=pから恒等零、actualPf=HFG。
一般parameterの実像がnamed成分に入り、その非零像の閉包も同じ成分に入ることは
`TWO_QUADRICS_CONFLUENCE_NAMED_MPS_PARENTS.md` と独立
`TWO_QUADRICS_TAN_SIMPLE_SNC_SEPARATE_REVIEW.md` の既受理入力である。
以下ではepsilon非零のactual PGL像だけへこの定理を使い、中心は全括弧正則性と閉性で戻す。

## 3. B1の全六parameter：T17への実曲線

元の自由係数はn0,n3,n4,n5,n8,n9で、n0 n5 n8非零。

    F=n0 h²-2n5 xy,  α=-n8/(2n5),  β=n9,  ζ=z+(n3/n5)h。

元全六括弧は正確に

    Pxy=F,
    Pxζ=-x(n5ζ+n4w),  Pxw=-n5 xw,
    Pyζ= y(n5ζ+n4w),  Pyw= n5 yw,
    Pζw=αF+βw²。

αのbindingは元表のPzwから読み、n8Fと取り違えない。
β零なら全actualDはhF²の定数倍となり相異二quadricではないので、その補集合は対象外。
従って本対象ではβ非零である。

parent T17へ次の全periodとactual frameを入れる。

    H=h, A=x, B=y, C=w-epsilon ζ, D=w+epsilon ζ,
    p=-2n5, u=n0, q=-2epsilon β, v=-2epsilon α,
    a=-n5+epsilon n4, b=-n5-epsilon n4。

frameのdetは-2epsilon、逆はζ=(D-C)/(2epsilon),w=(C+D)/2及び元zへのtranslation。
空間点又は未知source係数を反転せず、epsilon非零と元n5unitだけを使う。
全sourceのG=qCD+vFを戻すと、targetへのinverse pushは正確に

    Pxy=F,
    Pxζ=-x(n5ζ+n4w), Pxw=-x(n5w+epsilon²n4ζ),
    Pyζ= y(n5ζ+n4w), Pyw= y(n5w+epsilon²n4ζ),
    Pζw=αF+β(w²-epsilon²ζ²)。

全六括弧がepsilon多項式で、epsilon零は元全sourceそのもの。
全actualPfもhF[αF+β(w²-epsilon²ζ²)]であり、中心が零になることは無い。
よってB1全対象点はclosed T17に入る。

## 4. B2の全二branch

n0零branchの元Pfはxyで割れるため、相異二既約quadricという本対象では空。
これは必要supersource上のparameter branchを因子条件で正当に戻したもの。

n0非零branchではn6 n9非零を保ち、実frame

    w'=w+(n5/n6)h, z'=z-(n4/n6)w,
    d=n3-n4 n5/n6, p=-n6,
    F=n0 h²-n6xy, α=-n9/n6, β=n8

を使う。元全括弧は

    Pxy=F, Pxz'=-d hx, Pyz'=d hy,
    Pxw'=p xw', Pyw'=-p yw', Pz'w'=αF+β h w'。

β零は再びhF²なので対象外。全n3,n4,n5は上のglobal affine PGLへ吸収し、逆は明示的に解ける。
parent T17を

    H=h,A=x,B=y,C=h+epsilon z',D=w',
    p=-n6,u=n0,q=epsilon β,v=epsilon α,
    a=-epsilon d,b=p+epsilon d

とする。detはepsilon、逆z'=(C-H)/epsilon,w'=D。
全六括弧は

    Pxy=F,
    Pxz'=-d x(h+epsilon z'), Pyz'=d y(h+epsilon z'),
    Pxw'=(p+epsilon d)xw', Pyw'=-(p+epsilon d)yw',
    Pz'w'=αF+β(h+epsilon z')w'。

全epsilon多項式で、中心は元全B2対象源、actualPfはhF times最後の二次式。
従って全B2もT17又は対象外factor branchである。

## 5. B3・n1非零branch：T17への三次pole完全相殺

元unitはn1 n6 n7非零、n3,n4任意。

    V=w²-2hz, Q=n4 h²+n6 V, Qε=Q-epsilon² n6 z²。

元の全括弧は

    Pxy=(n1/n6)Q-n7xy,
    Pxz=-(n3+2n6)wx, Pxw=-(n3+2n6)hx,
    Pyz=n3wy, Pyw=n3hy, Pzw=Q。

T17のframeと全periodを次に固定する。

    H=h,A=h+epsilon w+epsilon²z,B=h-epsilon w+epsilon²z,C=x,D=y,
    p=2epsilon n6, u=-2epsilon n6-2epsilon³n4,
    q=-n7, v=-n1/(2epsilon³ n6),
    a=epsilon(n3+2n6), b=-epsilon n3。

detは-2epsilon³（座標の並べ替えを含む符号は同じ）。
逆はw=(A-B)/(2epsilon),z=(A+B-2H)/(2epsilon²),x=C,y=D。
source Fは厳密に-2epsilon³Qε、Gは厳密に-n7xy+(n1/n6)Qεとなる。
特にvuH²を勝手に別parameterへ緩めてはいない。

全target括弧は上の元表のQをQεへ、Pxw/Pywのhを(h+epsilon²z)へ替えたもの。
他の四括弧もそのままであり、全てepsilon多項式となる。
actualPfはh Qε[(n1/n6)Qε-n7xy]。
三次のperiod poleとinverse frame poleは六括弧全てで正確に相殺している。
中心は全n1非零branchそのものであり、従って全T17。

## 6. B3・n1零branch：T16への実曲線

元unitはn0 n6 n7非零、n3,n4任意。元全表は

    Pxy=n0 h²-n7xy=F,
    Pxz=-n3wx, Pxw=-n3hx,
    Pyz=n3wy, Pyw=n3hy, Pzw=Q。

T16を

    H=h,A=x,B=y,C=h+epsilon w+epsilon²z,D=h-epsilon w+epsilon²z,
    p=-n7,u=n0,q=2epsilon n6,v=-2epsilon n6-2epsilon³n4,
    k=-epsilon n3

とする。detは-2epsilon³、逆は前節と同じC,D版。
G=-2epsilon³Qεで、Cartan field C∂C-D∂Dのinverseは

    (w/epsilon)∂z+(h/epsilon+epsilon z)∂w。

従って全六括弧は元表のPzwをQεへ、Pxw/Pywのhを(h+epsilon²z)へ替えるだけで、
全epsilon多項式、中心は元全branch、actualPfはhFQε。
よってこの全branchはT16に入る。

## 7. 全source・全補集合への結論

B各型の全normalizer inverseは旧全五型入口と上記explicit PGL inverseで戻す。
新たに反転したsource係数は旧unit n5又はn6だけ。相異性からβunitを得た零補集合も明記した。
全磁気/自由係数、境界n3=n4=0、B3のn1零/非零を含み、一般sourceのみの主張ではない。
全epsilon非零の点は真completed T16/T17 sourceのactual PGL像、中心のtensorは非零。
projective曲線は中心でも正則なので、closed named成分により全対象点が被覆される。

これは全像≤23という旧G上界より強いwhole E被覆である。
joint半単純rank零のpure-nil factor-specific源、他factor、重複H、全E又は総成分数は別義務。
新成分増分零。旧証明書と旧未完記述のbytesを改変せず、この新実曲線で閉じる。

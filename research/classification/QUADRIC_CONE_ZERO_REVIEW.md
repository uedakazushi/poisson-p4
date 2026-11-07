# Quadric cone全rank零枝のG非優越：別作者の全源・Ore・minor再構成

2026-11-07。主義務G。**独立再現は全項目PASS**。
quadric cone上で零となり、実五次Dが既約被約である全rank四支持のPGL5飽和射影像は23次元以下。
従って未知の25〜27次元G優越成分の一般元にはならない。
E帰属、rank四不存在、他の曲面又はG全体の完成は主張しない。

独立実装は `code/classification_quadric_cone_zero_review.py`。
producer二本、SymPy、native CASはimportも実行もしていない。
既済独立scrollレビューの標準ライブラリFraction・疎多項式・有理消去・行列式だけを再使用した。
Jacobi、全源、Ore恒等式、principal matricesは今回の別ループで組み直した。
generic complete kernelや新GBはこの上界に不要である。
契約は `jobs/quadric_cone_zero_review.json`、結果は
`CLASSIFICATION_QUADRIC_CONE_ZERO_REVIEW_20261107.json`。
事前gate後1.454秒、最大RSS18,124 KiB。

## 全六十源と全補集合

座標 `(u0,u1,u2,z,w)`、`q=u0u2−u1²`、`S=(w=q=0)` とする。
Koszul/Eulerから全sectionが `Π=wL+qA` に入るという幾何的入口を読み直した。
全L五十、全A十を実150係数へ独立に写すとrank60、Euler鉛直は24、合併rank79。
従ってprojective sectionは55次元で、重複五方向はL内のEuler∧constant vectorそのものである。
Aにはこの冗長性が無い。

prime CI squareのprimary性と、既済のgeneric abelian conormal条件を適用すると
`{w,q}∈(w,q)²` が全点の三次係数条件になる。
三十五三次単項式、六membership係数を残した35×66行列を全raw源から組み直した。
rank17、保存された49本のkernelの独立性・全積零を検査し、kernelの完全性を得た。
raw六十への全十一条件と、そのA像七次元も再構成した。
許されるAは正確に `A^{w,ui}=0` の三条件を満たす。
conormalの仮定は既済のnonplane multiplicity-two/Frobenius first-jet入口に依存し、
今回clean、ordinary node、vertical divisor無しを追加仮定していない。

`λ=A^{w,z}`、`b_i=A^{ui,z}`、base skewをCと書く。
H=(w=0)上の実五次は

    D|H=−λ q²(C12u0−C02u1+C01u2)

なので、λ=0又はC=0の全係数枝はH因子を持つ。
これらを既約Dの仮定だけで除外し、零baseを黙って捨てていない。
λ非零のnormal shear `ui'=ui−b_i w/λ` の逆と全constant Aの変換を、
λの分母を払った普遍多項式行列で独立検査した。
qの差がw*linearであるため、全Lに吸収できることも確認した。

dual linear formのq-normは `Δ=4C01C12−C02²`。
非零normのSO3軌道と、非零isotropic coneのSO3軌道は、直交basis又はWitt basisで全て覆える。
ここには有理Gram chartの未処理分母は無い。
非isotropic側では複素tensor scalarの平方根を取りΔ=−1へ正規化し、
vertex scaleでλ=1を戻す。
従って全Aはsemi `∂w∧∂z+∂0∧∂2` とnil `∂w∧∂z+∂0∧∂1` の二軌道。
元のallowed A strataの次元は7と6である。

## 全mod-H三十kernelと、二十五変数への可逆橋

両Aについて、元五十Lの全projective Jacobiを独立に展開し、wを含まない全係数を残した。
L²項はwで割れるという構造も保ち、conormal十一式と合わせたrankは20。
保存されたparticularと三十kernelの全積を検査し、kernel独立性とrankから完全性を確認した。
Euler五方向を加えても条件を変えず、非零五minorによるliteral gauge後のkernelは二十五。

Ore表示の25方向とEuler5方向を元の五十L係数へ独立に展開し、
保存された30×50 recoveryがleft inverseであることを全成分で検査した。
さらに旧三十kernelと旧particularとの差を全て回収できる。
従って新25表示は旧全必要sourceの一部ではなく、Eulerを除いた全sourceである。

chart w=1で、元のprojective bivectorからEuler補正したaffine representativeを直接計算した。

    P=B+(Y+zT)∧∂z,
    B=qC+u∧X+D−u∧a,  X=C*dq,
    T=a−d−cu−X,
    Y=−u q−u(b·u)−Mu−e.

25parameterを全て独立変数として、元のE∧Jacobiとこのaffine Jacobiの恒等式、
さらに

    Jac(B)=0,  Lie_T B=0,  Lie_Y B=Y∧T

への全係数同値を確認した。
特定parameter又は有限jetだけのチェックではない。

## 全B,T支持とc=0補集合

Jac(B)の空間次数一・二の係数は、D三係数とa三係数のrank三の有理線形条件。
その完全kernelは `a∈image C, D=δC` であり、代入後の残Jac(B)も恒等零。
このためBのPoisson全支持は三変数で尽くされる。

Tのtranslation係数にはconstant pivotsを用い、残る全多項式のQQ-rowspaceを比較した。
semiでは

    t=((c+1)p0,0,(1−c)p1),  c(δ−p0p1)=0.

従ってc=0ではp0,p1,δの三次元、c非零ではδ=p0p1の三次元で、
固定非零cでは二次元以下。

nilでは

    t1=p1−cp0/2, t2=0, t0自由,
    cp1=0,
    −4cδ+cp0²−2p0p1+2p1t0=0.

c非零ならp1=0、δ=p0²/4、自由変数p0,t0,cの三つ。
c=0では `p1(t0−p0)=0` の両閉枝とδ自由を残し、各々三次元以下となる。
δ、nilpotent base、p1=0、c=0を反転して捨てていない。
これらは全LTの係数理想の両方向をconstant pivotとQQ-rowspaceで検査した結果であり、
generic frame又は任意点から全支持を推測していない。

## 全Y fiberのrank下界と全例外

Yのradial quadratic三係数、linear九係数、constant三係数をこの順で解く。
空間次数三、二、一のdiagonal blocksを独立に再構成し、全matrix entryを保存値と一致させた。
B,Tの低次数係数及び前段の解はaffine右辺だけを変える。
従って各段のnullityの和は、全base値でY fiberの次元上界である。
fixed cubic項や未使用の最高・最低空間次数はfiberを小さく又は空にするだけである。

generic minorはQQ[c]上のfraction-free Bareissで厳密に再計算した。
既知kernelの再計算や有限素数でのgeneric rank推測はしていない。
全表示minorを、記載された有理根に対応するlinear factorで余り無しに割り切り、
残りが非零定数であることまで確認した。
従って複素数体上の未記載例外根は無い。

semiではgeneric rank下界 `(2,8,3)`、Y fiber≤2。
全例外c=−2,−1,1,2では保存された別minorも使って `(2,7,3)` を得、fiber≤3。
c=0も別に評価し `(2,8,3)`。
generic base3+fiber2、非零例外base2+fiber3、c0 base3+fiber2は全て5以下。

nilではgeneric rank下界 `(2,8,3)`、minorの全例外はc=0のみ。
c=0では `(2,7,3)`、fiber≤3で、normalized full source≤6。
generic側は3+2=5以下である。
これらの評価はrank下界の非零minorだけで十分であり、generic kernelの完全再計算は要求しない。

## 全PGL像へ戻す

Aut(S)はGO(q)、vertex変更、normal shear、w scalingから十二次元、tensor scalarを加えて十三次元。
normalized Poisson fiberの上界5、6に、全A軌道7、6を戻すと全固定S affine sourceはいずれも十二以下。
PGL5×非零source上ではAut(S)とscalarが十三次元の自由なfiber作用を持つ。
最初のPGL factorがAut(S)の元を復元するので、個別Poisson tensorのstabilizerが零である必要はない。
従って全射影飽和像は

    24+12−12−1=23

以下となる。
全constant Aの二軌道、全B,T分岐、全c例外、任意の十五Y係数と全Euler冗長性を保持した全支持上界である。
主義務Gのこのquadric-cone rank零枝を閉じるが、別の曲面・rank二枝・全G完了を含めない。

producer bytes、既存manifest、旧証明書は全て不変保存した。
checker初回は根0の単項式c8を「項数一なのでconstant」と誤判定してreport生成前に停止した。
次数で判定する実装へ直し、再gate後の最終実行で全項目PASSした。
失敗した実行を証明として用いていない。

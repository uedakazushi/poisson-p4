# smooth cubic零点を持つrank二planeの全normal weight源：既約quinticからの排除

2026-11-07、主義務G。S=P2⊂P4を実退化因子DのSingDに入るPoisson invariant planeとし、
π|Sはgeneric rank二、その零divisor Rはsmooth irreducible cubicと仮定する。
**結論。** rank四πの実Dは既約・被約quinticにならない。
横断ordinary double、normalization又は導手のsmoothness、vertical divisor不存在を仮定しない。
以下は全normal-weight源に対する幾何・代数証明であり、chosen tensorの接空間評価ではない。

## 1. conormal接続はcanonical scalar接続と完全に一致する

`PLANE_SMOOTH_CUBIC_CONORMAL_CONNECTION.md` の仮定・証明を使う。
SingD仮定により `{I,I}⊂I²` であり、cotangent algebroid接続は全Sで定義される。
smooth cubic R上の一位極とglobal meridian monodromyから、tracefree projective接続は零である。

base homogeneous座標をx=(x0,x1,x2)、normal座標をz,wとし、
R={f(x)=0}、fはsmooth ternary cubicとする。
canonical liftを

`P_f=f0 ∂1∧∂2+f1 ∂2∧∂0+f2 ∂0∧∂1`

とする。これはstrict Poissonで、同じplane bivectorを誘導し、conormal上にscalarなflat接続を与える。
元の接続とこのcanonical接続の差は、projective partが零なので `α*identity`。
αは全球 `H⁰(P²,Ω¹(R))` の元で、両接続のflatnessからdα=0である。
閉一形式の一位極はsmooth divisor上でlogarithmicである。
同noteの `H⁰(Ω¹(logR))=0` によりα=0となる。
従ってprojectivizationだけでなく、conormal接続全体がcanonical接続と等しい。

## 2. 全global tensorの有限normal filtration

normal dilation `(z,w)→(t z,t w)` によるweightを使う。
global bivectorは全二次homogeneous係数のliftを持ち、Euler∧linear fieldを商する。
Sのinvarianceと `{I,I}⊂I²` により、base-normal係数はI、normal-normal係数はI²に入る。
従って存在するweightは0,1,2だけであり、higher formal seriesは無い。

weight0のbase restrictionとconormal接続は§1で固定されるので、
その全projective sectionは

`P0=P_f+K(z,w) ∂z∧∂w`、K任意のbinary quadratic

となる。これはscalar mixed connectionを単に捨てたansatzではない。
base bivectorのlift差はEuler_base∧linear field、mixed normal一次係数はconormal接続そのものを表す。
両者を同時にglobal Euler gaugeで取り除いた後、残るweight0はnormal-normal二次係数だけである。

全weight1 sectionのraw係数は36、Euler gaugeは六次元である。
具体的にEuler∧`(zV+wW)` を使い、`z²∂i∧∂z` と `w²∂i∧∂w` を消す。
その結果、全三十次元sourceは一意に

`P1=zL+wM+z² D∧∂w+zw(B∧∂z+E∧∂w)+w² C∧∂z`。       (1)

L,Mはbaseのlinear係数bivector（各九係数）、B,C,D,Eはconstant base vector（各三係数）である。
ここでDはvector名であり、退化因子と同一視しない。
全weight2 sectionはbase constant bivectorをnormal quadraticで掛けた九次元source
`P2=z²A+zw A'+w²A''` だけである。

## 3. smooth cubicで使う唯一のinjectivity

constant base vectorXについて `U_X=[P_f,X]=-P_(Xf)` と書く。
smooth fは二変数だけのcubicではないので、`Xf=0⇒X=0`。
従って `U_X=0⇒X=0`。

linear base bivectorLを三次元volumeでlinear vector l(x)と同一視すると、
`[P_f,L]` は `grad(f) dot curl(l)` を係数とする三vectorである。
curl(l)はconstantなので、smoothnessから `[P_f,L]=0` はcurl(l)=0と同値であり、
Lはあるquadratic gのgradient bivectorP_gとなる。

以下の二正規形は複素数体上の全nonzero binary quadratic Kを尽くす。

## 4. K=a*z*w、a≠0：全weight1 kernelは座標変更六次元

normal二torusにより(1)は四characterへ分かれる。
`zL+zw E∧∂w` に対するstrict linear Jacobiは

`z[P_f,L]+zw(U_E−aL)∧∂w`。

Eulerとのwedgeの、二base方向とz,wを持つ成分は `z²w(U_E−aL)` なので、
projective Jacobiは `L=a^-1 U_E` を強制する。
これなら `[P_f,L]=0` は自動であり、これはまさにnormal→base field `z*a^-1 E` による座標変更である。
同様に `wM+zw B∧∂z` は `M=−a^-1 U_B`、field `−w*a^-1 B` の座標変更となる。

残る `w²C∧∂z` のEuler–Jacobiの三base方向とzの成分は
`w² Euler_base∧U_C`、すなわち非零定数を除いて `w² C(f)` である。
従ってC=0。`z²D∧∂w` も同様にD=0となる。
よって全closed weight1源は六つのnormal→base座標変更だけである。
genericなa又はfを取る議論ではなく、全a≠0と全smooth fに成立する。

## 5. K=a*z²、a≠0：重根も全く同じ六次元へ尽くされる

この場合、(1)のstrict linear Jacobiのnormal-normal成分は

`(-a z³ E+a z²w B+2a zw² C)∧∂z∧∂w`。

二base方向とz,wを持つEuler–Jacobi成分のz,w各係数は正確に

`z³: U_D−aL−a Euler_base∧E`,

`z²w: U_E−aM+a Euler_base∧B`,

`zw²: −U_B+2a Euler_base∧C`,

`w³: −U_C`。

最後から順にC=0、B=0、M=a^-1 U_E、
`L=a^-1 U_D−Euler_base∧E` を得る。
これらは全六parameter D,Eによるnormal→base座標変更そのものである。
実際 `zY` のorbit微分は `zU_Y+a z²Y∧∂w`。
`wZ` のorbit微分は `wU_Z−a z²Z∧∂z` で、後者の除去gaugeを施すと
`L=−a Euler_base∧Z, M=U_Z, E=aZ, D=0` となる。
従って全kernelを尽くし、他のJacobi成分も座標変更であることから自動的に消える。

## 6. nonzero Kではhigher weightも消える

§4–5によりP1は `X=zY+wZ` によるP0のorbit微分である。
Xはsquare-zeroなlinear matrixなので、`exp(X)=identity+X` は実際のunipotent PGL5変換である。
この変換でP1を消しても、global quadratic sourceの最高weightは二のままである。
従って変換後の全tensorは `P0+P2'`。

base constant bivectorAについて `[P_f,A]=0` なので、weight2 Jacobiは `[K,P2']=0` のprojective条件となる。
各normal quadratic Qに対し、Euler–Jacobiの二base方向とz,wの成分は

`Euler_normal∧K^sharp(dQ)∧A=−2KQ*A∧∂z∧∂w`。

nonzero Kは零因子でないため、P2'=0である。
よって全tensorはPGL5によりP_f+Kに共役し、その実quinticは

`D_actual=3*f(x)*K(z,w)`

（採用するPfaffian符号による全体scalarを除く）となる。
Kが相異二根なら可約、重根なら非被約であり、既約被約quinticにならない。

## 7. K=0もrank四点を持たない

(1)に対する二base方向とz,wのEuler–Jacobi成分は
`z³U_D+z²wU_E−zw²U_B−w³U_C`。
従ってB=C=D=E=0である。
残る式は `[P_f,L]=[P_f,M]=0` なので、§3により
`P1=zP_g+wP_h`、g,hは任意のquadraticとなる。
P2は初めからbase-base方向だけである。
よって全tensorP_f+P1+P2は三つのbase方向だけを持ち、projective rankは二以下である。
他のhigher-weight Jacobi方程式を解かなくてもrank四不存在が従う。

## 8. 結論と残す補集合

全K、全三十weight1、全九weight2を保持して、仮定を満たすrank二plane枝を既約被約quinticから排除した。
normal-normal quadraticの平方又は零をordinary-doubleと読み替えてはいない。
Rがsingular/reducible/nonreducedなplane、rank零plane、plane以外のSは補集合である。
新成分の個別E帰属や、既約非正規quintic全体の分類を主張しない。

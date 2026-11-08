# 不安定 elliptic 全源の実 T40 曲線：独立幾何レビュー

2026-11-08、PRIMARY E。**独立 PASS。**
任意の指定 rank四 source、任意の A₂⊕F₂,₃ と完全五 section frame を中心に、
actual ambient Poisson tensor の代数的曲線が構成され、中心は T40 に入る。
新CAS・数値・matrix・rank・point・jetの計算は無い。
以下は特に total image、relative flatness、conductor の全 fiber base change を独立に確認する。

## 1. 使う全入口と extension の全場合

既受理 genus一 recognition と独立reviewから、対象の全不安定正規化は
`Y0=P_B(A⊕F)`、B smooth elliptic、deg A=2、F stable rank二degree三。
L0=O(1)はample/GG、h0=5,h1=0、元の有限birational正規化写像は完全五 sectionである。
canonical／零 summand／nef非ample／cusp の補集合を旧入口で既に処理したことを確認した。
元 source の全 normalized tensor は conductor C0=-K_Y0 に沿って零となり、

```
H0(Λ²T_Y0(-C0))=H0(Ω¹_Y0)=p*H0(Ω¹_B)=C。
```

これは任意の閉点を覆う全source定理であり、一般bundleだけの主張ではない。

F⊗A^-1 は stable degree−1なのでH0零、RRよりH1一次元。
従って全ての指定A,FでExt¹(A,F)=C。
非零 cocycle ηの t倍で algebraic bundle family
`0→F→E_t→A→0` を B×A1 上に作れる。t0は指定split bundleそのもの。

t≠0のstabilityも全 subbundle で確認した。
destabilizing line M はdegree≥2。Fへ入ればstabilityに反し、Aへの非零mapなら
degree二でAと同型、extensionをsplitして矛盾する。
destabilizing saturated rank二 subbundleはdegree≥4、quotient line Qはdegree≤1。
Hom(F,Q)=0とHom(A,Q)=0によりE_t→Qは零となり矛盾する。
これでrank一／二の全destabilizing場合を尽くす。

## 2. 正規化写像の family と total Cartier image を先に作る

cohomology/base change と中心のH1零により、零を含むopen T上で
H0(E_t)はrank五 vector bundle。
L0のampleness/GGの開性を使い、全fiberでO(1)をample/GGにする。
そのsection bundleの中心の任意の指定frameは、Tを縮めてそのまま延長出来る。
frame変更による中心のPGL座標変更を勝手に許したのではない。

Y=P_(B×T)(E)はT上smooth/projectiveで、得られるf:Y→P4×Tはproper。
各fiberのmapはample L_tをpullbackするため曲線を縮約せずquasi-finite、従ってfinite。
同じ理由でtotal fもproper quasi-finite、finiteである。

各Y_tはintegral、完全五sectionの独立性によりimageはnondegenerate三次元。
degree formula `5=L_t³=deg(f_t)*deg(image)` と5の素数性から、
degree一image=hyperplaneは除かれ、全imageはdegree五、deg f_t=1。
この議論はgeneric fiberにも適用され、total fもbirational。
全f_tのsmooth source Y_tは、そのimageの正規化になる。

**次に**scheme-theoretic total image Dを取る。
Y integralからDはintegral horizontal prime divisor、ambient P4×TはregularなのでCartier。
TをA1のaffine openとしたためPic(T)=0であり、generic degree五から
I_D=O(-5)、一本のhomogeneous equation F∈O(T)[x0,…,x4]で定まる。
Fの全係数がt-t0で割れればtotal divisorにvertical componentが生じるので、
どのt0∈TでもF_t0は零多項式ではない。
それは既にdegree五と分かったirreducible image f_t0(Y_t0)上で零、かつ同じdegree五。
従って非零scalar倍のそのirreducible equationであり、**各scheme fiber D_t0は正確にそのreduced image**。

この順序が重要である。imageのbase changeやnormalization-flatnessを先に仮定していない。
Fは各ambient fiber上nonzerodivisorなのでD/Tはflatなrelative Cartier hypersurface。
total DはGorenstein/CM、各fiberはintegral Gorenstein hypersurface、
adjunctionからω_D/T=O_Dとなる。

## 3. MCM と conductor base change を全 local ring で確認

局所的に A=O_D、B=f*O_Y と書く。A⊂Bは同一fraction field内のfinite birational extension。
total Yはsmooth、従ってregular/CM。
finite morphismでの点のresidue field extensionは代数的、両total varietyのdimensionは同じなので、
Aの局所dimension dと、その上の各B局所環のdimension dは一致する。
Aのsystem of parametersは各B局所環でもparameter idealを作る。
各B局所環はCMなのでその全列はregular、従ってBはA上MCM。
単に「finiteだからflat/free」としていない。

AはGorensteinなので、MCM canonical-module dualityの正確な仮定の下で

```
Ext^i_A(B,A)=0 (i>0), Hom_A(B,A)=f*ω_Y/T
```

を得る。relative canonicalの同定はfinite dualityとω_D/T=O_Dによる。
total canonicalから書いても、同じbase canonical factorが両側で相殺される。

fraction field上のA-linear map B→Aはφ(1)によるmultiplicationである。
従ってevaluation at1は

```
Hom_A(B,A)=c={b∈B:bB⊂A}⊂A⊂B
```

という実conductor inclusionそのもの。
cはB-idealであり、dualizing line ω_Y/Tとの同定によりY上invertibleである。
このinclusionがanticanonical section sを与える。

τ=t-t0はAにもBにもnonzerodivisor。
0→A --τ→ A→A0→0へHom_A(B,-)を適用し、Ext¹零から

```
Hom_A(B,A)/τ = Hom_A(B,A0) = Hom_A0(B/τB,A0)
```

を得る。有限写像のfiberはaffine tensor productでB/τBそのものであり、
別のnormalization/base-change定理を仮定する必要が無い。
§2ですでにA0⊂B/τBをfinite birational正規化と証明したので、右辺は全fiberの正確なconductor。
この同型はevaluation at1と可換である。
従ってsの中心が余分なt倍で消えることは無く、全codimension二点を含むsheaf mapとして特殊化が正しい。
必要ならA0→B/τBのinjectivityとTor exact sequenceからB/AのT-flatnessも後で従う。

## 4. conductor係数の全descentと ambient Jacobi

非零closed α_B∈H0(Ω¹_B)を取り、PY=i_(p*α_B)s。
局所volumeで対応一形式はscalar倍のclosed α_Bなので、三次元Frobenius条件を満たしPYはrelative Poisson。
中心は上の一つの全source線を非零に張る。
定数scalarでα_Bを調整すれば、指定の元π_Y0と**同じtensor**に出来る。

df(PY)の局所係数は、conductor generatorにregular Jacobian係数を掛けたもの。
全係数がcB=c⊂Aに入るため、Λ²T_P4|Dの正則section ξへdescendする。
scalarだけ、又はnormalization上のrational tensorだけのdescentではない。
target tangent bundleのtransition matricesはA上にあるので、これら局所係数はglobalにglueする。

relative exact sequence

```
0→Λ²T_P4(-5)=Ω²_P4→Λ²T_P4→Λ²T_P4|D→0
```

とH0(Ω²)=H1(Ω²)=0からξは一意のregular ambient relative bivector Πへliftする。
これは各fiberだけのliftを後からglueする操作ではない。
中心のrestrictionは元π0と一致し、同じH0零性により中心Π0=π0。

relative smooth D上ではfは同型、ξはtangent Poisson。
ambient liftのtangencyはmod Fでnormal componentが零という条件なので、
Schouten bracketのrestrictionはintrinsic bracketと一致し、[Π,Π]|Dはdense smooth locusで零。
Dはintegral、restriction bundleはlocally freeなので全Dで零となる。
従って[Π,Π]はΛ³T_P4(-5)=Ω¹_P4のsectionであり、H0(Ω¹)=0から零。
これで全parameterの元ambient Jacobiを満たす実際の代数的曲線となる。

## 5. 実 Pfaffian、全 frame と T40 への閉性

smooth Dのtangent bundleはrank三なのでξ²=0、従ってactual PfaffianはFで割れる。
同じdegree五から `Pf_actual(Π)=c(t)F`。
global kernel O(5)(-D)=Oを使えばcはregular base functionであり、rational quotientのpoleを仮定していない。
中心Π0=π0は指定rank四なのでc(0)≠0。零を含むopenへ縮めればcはunit。
全fiberで実退化因子は構成したintegral reduced quinticそのものとなる。

t≠0の全E_tはstable rank三degree五。
旧独立genus一recognitionのstable bridgeは、Atiyahのcoprime determinant一意性と
degree五determinantへの全translationを使い、任意のこのpolarized normalizationを既知q5,2へ認識する。
完全五sectionのframe変更はactual PGL5、conductor-divisible normalized tensorの全spaceは一次元。
従ってdivisorだけでなくΠ_t自身のprojective classがT40に入る。
閉性により中心π0もT40に入る。

全A,F,Bおよび中心の全frameに対し同じ構成が出来る。
一つのmoduli universal bundleやgeneric中心への制限を仮定しない。
openへの縮小は各任意の中心を残し、tensorlineのscalarもその中心に正確に合わせている。
rank二/Pfaffian零、別degree、rank二restriction、他の全Eへの拡張は行わない。

## 6. 保存

独立確認は全て純証明であり、CASその他の新実行は不要だった。
全producer seven-file manifestの不変hashを照合し、本独立note/receipt/manifestを別に保存する。
旧G上界・normalization資料・producer proof・主稿・global・Gitを変更していない。

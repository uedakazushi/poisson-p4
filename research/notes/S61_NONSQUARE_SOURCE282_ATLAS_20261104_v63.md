# 37関係後の実source bundleと双対rank43 criterion

2026-11-04。基礎環は
`R=QQ[a,e,d,1/(a*d*F*B*C*J)]`、正則integralな次元3の環である。
元135式・同一Tn frameからの実319×239 mapを、sourceからtargetへ
`M:R319→R239` と書く。以下はその**同じcokernel**を扱う。

新たに行った算術は、sparse174×22 kernel matrixの構築と全22実恒等式の
一回の検査のみ。18.355秒、963非零entries、core numerator/denominator合計
14525 terms。minor/pivot search、quartic coefficient matrix、新点・field・GBは
計算していない。

## 1. 固定15 source rowsを実際に削除できる

既受理11 quartic pivotsと4 cubic pivotsから得た15関係は、319 sourceの15個の
座標への射影が `[[I11,*],[0,I4]]`、determinant1である。よって15関係を割った
sourceは、**その15 source rows以外の座標をbasisとするR304** である。

従って、保存された319 sourceから指定15行を単に削った固定304×239行列M304は、
元Mと同じimage/cokernelを持つ。新しいquartic係数を展開する必要はない。
15行は全てgenuine cubic/quarticなので、最初のquadratic174 sourceは残る。
`atlas.json` に15 deleted indices、304 retained indices、全source definitionsを
保存した。全行に同じD/Qhat規約を使い、T34/Tnを混在させない。

## 2. actual174×22 matrixと282 source atlas

実quadratic source174は `6 quadrics ×29 weight2 monomials`。
旧22関係をこの174座標へ埋め込んだ行列を `K:R22→R174` とする。
KのH-cofactor rowsは零である。

今回Kを元六式から直接構築した。15 Koszulは元六式basisを用いる。
旧六式basis `(q0,q1,q3,q5,ellZ,ellW)` との変換detは−80、Koszul15 basisの
変換はそのexterior squareでdet `(−80)^5`、既知unitである。残る6 linear-relation
倍と新しい22本目は保存されたL,ell,TZ,TWを元六式へ戻した。

lightのinverse frameは

`(x0,x1,x2,x3,x4,Z,W)=(v0,v1,v2,v4,v7,3v0−v1+v6,v3)`。

全22の `sum_i cofactor_i*q_i=0` をliteral polynomial arithmeticで検査した。
`kernel22.json` はKの全963 entriesを、実174 source番号とQQ(a,e,d)の
numerator/denominatorで保存する。既受理のpointwise independenceとwhole-kernel
定理によりKはFULLU上subbundle inclusionである。

`E152=coker K` はlocally free rank152。従って同じcokernelを持つ実mapは

`N:P282=E152⊕R130 → R239`。

これが正しい「282 source」である。**固定22minorを全域unitと仮定した単一の
282×239行列は、今回まだ構築していない。**

具体的atlasは次の通り。174個のquadratic source座標の22-subset sigmaについて
`delta_sigma=det K[sigma,:]` とする。D(delta_sigma)上では、その22 source rowsも
M304から削った実282×239行列が同じcokernelを与える。全delta_sigmaのopensは、
Kのpointwise rank22からFULLUを覆う。今回どのminorも展開・探索していない。

global free presentationの**存在**も示せる。E152はprojectiveなのでKには
left inverse Pが存在する。M304をmap向きで `[A B C]` と分け、Aは174 source、
Bは残る130 sourceのうち22、Cは残り108とする。その場合

`[A+B*P, C]:R282→R239`

は実Nのglobal free presentationとなる。実際、
`R174 ≅ E152⊕R22, g↦([g],Pg)` の逆は `([g],u)↦(1−KP)g+Ku`。
これは具体的Pを与えた後の明示式であり、Pそのものは今回未計算である。
Bézout completionやglobal basisを計算済みとは呼ばない。

## 3. Reflexive rank87 targetを使う別の実表示

target239には `pure light quartics210 ⊕ H*Sym²light28 ⊕ H²1` という固定分解が
ある。quadratic174のimageは、pure-quadratic image E146と、H*q_iの6方向に
分かれる。前者の商Q64は既受理のreflexive rank64、後者の6方向は六quadricsの
全FULLU定数独立性により28内のsubbundleである。その商をF22と書く。

従ってquadratic174をtargetから割った実加群は

`Q87 = Q64 ⊕ F22 ⊕ R`、reflexive rank87、pd≤1。

残る130個の実cubic/quartic sourceの像による

`R130 → Q87`

も元のfullquartic cokernelを正確に与える。Q87の非自由点は既知の有限
quadratic-defect集合に限られる。これは「Q87が全域free87」という主張ではない。
曲線一般点ではQ87をfree87として扱える。この表示はsource282の大きなbasisを
作らずに曲線欠損を扱う候補である。

## 4. Dual Fitting criterionとcodimension2

以下Cを実fullquartic cokernelとする。固定M319 codimension-transfer定理から
`codim Supp C≥2`。Nはgenericに全射なので、双対は単射であり

`0→(R239)*→P282*→D→0`。

Dはrank43、pd≤1。高さ1ではNが全射なのでDはfree43。高さh≥2では
depth lemmaから `depth D_p≥h−1≥1`。Dにtorsionがあればassociated primeが
高さ1以上に現れるが、上の二つの性質がそれを排除する。よってDはtorsion-free。
このlength1 projective resolutionを双対に戻すと

`C ≅ Ext¹_R(D,R)`、`Ext^i_R(D,R)=0 (i≥2)`。

したがって正確な代替gateは

* fullquartic全射 ⇔ D locally free rank43;
* 欠損曲線なし（Cのsupportが有限） ⇔ D reflexive。

後者は、正則2次元局所環ではpd≤1 moduleのS2がfreeと同値であり、
高さ3では既にdepth≥2だからである。**反射性だけでは有限閉点の欠損を消せない。**
既存Q64の反射性と、このDの反射性は別の主張である。

free282 chartでは

`Fitt_43(D)=I239(N)=Fitt_0(C)`, `Fitt_44(D)=I238(N)`。

この同じidealをbundle chartごとに計算できる。Fitting idealだけでなくcokernel
そのものが保存されているため、nilpotent/module structureを捨てる還元ではない。

## 5. kernel43のabstract free性では解けない

これは具体的反例で確認できる。dummy regular ring S=QQ[x,y,z]上で

`Ntoy=diag(I238,[x,y,0,...,0]) : S282→S239`

とする。最後は44 sourceの一行である。cokernelはS/(x,y)、欠損はcodimension2。
kernelは `S42⊕S*(-y,x)` という**free43**である。それでも欠損曲線上では
最後のkernel vectorがsource内で零になり、subbundleではない。
37個の零sourceを加えると319→239で37個のconstant独立kernel directionsと、
abstract free80 kernelを持つ同じ反例になる。

高さ2の実欠損一般点でもC_pは有限長、R_pは正則2次元だから、
Auslander–Buchsbaumによりkernel43のfree性は自動である。従ってその性質を
証明しても曲線は除去できない。必要なのはkernelの**source内のsplit独立性**、
または上のDの局所自由性である。

さらに `[x,y,z,0,...,0]` を最後の一行とする反例は欠損が原点だけで、Dは
reflexiveだがnonfree、C=S/(x,y,z)である。曲線除去後にもpoint gateが残る。
これらは論理的countermodelであり、新しいPoisson点の計算ではない。

## 6. Cyclic chartの44 scalar criterionとHilbert–Burchの限界

free282 chartで238×238の実minor Deltaがunitであるとき、Nを

`[[A238,B],[c,d]]`

へ並べれば、同じCは44個のscalar entries

`g=d−c*A^-1*B`、すなわち `C≅R_Delta/(g1,...,g44)`

で表示される。Deltaを掛けた `Delta*d−c*adj(A)*B` はbordered239 minors。
これは具体的な小さいlocal Fitting criterionだが、Deltaは追加chart条件である。
rank≤237の欠損曲線が無いとはまだ証明していない。単一generic cyclic chartの
計算だけを全欠損の証明へ広げてはいけない。

一般にfiber corank sなら、unit blockを消したlocal残余は
`s×(43+s)`。s=1の場合だけ上の44 scalar表示になる。height2 perfect scalar
idealならHilbert–Burch resolutionを使えるが、proper height2 idealも同じ
resolutionを持つ。従ってperfection/CM性/abstract free43を全射性の代用にしない。
有効な新certificateは、実際の44 entriesのunit ideal、またはDのS2/localfree性を
証明する追加のものとなる。

既存C122 chartでは22 quadratic relationsを既に割っている。受理済みSchur
操作後のmapは `(p,z)↦(p+A*z,C122*z)`、pは175 coordinatesである。そのkernelは
zへの射影でker C122と同型で、逆は `z↦(-A*z,z)`。従ってD15は122 source内でも
pointwise独立に降り、rank107 source bundle→rank64 targetを得る。
**同じSchur coordinate transportをまだ係数評価していないので、固定actual
107×64行列を構築済みとは呼ばない。** 完成すればcyclic chartでは63pivot後の
1×44 criterionとなり、残り43方向という数と一致する。

今回得たのは、固定304 source・sparseK22・282 atlas・双対rank43 gateを正確に
固定したことである。欠損曲線、残る点、full239/C122rank64は未解決である。

# 低rank marked modular旗の外：全rank四sourceの射影次元上界25

2026-11-08、PRIMARY `E_multiple_hyperplane`、担当 generic_exhaustion。
**純証明の全case在庫。新しい点・Jacobi行列・minor・jet・GBは計算していない。**
既存の全source定理に、新しい nil22-image の全parameter fiber 上界を加える。
Gの非優越又はprimitive β像の上界を、全Poisson像の上界へ読み替えない。

## 1. 閉旗と主張

射影Poisson空間のgeneric rank四開集合を \(\mathcal P_4\) とする。
actual退化因子 \(D\ne0\) に対し、次の部分集合を \(\mathcal L\) と書く：

\[
 \exists[\ell]\in(\mathbf P^4)^*,\qquad
 \ell^2\mid D,\qquad \operatorname{rank}N_\ell\le1.
 \tag{1}
\]

ここで canonical divergence-free lift を \(\bar\Pi\)、Euler場をEとし、
column-sharp規約で \(N_\ell=5\bar\Pi^\sharp(d\ell)/\ell\) とする。
各D因子はPoissonであり、\(N_\ell\) は一意なtrace零homogeneous linear場である。
`RANK_ONE_LONG_LINEAR_AND_PRODUCT_ENTRANCE.md` §§2–3 とその独立査読は、
marked factor incidence上でのこの割算のregular性を証明している。
dual-projective hyperplane parameterはproperなので、(1)は **\(\mathcal P_4\) 内で閉**。
Pf零点へ任意のfactorを付けたincidenceは使わない。

**定理。** \(\mathcal P_4\setminus\mathcal L\) は、全射影像の閉包次元が25以下である
有限個のconstructible sourceの和で覆われる。従って、この開集合と交わる
全Poisson成分の射影次元は25以下である。

この定理は\(\mathcal L\)内の全像に25上界を主張しない。実際、CA/CBの28成分はその中に入る。
全成分の個数上界又は全Eの点別既知被覆も、この定理だけの帰結ではない。

## 2. actual Dが被約、又は重複Hを持たない場合

以下で用いるのは**全像**の既受理定理である。各sourceには全plane/surface/frame、
全normalizer・全period/zero complementとtensor amplitudeを戻す。

| actual Dの全型 | 次元上界に用いる全sourceの結論 |
|---|---|
| normal | 全点T39、射影25 |
| integral非normal | `IRREDUCIBLE_NONNORMAL_QUINTIC_WHOLE_POINTWISE_COVERAGE.md` と独立COMPOSITION査読により、全点が有限既知閉和へ入る |
| reduced 1+4 | 全cofactor入口の既知閉和又は全射影像≤24 |
| reduced 1+1+3 | `TWO_HYPERPLANES_CUBIC_GENERIC_NONDOMINANCE.md` の全source：MPS又は全射影像≤24 |
| reduced 1+2+2 | `TWO_QUADRICS_SUPPORT.md` の全normal場・全pure-nil源：既知閉和又は全射影像≤24 |
| reduced 2+3 | 全MPSと六23成分 |
| reduced 1+1+1+2 | 全配置の有限既知閉和 |
| reduced 1+1+1+1+1 | 全T01又はrank四なし |
| nonreduced、重複Hなし | 正確にH+2Q、Q prime、rankQ≥3。全点既知閉和、最大25以下 |

reduced可約三型の「全像≤24」は、
`GLOBAL_GENERIC_G_EXHAUSTION_CLOSURE_AUDIT.md` §1末尾、§2および各個別source証明の意味そのものを使う。
固定βの像上界だけを指していない。1+4の最高quartic-openの全源を使うのはこの行だけであり、
後のmultiple-H表へそれを転用しない。

integral非normalの現在のpointwise定理では、quadric-coneの欠落を修復したV2、
rank零scrollの修正版全L6、全A25、全planeの全m/全rootを戻した合成だけを用いる。
撤回済みquadric V1、plus-sign triple coverage、旧G-only上界を所属定理に戻さない。

既知閉和の各rank四成分は、CA/CBの28を除けば射影25以下である。
この既知次元在庫は固定prior129全表と、その後に別途受理された新25/23成分の証明により確認される。
CA/CBの全rank四閉包はそれぞれproperな因子条件
\(D=H^3L_1L_2\)、\(D=HL_1^2L_2^2\)を持ち、generic marked labelで\(N\)rank≤1である。
proper factor incidenceでその条件を閉包へ延長したため、どちらも\(\mathcal L\)内である。
generic rank二以下の既知成分は閉rank条件によりここへ寄与しない。
従ってこの表の\(\mathcal L\)外の全像は25以下。

## 3. 重複H：全affine入口と除く低rank行

Hを無限遠へ移し、唯一のH-Casimir quadratic liftを
\(P=C+L+Q\)、\(\operatorname{div}P=c+Ax\) とする。
ここでは原点がsymplecticという仮定を置かない。
重複Hは正確に\(\operatorname{Pf}Q=0\)と同値であり、全元Jacobi、
\([c+Ax,P]=0\)、\(Q\wedge Ax=0\)を保持する。
これは `REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` §3 のH-Casimir gaugeの全入口である。
固定J/primitive βの後段を使わず、全C6/L24/Q60を残した表示を使う。

\[
 N_H=\begin{pmatrix}0&0\\c&A\end{pmatrix},\qquad
 \operatorname{rank}N_H=\operatorname{rank}[c\ A].
 \tag{2}
\]

A零、又はrankA一でcがimage内の場合は(2)がrank一以下である。
従ってこの全行は\(\mathcal L\)内であり、source次元を上界評価する必要がない。
ここでA零の巨大unimodular源やCA/CBを、小さい部分sourceで置き換えていない。

残る全型はrankA≥二、またはrankA一・c非imageである。
traceA零なので、rank二の非nilpotent型はdiag(λ,−λ,0,0)、λ≠0だけであり、
nilpotent型は(3,1)/(2,2)。rank一はnilpotent(2,1,1)だけ。
rank三nilpotentは(4)だけである。cをimageで割った後の全coker-orbitも下表に保持する。

## 4. rankA≥三とrank二非nilpotent：全点MPS又は小さい全source

`MULTIPLE_HYPERPLANE_MODULAR_SUPPORT.md` §§2–6 と入口の §§6–8 による全MPS被覆は、
次の**全行**に適用する。

* rank四regular SS：非共鳴、pair-zero共鳴、pair-zeroなし共鳴。
* rank四repeated-eigenvalue SSと全非semisimple。
* rank三のsimple-zero：regular SS、double-eigenvalue SS、nonsemisimpleを含む。
* rank三のspectrum(a,−a,0,0)、zero Jordan block。

この列挙はpotential SNC点の存在だけではなく、全低次項を保持したactual polynomial pencilと
有限MPS閉和の閉性による全所属を使う。全零period、対和零、非零translation、
linear共鳴係数が零になる補集合を、generic diagonal開の外へ置き去りにしない。
rank二非nilpotentの全c=image/非imageは
`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §2 の全g零/定数/非定数・全resonance
actual pencilによって同じMPS閉和へ入る。
MPS40の各全成分は射影25以下なので、これら全行の像は25以下。

rank三nilpotent(4)は `MULTIPLE_HYPERPLANE_NILPOTENT_SUPPORT.md` の全sourceを使う。
c非imageの全sourceはPf零。c=imageのrank四全点は五つの固定source各四変数に尽くされる。
fieldのAff4 orbitは15なので、全affine像≤19、H4とprojective scaleを戻して≤22。
式の零境界から別の未処理源は出ない。

## 5. rank二nilpotent：quartic-open表でない全五型

全入口は `CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json` の
全90係数reductionと全元Jacobi/Pfである。
`HYPERPLANE_QUARTIC_NILPOTENT_*` のtype5/6/7をこの表の証拠に使わない。
この表では最高Pf零、すなわちmarked Hの重複を最初から保持する。

| A,cの全affine型 | 固定fieldでの全source上界 | Aff4 field-orbit | 全射影像上界 |
|---|---:|---:|---:|
| (3,1), c=image | 三canonical source、各10 | 12 | 25 |
| (3,1), c=長いchain top | 六source、各5 | 14 | 22 |
| (3,1), c=短いchain top | 六source、各6 | 13 | 22 |
| (2,2), c=image | 全十一canonical源から≤12、§6 | 10 | 25 |
| (2,2), c非image | 八source、各7 | 12 | 22 |

(3,1) image行は `MODULAR_NIL31_COMPONENT.md` のn7 unit十変数と、n7零に残る
half/thirdの全十変数二枝を全て数える。n7=0でn14等にunit行の条件を押し付けない。
long/shortはそれぞれ `MODULAR_NIL31_LONG_SUPPORT.md`、`MODULAR_NIL31_SHORT_SUPPORT.md`、
最後のtop行は `MODULAR_NIL22_TOP_SUPPORT.md` の全残式からの有限source定理を使う。
各rank四のunit条件は正確なPf式で判定され、unit零は別行又はPf零へ戻る。
個別のE帰属が未完でも、この全像上界はその完成を仮定しない。

orbit数はGL4 centralizer次元6/8、translation image二、およびcoker上の
long二/short一/(2,2)非零二の実orbitから得る。
全source affine像≤(固定source上界)+(field orbit)というfield写像のfiber上界を使い、
source stabilizerを余計に二重に引かない。

非零scalar倍は各nilpotent affine型をactual affine共役で同じ型へ戻す。
chainの座標に逐次整数重みを選び、Nの各矢印を同じscalar倍、残cのtopも同じscalar倍に出来る。
従ってこの全affine像はconeであり、projective scaleにより一次元引ける。
marked Hを動かす四次元を足し、表の最終列はsource+orbit+4−1である。
これらはsymplectic値の固定や三つのTaylor fiberを引いたβ上界ではない。

## 6. nil22 image：全十一源のCstar fiberによる≤12

`NIL22_IMAGE_SOURCES_AND_COMPONENTS.md` §§1–5、全90から独立に復元した
`NIL22_IMAGE_SOURCE_REVIEW.md` §§1–8 の必要十分な全表示を使う。
N=p1∂q1+p2∂q2で、全Pは

\[
 P=F\partial_{q_1q_2}+\sum M_{ij}(p)\partial_{q_i p_j},\qquad
 M=-\operatorname{adj}(S+\operatorname{Hess}D_3/2)+p\otimes(V_0+V_1p),
\]

であり、Fのq係数はdivから一意、残F0(p)は任意二次以下の全六係数。
同時変換(q',p')=(Aq,Ap)、A∈GL2はNを保つ。全変換則は

\[
 D'_3=\det(A)^2D_3(A^{-1}p'),\quad
 S'=\det(A)^2A^{-t}SA^{-1},\quad
 V'=AV(A^{-1}p'),\quad F'_0=\det(A)F_0(A^{-1}p').
 \tag{3}
\]

binary cubic三root型/零、Hessian全rankの全補集合を含む十一canonical源は、
六磁気係数込みで順に

```
stable_line_conic 6, stable_three_lines 6,
double_A 6, double_B 6, double_C 6, double_D 8,
triple_A 6, triple_B 6, triple_C 9,
quadratic_A 7, quadratic_B 9.
```

GL2四次元を復元すると、9変数二行以外は直ちに固定field像≤12。
残二行も、以下の**全parameterを保つactual Cstar作用**を持つ：

* triple_C：D3=x³,S=0,V0=(0,n),V1=[[-1/3,0],[c,d]]。
  A_r=diag(r²,r)、r≠0でD3とSを保ち、(n,c,d)↦(rn,c/r,d)。
* quadratic_B：D3=0,S=diag(1,0),V0=(0,n),V1=[[-1/2,0],[c,d]]。
  A_r=diag(r,1)でD3とSを保ち、(n,c,d)↦(n,c/r,d)。

F0は(3)により六次元全空間を可逆に動く。どのn,c,dを零にも特殊化せず、unit条件も増やさない。
GL2×Sからのtensor写像の各fiberは
\((A,s)\mapsto(A A_r^{-1},A_r\cdot s)\) の一次元orbitを含む。
第一factor上で作用はfreeなので、tensorにCstar安定化群が有るとの仮定は不要。
従って両9変数行の像は4+9−1=12以下。

D3=0,S=0だけはD=0でrank四の外。全十一源の有限被覆からfixed N全sourceの次元≤12を得る。
これは旧 `MODULAR_NIL_RANK2_FIRSTJET_BOUND.md` の安全な14上界を純粋に改善するが、
旧β証明/旧数値証明書は変更しない。新しい成分帰属もこの計算には使っていない。

## 7. rank一・c非imageの全long/short

`RANK_ONE_TRANSLATION_WHOLE_FINITE_EXHAUSTION.md` と独立FINITE_EXHAUSTION査読の
十七whole源を用いる。非定数g一、定数g/H3零一、定数g/H3非零二、g零十三の全源であり、
各固定long sourceの次元は12以下。全affine U、全F6、全H低次/三次、全λと全zero/unit復帰を含む。
全caseの成分性を証明済みと仮定してはいない。

long field X=∂x+x∂yのAff4 orbitは10、全sourceはscalar-stableなので、
H4とscale−1を戻して各全射影像≤12+10+4−1=25。
short全sourceは `RANK_ONE_SHORT_TO_LONG_WHOLE_POLYNOMIAL_CURVE.md` とその独立査読の
actual polynomial familyにより、同じg枝のlong源閉包の有限和へ入る。
従ってshortの全定数、全g零、全cubic resonanceも≤25に入る。
actual affine/PGL mapでないr=y−x²/2を成分同型に使わない。

## 8. 有限在庫から全成分の上界へ

actual quinticの因子分割は有限。重複primeのdegreeは一又は二だから、
§§2–3は全rank四tensorを尽くす。重複Hを持つ場合、rankAと全Jordan/coker型は§§3–7で尽くされる。
spectral parameterを持つ非nilpotent行は、全periodを保持した有限MPS閉和へ既に送られており、
新たに無限個のnormal form閉包を足していない。

従って\(\mathcal P_4\setminus\mathcal L\)は、次元≤25の有限個の閉像の和に入る。
この開集合と交わる既約成分Cは、その開部分が有限閉和に入るため、そのどれかの閉像に含まれる。
よってdim C≤25。ここではGの必要成分次元25–27、βの優越、点の接次元、
一般log型の三つの五点を一切使っていない。

## 9. M1/M11の全25像への適用

既独立純入口 `RANK_ONE_LONG_LINEAR_AND_PRODUCT_ENTRANCE.md` は、全sourceからのinjectivity、
Stab(X)=10、全action fiber、scalar stability、一般marked Hの一意性により、
Z1/M1とZ11/M11のPGL飽和像閉包が正確に25次元の既約閉集合であることを証明している。

それぞれの一般Dは4h*L、3h*L1*L2で、hは唯一の重複hyperplane。
そのN_h=h∂x+x∂yはrank二であるため、どちらの一般元も\(\mathcal L\)の外にある。
Zを含む全空間既約成分Cを取る。Cはrank四開と\(\mathcal L\)外に交わるので、
本定理よりdim C≤25=dim Z。従ってC=Z。

**帰結：M1/M11の両全25像閉包は、相異なる二つの全空間既約成分である。**
相互および既採録成分との相違は既入口の全marked-rank/factor監査を用いる。
両源の全boundaryも閉包に入るため、この成分性はwhole sourceの被覆に使用出来る。

最初の一点における接次元28/27および真正quadratic coneの正次元性は不変である。
滑らかさ、局所成分一意性、全schemeの非被約構造、横断multiplicityはこの純次元証明からは従わない。
新しい点、より高いjet、別normal sliceは実行していない。台帳採録と独立/root受理は統合担当が行う。

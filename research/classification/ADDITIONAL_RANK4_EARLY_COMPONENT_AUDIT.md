# 追加rank4成分の早期五族：被約成分性と相違

2026-11-07。分類義務AのDOUBLE_SMOOTH_QUADRIC、TRIPLE_HYPERPLANE_QUADRIC、
CA、CB、ZERO_LOG_5Hを扱う。全射影Jacobiスキームの被約既約成分についての監査である。

**結論：下記の既存全係数証明書を入力として、五族の全空間成分性・次元・相違を受理する。**
旧行列をこの回に再実行したとはしない。真の族の微分と全周囲方程式を接続する論理、
および一般の因子型から閉包全体の識別へ移る論理を確認した。
新たな成分数を加算するものではない。

## 1. 滑らかな証拠点を持つ四族

全二次斉次liftの150次元空間で、Dを全350式の射影Jacobi微分、Vを24次元の
鉛直部分、Mを鉛直・PGL作用・スカラー方向、Fを実際の族の全微分を含む行列とする。
各点でD F=0が厳密な係数恒等式として検査されている。
実際の既約像閉包の次元と局所次元は

\[
\operatorname{rank}F-25\le\dim\overline{\mathrm{image}}
\le\dim_{[P]}\mathcal P\le125-\operatorname{rank}D
\]

で挟まれる。下表では両端が一致するので局所環は正則で、実際の像閉包が全空間成分である。
固定退化因子や固定補助スライスの接空間を使っていない。

| 族 | rank D | rank M | rank F | 射影次元 |
|---|---:|---:|---:|---:|
| DOUBLE_SMOOTH_QUADRIC | 103 | 47 | 47 | 22 |
| TRIPLE_HYPERPLANE_QUADRIC | 103 | 47 | 47 | 22 |
| CA | 97 | 47 | 53 | 28 |
| CB | 97 | 47 | 53 | 28 |

前二族は実軌道閉包である。CA、CBはbinary quadratic fと四変数三次potential hを
用いた全19パラメータPoisson族をPGLで動かす。Jacobi恒等式は全パラメータ上で成立する。
したがってfamily微分は仮想変形方向でない。
有理数体上の非零minorとD F=0による階数の挟み撃ちが、有限体階数から
根拠なく標数零の上界を推測する問題も避けている。

## 2. ZERO_LOG_5Hの成分性に必要な弱い局所命題

元の二パラメータ族P(a,b)=h²J+Q(a,b)は全(a,b)でPoissonであり、退化式はh⁵である。
(a,b)=(5,7)では全Dの階数101、真の族・群・鉛直・スカラーの微分階数48。
一つの余剰法線Rを加えると49となり、全接空間を張る。
以下のB(U,V)は、JacobiのpolarizationにEuler場を外積した全350式の値を表す。
全方程式で二次補正Tが存在する一方、

\[
\operatorname{rank}[D\mid B(R,\ker D)]=101,\qquad
\operatorname{rank}[D\mid B(R,\ker D)\mid B(R,T)]=102
\]

である。従って二次補正の任意の変更を許しても、三次残余は微分像で消えない。

`ADDITIONAL_RANK4_LATE_COMPONENT_AUDIT.md` 第1節の形式横断切片補題を使う。
48個の真の族方向を形式基底座標とし、Jacobi線形部の101方向を陰関数で消去した後、
族方向を零と置いた真の切片は一法線の形式環である。
上の補正曖昧さを含む障害が、切片の少なくとも一つの残余式を非零にする。
一変数形式環の非零idealの商は零次元なので、元の局所次元は48以下となる。
真の族の48次元下界と一致し、射影化と鉛直商の25次元を引いて23次元成分を得る。

これは成分性に必要な支持論だけを取り出したものである。
全族上の普遍二次jetや正確な一般横断環C[ε]/(ε³)の決定を、主成分数の追加前提にしない。
元の強い局所環定理の証拠は保持するが、その精緻化を現在の分類のために再開しない。

## 3. 相違を全閉包で保証する幾何

rank4開上の写像[P]↦[P∧P]は退化五次式へのmorphismである。
Lq²、L³q、L⁵という形の五次式の集合は、それぞれ射影的な掛け算写像または
Veronese写像の像なので閉である。零式を含むaffine coneを引き戻せば、
rank2境界を含む全Poisson空間上でも閉条件として使える。

| 族 | 一般の退化式 | 閉包全体で使う条件 |
|---|---|---|
| DOUBLE_SMOOTH_QUADRIC | Lq²、rank q=5 | 平方の二次因子を持つ |
| TRIPLE_HYPERPLANE_QUADRIC | L³q、rank q=4 | 三乗の一次因子を持つ |
| CA | L³L_1L_2 | 三乗の一次因子を持つ |
| CB | LL_1²L_2² | 平方の二次因子を持つ |
| ZERO_LOG_5H | L⁵ | 一次式の五乗である |

全て被約SNC因子を排除するため、五族の閉包はMPSのSNC40成分と一致しない。
一般rank4により六rank2成分とも異なる。

五族内では次元22、28、23で大半を区別できる。
22次元二族では、滑らかなqのLq²は三乗一次因子を持たず、L³qの閉条件と異なる。
28次元二族では、一般CBの重複度(1,2,2)は三乗一次因子を持たず、CAの閉条件と異なる。
ZERO_LOGの23次元はこの四族と異なり、後期の23次元族とも一般因子型L⁵で異なる。
後期六族の比較は `ADDITIONAL_RANK4_LATE_COMPONENT_AUDIT.md` に接続する。

「一般点が非SNC」であることだけでMPS40との差を結論していない。
相違を示す閉条件が全像閉包に保たれることが必要であり、ここではそれを確認した。

## 4. 証拠の入口と全支持への追加の利用

| 族 | 原証明 | 独立有限証拠 |
|---|---|---|
| DOUBLE_SMOOTH_QUADRIC | `legacy/v03/p4_poisson_completeness_v3.tex`, prop:tangent | `legacy/v03/data/rigid_quadric_independent_audit.json` |
| TRIPLE_HYPERPLANE_QUADRIC | `legacy/v05/triple_hyperplane.tex` | `legacy/v05/data/independent_triple_hyperplane_audit.json` |
| CA,CB | `legacy/v06/ambient_families.tex`, thm:AB | `legacy/v06/data/independent_ambient_families_audit.json` |
| ZERO_LOG_5H | `legacy/v06/zero_log_component.tex`, thm:zero_log | `legacy/v06/data/zero_obstruction_independent_audit.json` |

有効な原producerと入力のSHA256は
`reports/CLASSIFICATION_ADDITIONAL_RANK4_EARLY_DEPENDENCIES_20261107.json`に保存した。
これは今回の新しい350×150行列演算の報告ではなく、依存監査の記録である。

DOUBLE_SMOOTH_QUADRICについては原稿のthm:quadricがさらに、
実際の退化因子H+2Q、Q滑らか、を持つ**全点**をこの一つの軌道で尽くしている。
この認識は接空間minorだけからは出ず、低次数有理二形式の全90係数空間と
conformal symplectic/orthogonal条件の議論に依存する。
この幾何的支持命題とrank4 quadricへの拡張は、別の現行二重二次因子ノートで扱う。

本監査は既知十一族の成分性・相違を主分類に使える形へ弱めるもので、
十一族以外の未知成分の不存在や一般S61ファイバーの尽くしは証明しない。

## 5. 独立査読

2026-11-07、後期六族担当による独立査読。原証明四件、上記独立証明書、
ZERO_LOGの `legacy/v06/code/audit_zero_obstruction.py` を照合し、
被約成分性・次元・五族の相違の結論を受理する。ZERO_LOGでは滑らかな実際の族源から
48方向を選べる。真の切片の暗黙補空間が選ぶ二次補正T'はTとの差がker Dに入り、
したがって三次係数B(R,T')の余核像は証明書の非零像B(R,T)と一致する。
さらに高次補正の三次への寄与はDの像なので、この一変数残余式は非零である。
これで零次元切片とKrullの次元上界が得られ、普遍jetや正確な横断長を追加仮定しない。
各因数分解条件の射影的閉性、一般点の重複度型、次元22・28・23による相違にも
論理上の不足は認められなかった。上でBの全350式への射影を明記したほか、
結論の修正は不要。旧行列の再実行は行っていない。

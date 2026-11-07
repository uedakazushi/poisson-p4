# regular・linearly normalな曲面上で零となる全bivector：ACMを使わない入口

2026-11-07、主義務G。S⊂P4をintegral・非退化なprojective surface、
`E=Λ²T_P4`、`r=h⁰(I_S(2))`、`l=h⁰(O_S(1))−5`、`q=h¹(O_S)` とする。
Sのsmoothness、normality、ACM性を仮定しない。

## 1. 全sectionに対するEuler upper bound

global外積Euler列とH¹(T_P4)=0により、全global bivectorは
全150二次係数のlift Πを持つ。
π|S=0なら、restricted外積Euler列のkernelから一意の
`τ∈H⁰(S,T_P4|S)` を得る。
τがambient H⁰(T_P4)へliftすれば、Euler∧Zを引いたΠの全係数はI_S(2)に入る。
従ってsnake argumentから

`0→I_S(2)⊗Λ²C⁵→H⁰(E⊗I_S)→coker{H⁰(T_P4)→H⁰(T_P4|S)}`

がexactである。
左端のinjectivityは、ambient projective linear fieldがSを点ごとに固定すれば零であることによる。
実際、そのmatrixのeigenvector locusは有限個のprojective eigenspaceの和であり、
integralなSは一つに入り、SがP4をspanするのでmatrixはscalarとなる。

restricted Euler列は

`0→O_S→O_S(1)⁵→T_P4|S→0`。

従って `h⁰(T_P4|S)≤5h⁰(O_S(1))−1+q`。
ambient H⁰(T_P4)はdimension24で、上の非退化性からrestrictionはinjective。
よって全linear sourceの上界は

`h⁰(E⊗I_S)≤10r+5l+q`。                                  (1)

これはPoisson条件より前の全sectionの評価であり、chosen pointやresolutionのBetti数に依存しない。

## 2. q=0かつlinearly normalなら完全quadratic sourceとなる

q=l=0ならτは必ずambient linear fieldへliftするので、上記は同型

`H⁰(E⊗I_S)=I_S(2)⊗Λ²C⁵`。                              (2)

したがって全sourceは `Π=Σ_j q_j A_j`、A_jは全constant skew matrixで尽くされる。
Euler商の未記録kernel又は「ACMなら」という追加仮定は無い。
この結果は全sourceの完全性であり、Poisson支持の方程式は必要な場合に全て保持する。

## 3. degree五以上の全rank零枝を排除

deg S≥5ならr≤1である。
もし二つの独立quadricsがSを含めば、各quadricはirreducibleである
（reducibleならprime I_Sにlinear factorが入り、Sの非退化性に矛盾）。
互いに素な二quadricのcomplete intersectionはdegree四なので、deg S≥5に矛盾する。

q=l=0では、r=0なら全bivectorが零、r=1なら `Π=q₂A` だけである。
後者の実anticanonical五次式は

`D_Π=q₂² Σ_i (−1)^i x_i Pf(A_without_i)`。

従って零又は非被約であり、既約・被約五次退化因子にならない。
この排除はJacobiを課す前に成立する。

**結論。** integral・非退化・regular（q=0）・linearly normalなdegree五以上の曲面は、
既約被約quinticを持つrank四Poisson構造のrank零曲面になれない。
deg S≤6という既存normalization-length公式と合わせれば、degree五・六の全該当枝を閉じる。
non-ACMなrational surfaceも(2)の仮定を満たす限り含まれる。

Sがnormalでresolutionのqが零なら、Lerayの単射
`H¹(O_S)→H¹(O_S̃)` によりq(S)=0を得る。
さらに `INTEGRAL_SURFACE_LINEAR_NORMALITY.md` で、特異Severi定理の射程を一次資料で確認し、
integral・非退化なS⊂P4はsmooth projected Veronese以外linearly normalであることを証明した。
この唯一の例外はdegree四なので、上のdegree五以上の結論ではl=0を追加仮定する必要が無い。
従ってnormalityやACM性に依存せず、**degree五以上・q=0の全integral非退化曲面**を排除する。

## 4. 残す例外と既済例との整合

smooth projected Veroneseはq=0だがl=1であり、(1)は五次元の可能性を保持する。
既存の完全rank零sourceが五次元であることと整合し、ここで誤って除外していない。
stable degree五楕円scrollはq=1、linearly normalであり、irregularな例外を保持する。
その既知T40 bivectorの認識は別証明である。

qが非零の例外及びdegree四以下の別軌道のPoisson支持分類は、本稿のEuler評価だけからは従わない。
linear nonnormalityについては上記の別補題が全integral曲面を扱い、projected Veroneseだけを残す。
同補題には特異曲面を含むSeveri定理という明記された外部入力があり、
smooth版をnormal又はintegral曲面に拡張して用いたわけではない。
既約非正規quintic全体やG全体の完了は主張しない。

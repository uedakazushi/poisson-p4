# 全Qa共通の三次Darboux/Cremona写像：独立Fraction監査

2026年11月3日。root作成の凍結packet `Q_QA_UNIFORM_DARBOUX_20261103_v63` を監査する。
対象producerは変更・importしない。

独立実装 `q_Qa_uniform_darboux_review_20261103_v63.py` は、
QQ[a,h,x,y,z,w]のmonomial辞書とPython fractions.Fractionだけで
和・積・微分・全多項式composition・Poisson bracketを実装する。
SymPyは保存された表示文字列を有理係数へ読むためだけに使い、
対象producerのmatrix微分・tensor積の検算を繰り返す方式ではない。
新規producer replayの全certificate byteとruntime以外の全report項目も一致させる。

原物理Qaを四つの非零homogeneous bracketから独立に組み立てる。
同じ三次shear

    X=x, Y=y, Z=z,
    W=w+6xz−3y²/2−9x²y/2

と逆を自身の係数辞書で作り、保存したF,Gの全係数を照合する。
両方向全5座標のcompositionは恒等座標のh8倍であり、
全10 bracketを両方向、計20本で直接比較して、raw homogeneous Poisson恒等式を確認する。

affine h=1でも全6 bracket、両inverse、最後のw微分が1であるtriangular Jacobian、
aに座標写像が依存しないことを別に検査する。
標準形は厳密に `a*∂X∧∂W−(2/9)∂Y∧∂Z` である。
どの演算にもaを逆にする分母は現れない。

F0=G0=h³なので共通因子はhだけが候補となるが、
両最後の座標のh=0での固定項はそれぞれ−9x²y/2、+9x²y/2である。
従って全parameter特殊化でgcd=1、bidegreeは正確に(3,3)。
三次写像の最小性やPGL同値性は主張しない。

a≠0で標準tensorは至る所rank4であり、連結なA4開全体が一つのsymplectic leafになる。
a=0では標準形のPoisson centerはC[X,W]で、各fiberは(Y,Z)のA2となる。
非零定数−2/9のbracketがあるため各fiber全体が一つのleafである。
独立実装は原Qa0のtensorでもx,Wと全4座標のbracketが零になることを直接確認する。

これはh≠0のaffine葉の完全記述であり、h=0での全境界葉、
全Poisson局所枝、既知成分への帰属や尽くし上げには拡張しない。

保存：

* `code/q_Qa_uniform_darboux_review_20261103_v63.py`
* `reports/Q_QA_UNIFORM_DARBOUX_20261103_v63_REVIEW.json`

全入力SHA256を最終receiptへ保存する。

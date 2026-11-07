# m4 two-line highrank211全空性の独立review

2026-11-07。別作者の全数学的読解とexact Fraction/RF算術。**PASS**。
実D被約既約quintic、rank二plane沿い位数四、generic normal quartic211、leading image rank二以上、
zero cubic R のsupportが二直線以下である全sourceは空である。
m2/m3、他支持又は全G/Eの完了を主張しない。

## 1. 全leading geometryと全actualD正常形

producerの `PLANE_M4_TWO_LINE_COLLISION_QUINTICS.md` と
`PLANE_M4_TWO_LINE_COLLISION_EXCLUSION.md` を全読解した。
既受理のfixed-double-root及びregular symmetric2-pencil分類を本射程へ適用出来ることを確認した。
rank三netはirreducible collision conicを要求し、support≤2linesには入らない。
semisimple unaligned marked rootは三つの相異used linesを要求するため除かれる。
残るaligned sem及びJordan型は二used linesを要求するのでR=3L高rankは空。
R=2L+Mの両orderingを固定せず、全independent A,B・全nonzero振幅を保つ。

semの全actualDはu4A+u2v2B+g5。
A/Bへのnormal-linear shearsはu5,u4v,u3v2,u2v3を全て除去する可逆GL5変換である。
残るαuv4+κv5でκ零ならDにu因子があるため、既約性はκ≠0を強制する。
nonzero normal/base/scalar補償でκ=1、α=0又は1へ戻る。両補集合を保持した。

Jordanの全actualDはu2v2A+u3vB+g5。
同じ全shearsでu4v,u3v2,u2v3を除去し、残るγu5+αuv4+κv5の両endcoefficientsは非零。
γ零はv因子、κ零はu因子を与えるからである。
両endcoefficientsを一へ戻した後、αを全complex τとして残すことが必要。
有限fifth-root ambiguityはあってもτを0/1へ過剰正規化していない。
全shear逆・尺度逆は全actualD源へ戻せる。

したがって必要な全quinticsは

\[
 F_{s,\epsilon}=u^4A+u^2v^2B+\epsilon uv^4+v^5\ (\epsilon=0,1),\qquad
 F_{j,\tau}=u^2v^2A+u^3vB+u^5+\tau uv^4+v^5\ (\tau\in\mathbb C).
\]

各式はAに関してprimitive linearであり、残るconstant termとA-coefficientのgcdは一。
ゆえにC上絶対既約。四gradientのconstant independenceからexactly essential-four。
これらはPoisson存在の仮定ではなく、全actualDを覆う必要正常形である。

## 2. full150源を別Fraction/RFで再構成

canonical divergencezero strict liftとPfaffian-Casimir identityを使い、
全源をL={quadratic Π:divΠ=0,Π#dF=0}に含めることは正当。
このidentityはessential-threeに限定されず、任意Euler代表のdivergenceを使う誤りとも異なる。

新checkerはproducer、producer helper、CASをimportしない。
既独立sparse Fraction engineとFraction(τ)算術だけを再利用し、各quinticのpartial derivativesから
全1050 degree六Casimir行とdivの全25行、全150quadratic係数をゼロから構成した。
producerの全1075×150 entry・basis・degree六monomial順と一致した。
vertex wを含むmonomialも最初から全て保持している。

| 全quintic | fullrank | 完全kernel dim | 独立検証するparameter範囲 |
|---|---:|---:|---|
| sem ε0 | 148 | 2 | 固定全源 |
| sem ε1 | 148 | 2 | 固定全源 |
| Jordan | 147 | 3 | QQ(τ)、全τ≠0 |
| Jordan τ0 | 147 | 3 | 原全matrixのτ=0 |

全maximal determinantsを独立に再計算した。
Jordan generic determinantは正確に1700051933833072276930560000000000 τ6で、
唯一のzeroはτ0。そのfull complementは別のnonzero147minorで覆われる。
generic kernelの全denominatorsもτのpowersだけであり、未検査のcomplex wallは無い。
各保存kernelについて全1075row代入とfree-coordinate identity minorを検証した。
nonzero maximal minorとindependent fullkernelの次元和が150なので、全kernelであることも確認した。
chosen generic pointのrankを特殊τへ外挿していない。

## 3. 共通vertex factorによる全rank四空性

全六base4×base4 pairに属する九十quadratic係数を、各fullkernelの全columnで検査した。
固定二源、generic Jordan、τ0の全てで九十行はliteral零。
したがって全必要源は

\[
 \Pi=\partial_w\wedge V(A,B,u,v,w)
\]

を満たす。w-dependenceを捨てたsliceではなく、全150係数から得た等式である。
rawrank≤2、projective quotientでもrank≤2。
rank四Poisson tensorを持つことは出来ず、全normalization逆でもこの矛盾は保たれる。
Jacobi、actualPfaffianのnonlinear eliminationや追加ambient tangent計算は不要。
producerの弱いdimension上界19/20/21を使わず、さらに弱く十分な共通因子支持だけで閉じた。

## 4. 独立receiptと残余

finite契約 `jobs/plane_m4_two_line_casimir_review.json` をgate後に一回実行し、
0.718079秒、RSS49,992KiB、全exact独立算術PASS。
gateそのものを数学の認証とはしない。

| artifact | SHA256 |
|---|---|
| code/classification_plane_quintic_casimir_review.py | 6029d24d0fba17ca9c50d8042217316b9922193a2ba9bea3f07a8aa0f46b5fe0 |
| reports/CLASSIFICATION_PLANE_M4_TWO_LINE_CASIMIR_REVIEW_20261107.json | b28507ba5334a4907a6491c02df7ffa22d7e4e32f38cd0beb7963ee2df3b48b2 |
| research/classification/jobs/plane_m4_two_line_casimir_review.json | b694f7c2255d11ae6c50d8a76871205943a39acb67781a4a3d57e4551d2a5787 |

全executed reviewer bytesと入力hashはreportへ保存した。
旧producer/report/currentjobs及び以前の全証明packetは変更していない。
nonredRのm4squarefree、他leading型、m2/m3等はそれぞれの独立証明と合わせて初めて全体を論じる。
本receiptだけで全m4/G/Eを宣言しない。

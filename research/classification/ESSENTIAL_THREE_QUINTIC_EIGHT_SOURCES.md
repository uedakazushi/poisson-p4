# 本質三変数・既約五次の全八actual Poisson源

2026-11-07、PRIMARY E。今回のproducer算術は全八源でPASSした。
全曲線の必要分類は `ESSENTIAL_THREE_QUINTIC_FINITE_CURVES.md` と参照する全linear-series証明による。
全canonical Saito定理・全磁気inverseは既に独立review済みである。
新しい八源の算術と全曲線分類の独立reviewは別receiptで確認する。
八個の全空間成分や新しい八成分を主張しない。

## 全閉点に対する有限入口

一般rank四Poisson構造で実退化五次Dが既約被約、本質変数が三個である全支持は、
以下八つのactual PGL5軌道の和になる。
各軌道は23次元であり、その有限個の閉包が対象全支持を覆う。

| curve | 全二次核次元 | det(E,U,V)/F |
|---|---:|---:|
| 一尖点(4,5) |2|−9/64|
| A12 |2|225/104|
| E8+A4 |2|2/9|
| E6+A6 |2|−175/144|
| E6+A4+A2 |2|−125/48|
| A8+A4 |2|125/72|
| 3A4 |2|−5/8388608|
| A6+3A2 |2|−75/3584|

表の値は公開inputの具体有理parametrizationと、reportのprimitive Fおよび完全kernel basisに対する値である。
各primitive F、全二次基底、GL2 frameと全六磁気式はreportにそのまま保存した。
異なる行はcurveの尖点configuration又は最大重複度で区別され、同じ軌道を重複させていない。

## implicit曲線と複素絶対既約性

各degree5 binary parametrizationの三sectionは共通basepointを持たない。
全degree四の平面式十五係数を代入する係数mapは全rank十五、
全degree五の二十一係数を代入するmapはkernel一である。
全行列と非零最大minorを保存し、その一核からprimitive quintic Fを得た。
元parametrizationへの代入は恒等零である。

P1のimageはC上で既約である。
もしFがC上可約なら、そのimageを含むproper factorはdegree四以下なので、
degree四のrelationを与えて上の全rankに反する。
従ってこれはQQ上の因数分解だけの検査でなく、複素絶対既約性の証明である。
五次Fの三偏微分の全係数行列と非零3×3minorも保存し、全行でessential rank三を確認した。

## 全十八係数と全Poisson式

各Fについて、全quadratic gradient syzygy十八係数と全degree六二十八行を構成した。
全28×18行列の階数は16、kernelは全て二次元である。
二本U,VがFを消すこと、det(E,U,V)=cF、c≠0を直接確認した。
従ってSaito criterionから全gradient syzygy moduleの自由基底であり、選んだ部分源だけではない。

全GL2 frameの四係数と全磁気Rの六係数を残したまま

    Π=U'∧∂u+V'∧∂v+
      (−div(V')u+div(U')v+R)∂u∧∂v

を構成した。
全十strict Jacobi、全五divergence、全五Casimir contractionを恒等零と検査し、
全五raw Pfaffianと実D=−c det(frame)Fを保存した。
局所化は固定非零Dのrank四源に必要なdet(frame)だけで、全磁気零・非零は残る。
Hilbert–Burch resolutionから全τ=12も従う。

全canonical essential3 theoremにより任意のΠがこの表示に入り、
全磁気inverseによってRは実shearで零へ移る。
全frameも実vertex GL2で移るので、各固定curveの全源は一つのactual PGL5軌道である。
curveをPGL3で移す操作もPGL5へ延びる。

## 三A4の有理代表と黄金比canonical代表の同一性

全曲線分類では、三cuspsを0,∞,1へ置いた全sectionのA4条件から

    A²+A−1=0,  B=−A−1,  C=−A−2,  D=1

を得ており、二つのconjugate choicesはp,q交換で同値である。
算術inputではMoe p.90に印刷された別の有理parametrizationを使った。
原PDF画像と係数を照合し、同じconfigurationとの記述だけを同一性の根拠にはしていない。

その有理mapのimplicit Fは絶対既約で、全18核が上記の通りSaito(2,2)を満たす。
二つのbinary微分vectorのcross productの全gcdは

    s³+4s²t−16st²/5

で、複素projective P1上三つの相異なる単根を持つ。
一般にbasepoint-free birational parametrizationでは、このgcdの各根の重複度は
image branchのmultiplicityから一を引いた値である。
従って三つの尖点は全て重複度二である。
free irreducible quintic定理は全特異点がcuspidalであることも保証する。
全必要分類で三つの重複度二尖点を持つ候補は三A4だけなので、
この有理Fはその唯一のclassに入り、黄金比canonical代表とPGL3同値である。
generic点の一致や特殊代表の取り違えによる推測ではない。

他の七入力は全曲線証明のcanonical representativeであり、
全binary differential gcdのcomplex root multiplicitiesも公開reportで一致を確認した。

## 証明書・停止範囲

- runner: `code/classification_plane_quintic_saito.py`、isolated `--eight-curves-input` mode。
- input: `inputs/essential_three_quintic_eight_curves_20261107.json`。
- contract: `jobs/essential_three_quintic_eight_sources.json`。
- report: `CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_SOURCES_20261107.json`。

最終算術は約2.588秒、RSS約62MiBで全証明書の保存までPASSした。
初回は同じ八入力の算術後にJSON型変換の保存エラーとなり、証明書は作られなかった。
そのfailed runを証明とは呼ばず、整数serializationだけを直して新gate後に再実行した。
以前の全functionはmain以外AST不変、旧report bytesと旧四契約は変更していない。

最初の(4,5)代表の全源は既採用の
`R4_PLANE_UNICUSPIDAL_QUINTIC_45_23`へ入る。
他の七軌道閉包が全成分か、より大きな既知成分の境界かは本jobで計算していない。
各残源に一つだけの全周囲・最弱支持又はknown-parent認識を接続することが次のE義務である。
新成分数を七又は八に固定せず、全P4分類やG一般五点の尽くしの完了とも呼ばない。

# generic-a cubic36 診断：入力・停止境界の独立監査

2026-11-04。対象は frozen diagnostic とその read-only preflight。
Singular、GB、std、field job、元行列生成は一切再実行しない。

## 実際に検査する算術

元 `after21` の SHA、36×70 header（第三整数21は完了pivot数）、使用21列・
残存36行の順序を V6 manifest と phase21 に結ぶ。2520 polynomial records を
全て読み、421667個の既約QQ係数・非負指数・厳密な降順を調べる。
各有理係数を32003で独立に還元し、分母がすべて単元、消える分子項が14、
残る項が421653であることを確認する。

一方、保存された Singular input を独自の限定文法parserで読む。a の各係数
多項式と e,d monomial を戻し、**全36×70 entry の全係数**を元の還元辞書と
比較する。producerの文字列生成・grouping関数やSingularは使わない。
書かれた非零entryは1631。行と列を入れ替えた検査ではない。

ring は `F32003(a)[e,d,z]`、column map は R^106→R^36。
F,B,C の定義と `rel=z*d*F*B*C−1`、全36 target componentへの対角追加列を
確認する。aは係数体の非零元であり、Jは反転しない。dは alternative ad57
row basis の比較条件として必要である。保存36×70のQQでのrow equivalenceと、
全prefix操作の32003整性を再検証したとはしない。今回の独立還元は保存された
polynomial presentation自体に対するものである。

保存sanityの2×3/3×2積を記号的に確認する。sanity logの成功はその小例に
限り、main logはINPUT_READYだけでSTD/PROFILE完了markerを持たない。
runnerは準備前からの60秒deadline、子2GB、全artifact32MiB、単一main実行と
process-group停止を実装している。記録はdeadlineによるreturncode−9で、
receipt出力を含むwalltime60.106秒。これは完了GB・OOM・数学的rank失敗の
証拠ではない。gzip復元SHA、全log/input/hash、artifact上限も照合する。

従って **unit/nonunit、dim、vdim、標準基底、annihilator、有限性の計算結果は
一つも得られていない**。QQへの移送、新しい成分、fullquartic全射も結論しない。

## 証明へ進むための有限次元縮約（条件付き）

二つの既証定理を混同しない。

* cubic36については、別の受理済み full cubic codimension theorem が Z3 の
  core次元≤1を与える。adFBC開での実際の57行基底と21単元pivotにより、この
  cokerの支持と比較できる。
* fixedM319の余次元≥2は **quartic coker** の定理である。それだけからcubic
  cokerの支持次元は従わない。また Z3 だけでは quartic bad curve cover の
  Z_any65 部分を尽くさない。

各定理に対応するcokerを QQ(a) へ移すと、支持は零次元以下となる。
一次元core成分のうちaが定数のものは消え、aを優越するcurveの一般fiberは
零次元だからである。従って有限生成cokerは QQ(a) 上有限次元である。
これは有限basisの**存在定理**であり、今回のmodp profileで有限性を観測した
という意味ではない。generic-aで消した有限個のa-special componentsも残る。

移送用の明確な小目標は、DVR

    V=Z_(32003)[a]_(32003),  Frac(V)=QQ(a),  residue(V)=F32003(a)

上での実coker M に対し、有限standard monomialsへの **exact spanning
identities** を得ることである。例えばe,d,zのmonic annihilatorと実module
membership証明、または単元leading coefficientsを持つ有限module reductions
で十分。これにより M は有限V-moduleとなり、その同じpresentationのspecial
fiberが零と証明できればNakayamaでM=0となる。free性やflat性は不要である。

別案としてe,dのmonic必要方程式から有限V-algebra A0を構成し、その全体で
dFBCが単元であることを示せばzもその有限代数内に回復できる。special fiberで
の単元性からの移送は、この **先行する有限性** があって初めて正しい。
inverse変数を追加しただけでは有限にならず、V[x]/(x−p)でxを反転すると
Frac(V)へ逃げる例がある。QQ上だけのmonic関係の係数がV整であるとも限らない。

現在はexact有限spanning/monic証明、全localizerのno-escape、完了special-fiber
検査のいずれもない。したがって同じsolverの再試行を提案・実行せず、この
有限証明の不足を次の検討対象として保存する。quarticには同じ論理を適用できるが、
元M319または完全なbadcurve coverに対する別の実certificateが必要である。

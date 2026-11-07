# 非被約R・m四squarefree leadingの全源は本質三八成分へ

2026-11-07、主義務G及びpointwise E認識。producer全source証明とexact full150 certificate。

**命題。** rank二plane S⊂SingD、actual D被約既約quintic、S沿い位数四、generic normal leading quartic squarefree、
πSの零cubic Rが非被約、という全rank四Poisson sourceは、既受理のessential-three八成分の閉和に入る。
特に全像は次元二十三以下でG非優越。
これはgeneric点だけの上界ではなく、指定した全sourceのpointwise recognitionである。

## 1. 全normal leadingからrank一又は六actual Dへ

`PLANE_M4_NONREDUCED_SQUAREFREE_QUINTICS.md`の全root/pencil/normalization証明を使う。
parallel tensor G4/fはJ²/I³をconstantにする。
I,J両非零ならΔ=cL6となり、**P2\L=A2全域**のfinite étale root coverから四rootがconstantとなる。
元のRが2L+Mで接続がM上に高位極を持っていても、このcoverはG4自身のdiscriminantだけで定まりM上もétale。
接続をMへ延長するという仮定を使わず、G4=L*g4、image-rank一が従う。

I0/J0の全stable linear spacesはpure Witt/Hankel/symmetric-pencil論証でrank二以下。
rank一は同じ枝へ送る。
rank二で三used linesを持つ型は非被約Rのsupportに入らず、二used linesの全三pencilsだけを保持する。
R=3Lではそれらも不可能。
R=2L+Mの全actual Dはnormal/base actual shears、diagonal scales及び有限swapの下で次の六つとなる。

| 全normal form | actual quintic F |
|---|---|
| J0 secant、残係数両零 | Au4+Bv4 |
| J0 secant、一方非零 | Au4+Bv4+u3v2 |
| J0 secant、両方非零 | Au4+Bv4+u3v2+u2v3 |
| J0 compression | Au3v+Buv3+u5+v5 |
| I0 mixed、残係数零 | Au4+Buv3+v5 |
| I0 mixed、残係数非零 | Au4+Buv3+u3v2+v5 |

normalization前の全nonzero amplitudes及び全零係数strataは含まれている。
compressionの両端係数非零等はactual Dの共通u/v因子排除から従い、chosen generic点で固定したものではない。

## 2. 全canonical150源のactual Pfaffianは恒等零

各Fの全canonical必要源

`LF={Π: 五変数(A,B,u,v,w)の全quadratic係数150, divΠ=0, Π#dF=0}`

を直接構成する。canonical strict liftとPfaffian–Casimir identityにより、
actualD class[F]を持つ全projective Poisson tensorはP(LF)に入る。
vertex w dependenceを先に捨てず、full1050 Casimir coefficientと25 divergence coefficientを全て用いる。

`CLASSIFICATION_PLANE_M4_NONREDUCED_SQUAREFREE_CASIMIR_20261107.json`は全六matrix、完全kernel、
非零maximal minor及び全basisを保存する。
六LFのdimensionsは順に5,3,2,1,1,1である。
それぞれ完全kernelの全linear combination Π(k)を作り、全五principal Pfaffiansをpolynomialとして計算した。
結果は**全六源で全五Pfaffianが恒等零**。
必要sourceの一部のbasis vectorだけでなく、全quadratic k_i k_j cross termsも含めた恒等式である。

従って全Π∈LFのraw rankは高々二であり、projective quotientでもrank四にはならない。
全actualD正常形の逆PGL変換はrankを保つため、rank二leading-pencilの全枝はrank四空。
full nonlinear Jacobiを未検証のまま満たすと仮定したのではなく、それを課す前の必要源全体で空性を得た。

## 3. rank一の全sourceを八成分へ戻す

残る全leadingはG4=L*g4で、D=L*g4+g5はL,u,vの三変数だけに依存する。
被約既約quinticなのでessential-variable数二以下はbinary factorizationに反し、exactly三となる。
`ESSENTIAL_THREE_QUINTIC_EIGHT_COMPONENTS.md`及び全source/curve/ambient別作者reviewにより、
その全Poisson sourceは既知八23成分の閉和に入る。
これで本命題のpointwise recognitionとG上界が従う。

## 4. 検証状態と射程

180秒/1536MiBのfinite job `jobs/plane_m4_nonreduced_squarefree_casimir.json`をgate PASS後に実行し、
18.472秒、最大125,776KiBで全exact QQ検査PASS。
gate自体を数学の検証とは呼ばない。
reportはfull150/1075、全kernel、universal全Π、全五Pfaffian vector及びactual Euler-Pfaffian零恒等式を含む。
以前の全二十九FunctionDef ASTと全immutable source/report/manifestを不変に保ち、実行producer bytesを埋め込んだ。
全root-cover/normalization数学と全六matrix/Pfaffianは、producer import無しの別作者reviewを別receiptに固定する。

generic repeated leadingの非被約Rについては、image-rank一/22/4は既八源、31は別の全-R theorem、
211高rankは`PLANE_M4_TWO_LINE_COLLISION_EXCLUSION.md`で全rank四空となる。
これらと本命題を合わせれば非被約R全m四がG非優越。
さらにreducedRの全m四既存packetを合わせた全R結論は、各独立receiptの受理後に採録する。
m二/三、全rank二plane及び全Gの完了は本命題から推論しない。

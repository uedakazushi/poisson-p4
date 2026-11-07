# nonredR・m4 squarefree全source認識の独立review

2026-11-07。producerとは別作者の数学的読解とexact Fraction算術。**PASS**。
実D被約既約quintic、rank二planeのzero cubic Rが非被約、plane沿い位数四、
generic normal quarticがsquarefreeである全sourceは既essential-three八成分で被覆される。
高rank special-pencil全源は空、残るrank一は全actualD essential-threeというpointwise認識である。
m2/m3や他のrepeated-leading packetをこの結論に混ぜない。

## 1. 非例外jのroot-coverを使うopenは元接続より広い

`PLANE_M4_NONREDUCED_SQUAREFREE_QUINTICS.md` と
`PLANE_M4_NONREDUCED_SQUAREFREE_RECOGNITION.md` の全数学を読解した。
Hamiltonian parallel tensorからinvariant rational function J²/I³はconstant。
I,J両非零ならUFDとbase degrees二/三からI=aL²、J=bL³、Δ=cL6、abc≠0である。

ここで使うcoverは **G4そのもののroot incidence** である。
Δのzero supportがLだけなので、P2\L=A2の全域でfinite étale degree四となる。
元R=2L+MのconnectionがMで高位極を持つ場合も、G4のcoverはM上でétale。
元connectionをMへ延長することやlog/simple-pole性は一切仮定しない。
このopenの区別が本拡張に必要な点であり、reducedRのpole argumentをそのまま流用していない。

A2(C)のsimply-connected性によりcoverは四copyへsplitし、Pic(A2)=0から四rootはpolynomial vectorsで表せる。
各pair determinantとoverall factorは nowhere-zero polynomial、従ってnonzero constant。
base polynomial degreeが高々一なのでroot-factor degreesの和も高々一。
少なくとも三rootがconstantであり、二つのindependent constant rootsとのdeterminantsがconstantなら残るrootもconstant。
よってG4=L*g4。この証明は全nonexceptional sourceを覆い、chosen connection matrixに限定されない。

I0/J0では既全Witt/Hankel/symmetric-pencil補題を、R multiplicityに依存しないpure algebraとして適用出来る。
rank一は同じessential-three枝へ。
rank二のI0三相異型は三used linesを必要とするため非被約cubic supportには入らない。
残る三stable pencilsは二相異used linesを要し、R=3Lにも入らない。
R=2L+Mで残る全三pencilsはsecant、compression、I0 mixedであり、両line orderingを保持する。

## 2. 全actualDの六正常形・全零非零補集合

全shearsはbaseA,Bへnormal linear formsを加える可逆GL5変換。
secantはu5,u4v,uv4,v5を除いた後、残るu3v2,u2v3係数の両零・一方非零・両非零を保持する。
一方非零の二orderingはnormal/base swapで同じ、両非零は二normal ratiosとfinite root choicesで一へ戻る。
compressionは両endcoefficients u5/v5が非零でなければ共通u/v因子を持つため、両方一へ戻せる。
I0 mixedはv5 coefficient非零が既約性から強制され、残るu3v2 coefficientは零/非零を別に保つ。

従って全actual quinticsは次の六つで尽くされる。

\[
\begin{gathered}
 Au^4+Bv^4,\quad Au^4+Bv^4+u^3v^2,\quad
 Au^4+Bv^4+u^3v^2+u^2v^3,\\
 Au^3v+Buv^3+u^5+v^5,\\
 Au^4+Buv^3+v^5,\quad Au^4+Buv^3+u^3v^2+v^5.
\end{gathered}
\]

各FはAに関してprimitive linearでC上絶対既約、四gradientはconstant independentなのでexactly essential-four。
全higher sourceも同じlinear inverseで戻せ、sourceの一部だけをnormalizationしたものではない。
secant両残係数零を既約性の名目で落としていないことも確認した。

## 3. full150と全五Pfaffianの独立再構成

canonical strict divergencezero liftのCasimir identityから、各固定actualD classの全tensorは
全150quadratic coefficientsを持つL={divΠ=0,Π#dF=0}に入る。
vertex wを含む全monomialsを保持し、任意Euler代表をcanonical liftと取り違えていない。

新modeはproducer/CASをimportせず、旧独立Fraction engineだけで、各Fの全1050Casimir行と25div行を再構成した。
全六1075×150 matrix、全basisとmonomial順がproducerと一致。
maximal nonzero determinantsとcomplete kernelの全original-row代入、free-coordinate identity minorを検証した。
独立kernel certificateによりdimensionsは順に5,3,2,1,1,1。

各complete kernelの全basis tensorと全pair sumsについて、全五principal4×4 Pfaffiansを直接計算した。
全diagonal/cross quadratic coefficientsは全六sourceで恒等零。
これにより任意のlinear combinationでも全五Pfaffiansが零、rawrank≤2、projective rankも≤2である。
secant両残係数零sourceにはnonvertex bracketsもあり、共通vertexを仮定するshortcutは使っていない。
全Pfaffian polarizationsを直接検証したので、この場合も全sourceの空性が従う。
full nonlinear Jacobi又はambient tangentを追加する必要は無い。

## 4. 残る全rank一sourceのpointwise認識

G4=L*g4ならactualD=L*g4+g5はL,u,vだけに依存する。
既約quinticは二変数以下にはなれないのでessential-variable数はexactly三。
既受理の `ESSENTIAL_THREE_QUINTIC_EIGHT_COMPONENTS.md` と全source/curve/ambient別作者receiptにより、
全Poisson sourceは八つの閉23成分に入る。
高rank六正常形が空であることと合わせ、指定したnonredR,m4squarefreeの全点認識とG≤23が従う。
既知点だけの帰属をfamily全体へ外挿する議論ではない。

## 5. 独立receiptとscope

finite契約 `jobs/plane_m4_nonreduced_squarefree_review.json` をgate後に一回実行し、
1.057852秒、RSS57,460KiB、全exact独立算術PASS。
旧reviewerの全七FunctionDef ASTは不変で、isolated modeだけを追記した。
旧reviewer sourceと契約exactbytesは前receipt/manifestから復元可能。

| artifact | SHA256 |
|---|---|
| code/classification_plane_quintic_casimir_review.py | ed194c995b2986370d989fa32a41b02f64a56f314c782ce63df7aee9ecee9548 |
| reports/CLASSIFICATION_PLANE_M4_NONREDUCED_SQUAREFREE_REVIEW_20261107.json | 86a652dfa88491811cec1f4cb9de0c336eca5bdd6fed749ab0e2e254268e6e14 |
| research/classification/jobs/plane_m4_nonreduced_squarefree_review.json | 03f8fd9980fa2e77a6eea4cc131648da96ef2f76d62351541001533153ce3300 |

gate自体は数学を認証しない。全入力hash/実行reviewer bytesはreportへ保存した。
nonredR,m4のrepeated31/211等はそれぞれのreceiptと合わせて初めて全体を論じる。
本receipt単独で全m4、m2/m3、全G/Eの完了を主張しない。

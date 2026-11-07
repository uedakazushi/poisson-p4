# rank二plane・被約可約cubic・m4非特別leading四次：独立レビュー

2026年11月7日。主義務G。
**transverse conic+line、tangent conic+line、triangleの全指定枝について、
射影Poisson像の次元上界24、23、23を受理する。**
指定枝とはplane沿い位数m=4、leading binaryquarticがsquarefree、classical I,Jがともに非零の場合である。
concurrentの上界31、I=0又はJ=0、repeated-root、m2/m3、非被約Rは未処理として保持する。
P1二次compatibilityを解決したという主張も含まない。

独立checkerは `code/classification_plane_reducible_cubic_m4_review.py`、
証明書は `reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SCALAR_REVIEW_20261107.json`、
有限契約は `jobs/plane_reducible_cubic_m4_review.json`。
producer又はproducer helperをimportせず、以前の別作者標準Fraction疎多項式・Schouten・有理行消去だけを使う。
CAS、SymPy、NumPyは不要。gate後0.637秒、最大RSS18,868KiBで全検査PASS。
gateは数学の証明と区別する。共有台帳・TeX・Gitは変更していない。

## 1. 全接続源とnonlog補集合を別に再構成した

Rの被約四型は `z(xy−z²), z(xz−y²), xyz, xy(x−y)`。
simple poleを許すtracefree接続はH⁰(Ω¹(3))⊗sl2の全24係数である。
Euler列からH⁰(Ω¹(3))は八次元であり、tracefree linear Vに対する
`i_E i_V vol` は八次元の単射、従って全源である。

checkerは八個のtracefree matrixを別に構成し、H,A,Bの全24係数を保って

    H(f)+det(E,A,B), A(f)+2det(E,H,A), B(f)−2det(E,H,B)

の全30三次係数を再構成した。
さらにdΩ+Ω∧Ωの三×三の全raw曲率係数を直接展開し、これらの式と一致することを確認した。
全sourceの30式を解いたとは主張しないが、非log sourceを初めから除外してはいない。

四型の各々について、公開Mから
`M⁻¹dM−(1/2)dlog(detM)I` の全numeratorを独立に作った。
全曲率零・horizontal・simple poleを確認し、あるcomponent g上で
`dg∧numerator` がgで割れないことを厳密多項式除算で確認した。
四型全てに真正nonlog flat witnessがある。
従ってirreducible cubicの場合のπ1による一律log化を、この結果へ流用していない。

## 2. parallel leading tensorから全接続が固定される矢印

plane上rank二、U=P²\RでπSが非退化とする。
normal-order≥2によりconormal作用はbase関数のliftに依存せず、Jacobiからflatである。
Hamiltonian flowはπ²とπSを保存する。
I-adic associated gradedへ制限すると `γ4=G4/f` がparallelとなる。
Hamiltonian場はTUを張るので、これは全方向についてのparallel性である。
log性、ordinary double、formal product splittingはこの議論に使わない。

相対GL2作用でI,Jのweightは二、三だからJ²/I³は不変でありconstant。
I,Jは非零homogeneous polynomialでbase次数二、三なので、
UFDでdivisorを比較すると `I=aL²,J=bL³`、Lは一次式となる。
squarefree性からdiscriminantは非零定数倍のL⁶。
parallel tensorのdiscriminantはU上で消えないため、LはRの一次componentである。

L以外のRの一般点ではG4の四rootが相異なる。
三rootを局所holomorphicに標準化するprojective frameを選べるので、
constant crossratioを持つ[G4]のprojective接続はその点で正則である。
四rootのprojective stabilizerのLie algebraは零であり、この接続は一意である。
γ4のscalar因子1/fはprojective接続に影響しない。
従って元のtracefree接続ΩはL以外のdivisorial poleを持たない。
simple-pole仮定とvector bundleのcodimension二の延長により
`Ω∈H⁰(Ω¹_P²(L))⊗sl2`。
Euler列のH⁰(Ω¹(1))=0からΩ=0。
これは指定leading-symbol枝の全nonlog接続も含む結論である。

projective接続零から[G4]はconstant、従って `G4=L*g4`。
scalar部分はparallel性から

    dlog(L/f)+2β=0,
    β−β_can=−dlog(L³/f)/6, β_can=dlogf/3

に固定される。右辺はEulerにhorizontalである。
m>2ではnormal quadratic Kは、非零πSとのPfaffianから零を強制される。
よってP0は一つの `P_f+Euler_normal∧V` に固定される。
Rの同型類とline Lの選択は四型それぞれ有限automorphismで尽くされ、隠れた連続parameterはない。

## 3. 全30/9行列と三角fiber上界

scalar一形式のnumerator `df/6−(f/L)dL/2` を、全八linear fieldのnumeratorへ展開してVを独立に再構成した。
base chartでP_fのprojective係数はEulerの公式により3fなので、この対応のfactor三も保持した。
公開scalar matrixと全係数が一致する。
P0自身の全350 Jacobiとactual Dが零であることを確認した。

全30のweight一basis、全九のweight二basisを別に作った。
各型について全350×30 Jacobi、D3=0の全24行、全350×9 weight二operatorを再構成し、
公開行列と全係数が一致することを確認した。
必要sourceの完全kernelも別の有理行消去で一致し、全rankには独立に計算した非零minorを付した。

| R | dimker(weight1 Jacobi,D3) | rank weight2 | fiber上界 | fixed-P0全源上界 |
|---|---:|---:|---:|---:|
| transverse conic+line | 8 | 6 | 3 | 11 |
| tangent conic+line | 8 | 6 | 3 | 11 |
| triangle | 8 | 6 | 3 | 11 |
| concurrent | 14 | 3 | 6 | 20 |

weight二の式は、任意の固定P1に対して同じD2のinhomogeneous linear systemである。
従ってその解は空又はkerD2のtranslate。
これは未知fiberのgeneric rankを全P1へ一般化する議論ではなく、同じ全operatorを使う全点上界である。
高weightの非線形Jacobiはこのsupersourceを狭めるだけである。

V2のG4=L*g4に対応する全十leading-shape行も別に再構成した。
追加後のP2 rankは前三型六、concurrent三で変わらず、全型で九未満である。
proper quadratic compatibility certificateは公開reportでもNoneのままであり、checkerもこれを保持した。
P1の二次compatibilityを解いた、又はconcurrentを閉じたとは扱わない。
受理する前三型の弱い上界には、この追加polarization試行の成功を仮定する必要がない。

## 4. 全PGL像への復元と補集合

baseのconnected Aut(R)次元dは順に1,2,2,3。
前二つはsmooth conic上の二点又は接点のstabilizer、triangleはdiagonal torus、
concurrentは三directionを固定するpencilのstabilizerから直接得られる。
planeとRを保つPGL5のconnected subgroupは、normalGL2の四とnormal→baseの六を合わせてd+10次元。
connected groupは離散的なline Lも保つ。
scalar接続は(R,L)の内在式であるから、f係数を一へprojectively規格化した全source Xはこの群で保たれる。

写像 `PGL5×X→P H⁰(Λ²T)` のfiberにこの群が入る。
PGL5因子への右作用が自由なので、個別tensorのstabilizerに関係なくfiber次元は少なくともd+10。
従って像の上界は `24+11−(d+10)=24,23,23`。
concurrentでは同じ弱い議論は31しか与えない。

以上により前三型の全m4・squarefree・I,J非零枝のG非優越を受理する。
concurrent、special j、重根、m2/m3、非被約cubic、P1二次compatibilityは明示補集合として保持する。
新成分、全E帰属、全plane分類の完了はこのレビューの結論ではない。

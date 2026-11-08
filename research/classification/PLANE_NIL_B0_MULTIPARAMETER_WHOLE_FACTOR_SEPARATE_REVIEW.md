# normalnil b零 mixed／Cartan 全源：別作者の whole necessary-kernel 監査

2026-11-08、PRIMARY E。producer
`PLANE_NORMAL_NIL_ZERO_B_MULTIPARAMETER_WHOLE_FACTOR.md`
(SHA256 `1c6be4d9853356fec10431811e80c85543b90d9a8a78f2a13263042d720ac5b0`)
を独立に読解し、全 source と全補集合を保ったまま **PASS** とする。
新成分は無い。実 D が被約既約 quintic である独立 T,V の対象枝は空である。
可約 D の全分類又は他 lane の完了は主張しない。

## 1. 対象入口と全 parameter

既独立受理 `PLANE_RANK_ZERO_M2_COMPLETION_SEPARATE_REVIEW.md` §§1–2 の
actual GL normalizers と全 35/30/9 source 表示を適用する。
非零 global V を持つ normalnil b零で T,V が独立な残二族は

```
mixed: T=H+tau*N,V=H+N,K=eta*z²,tau!=1,
Cartan: T=diag(-1-tau,2-tau,2tau-1)/3,
 V=diag(-1-eta,2-eta,2eta-1)/3,K=kappa*z²,tau!=eta.
```

H=diag(1,1,-2)、N=b∂a。normalizer は元の全 higher coefficients に作用する可逆写像である。
leading actual D は順に `-3b²c(tau-1)z²`、`-abc(eta-tau)z²`。
従って独立枝には higher multiplicity を minor の一般性で捨てた補集合はない。
tau=1 又は tau=eta は全 dependent source として
`PLANE_DEPENDENT_GLOBAL_VECTOR_POINTWISE_COVERAGE.md` へ戻る。
その新 H 二源は既別作者 `PLANE_DEPENDENT_H_WHOLE_EXCLUSION_SEPARATE_REVIEW.md`
で受理済みであり、旧 Ore constant-factor は現在の全 essential-three 八成分 theorem へ戻る。
G 上界だけを E membership に置換していない。

## 2. 今回新たに認証すべき必要 kernel

元 30 directions は basepairs の z/w 線形十八方向の後に D,B,E,C の順であり、
C は rows27..29 の `w²∂base∧∂z` 三方向である。
FULL9 weight二 source は basepairs×(z²,zw,w²) のみで、normal columns を持たない。
Π0 の z-column は z で割れる。
従って全 weight一 Jacobi 必要解で C=0 なら、実際の Poisson tensor の z-column は
full inhomogeneous weight二を戻しても z の倍である。
z=0 上では raw P も Euler も四次元 hyperplane 内にあるため、全 signed Euler-Pfaffian は z で割れる。
D=0 又は可約であり、対象の rank四・被約既約 quintic は存在しない。
この論証は J2 の injectivity、particular P2、higher compatibility の一般的可解性を仮定しない。

旧独立 m2 reviewer は必要な full operator、全 selected nonzero minors と全複素壁の被覆を検証したが、
G 上界に不要な全 kernel columns の exactness は一律には主張しなかった。
今回この差を埋めるのは、保存済み selected-row A_I の kernel だけで十分である。
r×r block A_(I,J) が非零で、K の J 補行が単位行列、A_I K=0 なら、
K は A_I の全 kernel を尽くす。full A の kernel はその部分集合なので C(K)=0 が結論を与える。
full A K=0 又は新 rank/operator を余分に要求しない。

## 3. 隠れた rational poles と minor 開集合

式として表示された K の分母を新しい除外条件にしてはいけない。
元の parameter chart の old nonzero minor 開集合では

```
K_J=-A_(I,J)^(-1) A_(I,free), K_free=identity
```

が adjugate による正則な一意解を与える。
保存 K との rational identity は free identity と A_I K=0 から従う。
従って表示 K に見える余分の pole はこの開集合上で可除である。
C 行の rational zero identity は、この正則延長でも全点で零となる。
old minor の本当の零集合は元の全 wall hierarchy に残してあり、generic K の外挿で置き換えない。
この議論は各 irreducible rational graph/curve と保存済み有限根にも個別に適用する。
非被約 source の nilpotent 構造を除外する必要はなく、E の全複素点結論として十分である。

## 4. 全階層の別 Fraction readback

新 isolated inline checker は producer/CAS を import しない。
既 m2 別作者 reviewer の Fraction sparse rational-polynomial field helpers のみを再利用し、
元四 report の保存 350×30 行列を読み、元 selected-row kernel を直接掛けた。
元 matrix はその保存 ancestor から辿り、各 inherited map を**順に、各段では simultaneous**に合成し、
最後に元 rational root を入れる。producer V1 の inherited line map omission は不変失敗記録として保持する。

finite job `jobs/plane_nil_b0_multiparameter_selected_kernel_separate_review.json` を gate PASS 後、一度実行。
`CLASSIFICATION_PLANE_NIL_B0_MULTIPARAMETER_SELECTED_KERNEL_SEPARATE_REVIEW_20261108.json`
は **1.094494202秒 / RSS205052KiB** で PASS。

| 元 report | 元 node 数 |
| --- | ---: |
| TWOPARAMETER / NIL_B0_MIX | 10 |
| THREEPARAMETER / NIL_B0_CARTAN | 226 |
| HYPERBOLA | 3 |
| CARTAN_SURFACE | 38 |

全277 path と hierarchy を独立収集して V2 receipt と一致確認した。
274必要 supersources で、元 I/J、nonzero minor 文言、自由列順、元 parameterization chain/root、
全自由行 identity、全3C zero を確認し、**40,521個の selected-row product identities** を Fraction で零と確認した。
新しい M1、operator、rank、minor、kernel、wall、点又は normal profile は生成していない。
旧全 nonzero minor と全 factor/root 被覆は既独立証明へ正確に依存し、再計算したとは主張しない。

## 5. domain と dependent 全補集合

mixed の残 hyperbola は `eta*tau-tau+3=0`、eta=1-3/tau。
tau=0 に点は無く、tau=1 は dependent へ戻る。
Cartan の二 surface は

```
F1=eta*tau-3eta+3kappa*tau+tau,
F2=eta*tau+2eta+3kappa*tau-3kappa-2tau-1.
```

actual base swap a↔b と tensor -1 は
`(tau,eta,kappa)->(1-tau,1-eta,-kappa)` を全 higher を含め可逆に実現する。
F1 の tau!=0 graph は kappa=-eta/3+eta/tau-1/3。
tau0 は eta0 を強制する dependent vertical kappa-line。
F2 の対応する tau1 は eta1 の dependent vertical line である。
従って graph 外を空と誤認せず、既 whole dependent theorem に戻す。

38件に含まれる三つの rank零 cleared matrices は
`all_whole_curve_walls/{0,6,7}/analysis/all_exceptional_roots/4` の tau0 だけである。
元 curve は順に eta0、eta=tau²、eta=tau²/(2tau-1) で、いずれも F1 graph の domain tau!=0 にある。
人工的 cleared tau0 を実 independent source の rank零 kernel として用いていない。
元 surface に戻せば上の dependent vertical line のみである。
元の他の全 line/curve/root、graph denominator complement は old factor/gcd 証明どおり保持する。

以上で独立 mixed/Cartan の全実 Poisson source で D の z 因子が従う。
dependent 全補集合の既 E 帰着と合わせ、この lane に未処理 parameter を残さない。
新 source・未知親の探索や nonlinear jet 延長は不要である。

## 6. 凍結と検証範囲

実行 inline source 全文、input hash、277 records/domain、274 exact identities の結果は新 report/job に埋めた。
producer の V1/V2、旧元四 report、旧別作者 receipt、全 runner/job は不変である。
新しい独立 note/report/job/manifest のみ作成した。台帳・原稿・Git は変更していない。

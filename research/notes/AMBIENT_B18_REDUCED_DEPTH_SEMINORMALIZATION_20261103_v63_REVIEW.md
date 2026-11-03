# generic B18 の被約交叉：exact depth と形式的半正規化

2026-11-03。独立した概念的検証。対象は D_red の generic B18 点での形式的完成環である。
元の非被約 D の深さ、元の D_red の conductor、収束的同値、global な半正規化と完成の
交換は主張しない。新しい CAS や物理ジェット計算は用いない。

## 被約交叉の深さ

係数体を k=C と書き、既受理の同時 G19 分解を取り除く。
N slice の全環を A、R slice を B=k[[T]] と書く。
全 N completion により `dim A=9, depth A=1`、normalizer は正則9次元 B_N であり、
`A⊂B_N` は有限 subintegral である。
独立した全150三次 graph 検証と h=4 により、実際の交叉商は C=k[[T]]/(T⁴) である。
被約 slice は、元の ambient から誘導される実際の全射に沿う fiber product

    F=A ×_C B

である。この式で A→C の全37座標の写像を別の抽象的な長さ4写像に置き換えていない。
その全写像を展開した full normal form も本稿では主張しない。

全射の差を使う完全列

    0 → F → A⊕B → C → 0

で、中央の深さは1、右の深さは0である。depth lemma の異なる深さの場合により
`depth F=1` が**正確に**従う。同時 G19 分解を戻すと D_red cone の深さ20、
smooth scalar を除くと projective の深さ19となる。
これは actual D の深さ19–20という別の二択を決定しない。

## 半正規化の候補と有限 subintegrality

slice の候補は

    S=B_N ×_k B,   C_red=k.

N の正規化と R の正規化はそれぞれ B_N、B なので、F の全正規化は
`M=B_N⊕B` である。F→S→M は自然な包含である。
M は F 上有限で、S はその F-submodule だから S も F 上有限である。

F の素点は、A と B の素点を交叉商 C の被約一点で貼り合わせたものである。
S では B_N と B を同じ一点で貼り合わせる。
A→B_N が subintegral なので、素点と各剰余体は変更されない。
従って F→S は有限 subintegral である。
19変数を戻した場合は、貼合せの一点を `P=k[[g19]]` に置き換える。
交叉商は `P[[T]]/(T⁴)`、その被約環は P であり、同じ素点・剰余体の議論が成り立つ。

## 半正規性の直接証明

S の全商環は `Frac(B_N)⊕Frac(B)` である。
z²,z³∈S を満たす z を全商環から取る。
各成分の二乗が normal な B_N,B に入るので、各成分はそれぞれその環上整であり、normality により
z 自身が M に入る。二成分の P への像を a,b とすると

    a²=b²,  a³=b³.

P は domain なので a=b である。従って z∈S。
square/cube criterion により S は seminormal である。
F→S の有限 subintegrality と合わせ、S が F の形式的半正規化そのものである。
非被約な D 自身を reduced ring として扱ったり、その nilradical を半正規環に残したりしていない。

## 明示環、conductor、数値

P に整合する正規化座標を選ぶと

    B_N=P[[x1,...,x9]], B=P[[T]],
    S=P[[x1,...,x9,T]]/(T*x1,...,T*x9).

実際、両成分に共通する係数は P、残りは x のみの項と T のみの項であり、混合項だけが零となる。
この ideal は `(T)∩(x1,...,x9)` であり、従って S は reduced で、二つの次元の異なる成分を保つ。

S→M の conductor は正確に

    (x1,...,x9)B_N ⊕ T B

である。必要性は M の idempotent (1,0),(0,1) を掛ければ従い、十分性は両剰余値が零になることから従う。
S 内では `(x1,...,x9,T)`、商は P である。従って conductor の台は reduced smooth G19 である。
**これは S→M の conductor** であり、F→S や F→M の conductor の計算ではない。

slice の埋込次元は10、次元は9、深さは1である。
深さは `0→S→B_N⊕B→k→0` からも、`T*(x9)` の Koszul 分解からも得られる。
graded Hilbert series は

    1/(1-z)^9 + z/(1-z) = [1+z(1-z)^8]/(1-z)^9.

従って top-dimensional multiplicity は1である。低い次元の R 成分を multiplicity から消したわけではなく、
Hilbert--Samuel 最高次係数に寄与しないことを使っている。
cone では `(dim,edim,depth,e)=(28,29,20,1)`、projective では `(27,28,19,1)` となる。
重複度1から regularity を推論してはいけない。この環は異次元成分を持ち、CM ではない。

コードと flat receipt はそれぞれ
`code/ambient_B18_reduced_depth_seminormalization_review_20261103_v63.py`、
`reports/AMBIENT_B18_REDUCED_DEPTH_SEMINORMALIZATION_20261103_v63_REVIEW.json`。
入力ゲートの hash と、monomial intersection、Hilbert 係数、Koszul Betti 列を固定する。

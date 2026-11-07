# 二本の節点接線を持つ三次曲線の全磁気源：独立確認

2026-11-07。課題 `E_reduced_nonnormal`。
`PLANE_CUBIC_TORIC_COLLISION.md` の全源被覆を、producer を import しない
別実装と数学的読解により確認した。これは他の三配置や全分類の完了ではない。

## 独立に再構成した算術

実装は `code/classification_plane_cubic_toric_collision_review.py`、契約は
`jobs/plane_cubic_toric_collision_review.json`、証明書は
`reports/CLASSIFICATION_PLANE_CUBIC_TORIC_COLLISION_REVIEW_20261107.json`。
有限契約の gate を通過後、QQ 上の恒等式を実行した。
計算には SymPy を使うが、producer の関数・出力括弧を生成に用いず、
十個の独立係数を持つ任意の三変数三次式から始めた。
実測0.979秒、最大RSS 53,512 KiB。

五変数の順序を `(H,A,B,C,W)` とし、H 行を零にした

    {A,B}=f_C, {A,C}=−f_B, {B,C}=f_A,
    {A,W}=−λWA, {B,W}=−λWB, {C,W}=−λWC

について、全十個の五変数 strict Jacobi 式が零、
`Pf(P)=−3λWf` を直接再展開した。
三変数 Jacobi bivector が二次斉次で Euler 場に保たれるという幾何説明とも一致する。

六係数二次式 `R(A,B,C)` の原始式は、producer の表示を転記せず
`∫3R dB` で作った。従って全六係数を独立に保持したまま `G_B=3R` となる。
`f=−λABC+g(A,C)+εG` を代入し、五次元斉次座標行列

    (H,A,B,C,W)=(h,x,y,h+εz,x+εw)

を実際に逆行列化した。この順序で行列式は `ε²`。
逆行列による tensor の変換と ε 倍から全係数を生成すると、
非定数分母は完全に消え、全括弧は五座標の斉次二次式となる。
producer の六括弧との一致、全十 strict Jacobi、
`Pf=−3λ(x+εw)f(x,y,h+εz)` を確認した。

`λ=−1, g=−(A−C)³/8` とし、
`Z=−8z/3, V=(z+w)/8`、tensor の24倍を別に実行した。
ε=0 で得る全5×5反対称行列は、既存の
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json` にある
`node_two_node_tangents` の全行列と完全に一致した。
とくに六個の磁気係数による `(Z,V)` 成分の微分は正確に

    24h², 24hx, 24hy, 24x², 24xy, 24y²

である。源の一般点だけ、又は磁気項を零にした点だけの一致ではない。
全源の Pfaffian は `−576x[hxy−(x−h)³/8]` で、全磁気係数によらず非零。

## 全源所属へ進む幾何的論理

任意の三次式 f と `λ≠0` の補助源は既約であり、その像閉包も既約である。
独立に `f=ABC, λ=−1` から再生成した点の射影退化式は
`3HABCW`。五因子は独立座標なので、これは一般階数四の SNC 点である。
補助像閉包を含む全 Poisson 空間の既約成分を一つ選べば、その成分はこの SNC 点を含む。
`KNOWN_COMPONENT_AUDIT.md` 第4節で採用している MPS Theorem1.4 により、
SNC log-symplectic locus と交わる成分は既知40成分に限られる。
従って補助源の全像は既知 MPS 有限閉和に入る。
この議論に、選んだ toric 点で成分が一意であるという追加仮定は不要である。

ε≠0 の衝突族は補助源の PGL 変換と非零 scalar 倍である。
連結な PGL5 は各既約成分を保ち、有限閉和も保つ。
最終的な tensor は ε の多項式なので、ε=0 の全点もこの有限閉和に入る。
全六磁気係数が独立に残る完全一致により、対象源の全点が被覆される。
対象の r 係数を一つも反転しておらず、その閉補集合の追加義務は生じない。
唯一の反転 ε は PGL の解釈のためであり、その零集合そのものが被覆対象として保持された。

ここで再証明したのは明示変形と全係数恒等式・上記の閉包論理である。
外部の MPS 定理自体と、既存十三配置への全支持還元は依存入力として区別する。

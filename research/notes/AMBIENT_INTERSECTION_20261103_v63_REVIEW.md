# S61 全被約交叉と実交叉非被約性：独立監査

2026年11月3日。rank≥2、全非零 rank1、実 target 接空間の三つの frozen producer を
別ディレクトリで QQ 再実行し、秒数以外の全証明書が一致した。
独立コード code/ambient_intersection_review_20261103_v63.py と
reports/AMBIENT_INTERSECTION_20261103_v63_REVIEW.json が追加監査を保存する。

受理する全域の結論は **S61 内の** (CA∩CB)_red=N である。
元の射影 Poisson 空間へ戻すと generic Poisson rank4 開部分の主張となる。
generic Poisson rank≤2 の境界だけにある交叉成分の有無や、射影成分全体の R1 は決めない。
一次ジェット Kβ の行列階数と Poisson の generic rank は異なる。
AMBIENT_INTERSECTION_SYMPLECTIC_SCOPE_20261103_v63.md を明示的な入力とした。

rank≥2 では非臨界全旗ファイバーは受理済みの被約一点である。
残る τ=Δ=0, rank Kβ=2 では H が非零 rank1であり、Levi SL2、shear、torus によって
B=C=1 の代表へ移る。shear と torus の全二形式係数変換を独立に検査した。
全旗 GB=(r,q,w²) と chart 外の排除から幾何学的旗は一つである。
従って二つの成分の持ち上げ旗は同じとなり、受理済みの全固定旗 affine10平面の交叉が N を与える。

rank1 の非零閉 primitive 形式は ell(y)(ell∧m) で {ell,m}=0 であり、
Sp4 の一軌道 x dx∧dz となる。真の追加収縮条件を使うと許容旗の全幾何学的台は
P(span(ey,ew)) である。全 λ chart と残る一点の CA/CB 四族について、独立に保存した
実 quadratic tensor の divergence、r=0 の共通旗因子、全61座標の β 一致・q 原点零を検査した。
この frame の全旗で perspective は零なので、同じ target q は同じ実 quadratic tensor を表す。
δQ≠0 なら像が両方の旗を回復し、同旗交叉の r=0 と矛盾する。δQ=0 なら一方の旗で r=0、
そこで CA=CB=N の固定旗族に入る。零 β 層は受理済みの (M0∩M1)_red=R0 が尽くす。
これで有限型 C スキームの全閉点から S61 全被約等式が従う。

実交叉の厚みについては、正規化の通常の微分だけを使っていない。
各 full pinching 環 Ac の極大 ideal は regular target K,s と p²,p³,I*p で生成される。
独立コードは保存した全61出力から、選択した **実 target K を全次数で固定する正確な逆行列** を
構成し、その全式の係数・mixed derivatives から全 cotangent を作り直した。
producer の source を三次まで再帰する手順とは独立である。選択 block の determinant は非零定数、
従って逆行列の分母による未知の除外点はない。

再構成した三つの接空間は保存行列と QQ(c0,c2,c5) 上で列空間が正確に一致し、
slice の次元 N19,CA21,CB21、union22 を得た。
全行列成分が三つの自由係数の多項式であるため、union の全23次小行列は恒等的に零であり、
全特殊化で union≤22 が成立する。個々の次元は full Ac 表示から全非正規点で一定である。
ゆえに実 scheme intersection の接空間は slice で少なくとも20、共通群因子7を戻して少なくとも27。
被約台 N の接空間26より大きく、交叉は指定 a=1 Q と共通非正規10次元層の全点で非被約である。

これは正規化の tensor product の余分な nilpotent を target に移した議論ではない。
実 ambient61座標内の二つの full target 接空間の交叉を使っている。
交叉 ideal、nilpotence order、primary decomposition の全計算はまだ行っていない。

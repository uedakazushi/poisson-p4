# cross-block edge 恒等式と NEW source B 次数除外の独立検証

2026-11-04。二つの凍結 packet をまとめて独立検証する。著者 producer を
実行せず、一般4×4行列 C と D=diag(x,−x,y,−y)、受理済み graph の c,d,f
から Newton の pencil 係数 p2,p3 を直接計算する。さらに定義の trace word
から N12=2(x²−y²)²(gj−hi) を全係数で計算する。従って edge 用の u,v を
source u,v と混同しない。

これらより K16 と G32 の cofactor を照合し、D0p4≠0 の既知 chart 上で

```
G32=64*p4*D0²*(g*j)*(h*i)
```

を確認する。抽象的な F4,p2,p3,p4,N12 の五変数では36項を全て検査する。
primitive multiplier は正確な整数として計算し、旧 certificate の末尾 .0
の文字列も有理数として読み、その正確な値との一致を検査する。浮動小数点
の近似検査にはしない。

次は受理済みの低次数 source 多項式の**B 次数と分母指数だけ**を使う。
K16 の B8 H16 clearing 後の各項の上界は2,1,2,2,2,1,0。G32 の B16 H32
clearing 後では D0*N12² の次数が正確に5、残り二項は高々4。
次数5の係数に現れる num(p0) と leadB(num(N12)) の全保存係数を照合し、
両方非零であることを確認する。QQ(a,Y,W) は整域なので積は非零であり、
積自体を展開しなくても他の項との非消去が従う。source32の展開・source
多項式同士の乗算・新しい点の評価はしない。

従って G32 は実際の NEW source 上で非零である。NEW の既約性、受理済み
octic base 優越性により D0p4 が一般に単元であること、既存点で確認済みの
bc*ef≠0による b=e=1 gauge の非空性を合わせる。有限 ordered frame の一般
開上で ghij≠0 が成立し、既知の Laurent chart は一般 NEW に交わる。

これは G32 を NEW 方程式として除外し、edge 境界への全体包含を除くもの。
残余八次式の L が NEW 上でも単元であること、NEW multisection の次数、
種数、新しい定義式や次数下界は結論しない。fullInv24の41/67の区別も変わらない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_edge_boundary_review_20261104_v63.py --output work/new-family-edge-independent
```

60秒・2GB上限。小さな抽象恒等式と保存済み係数 metadata の検査のみ。

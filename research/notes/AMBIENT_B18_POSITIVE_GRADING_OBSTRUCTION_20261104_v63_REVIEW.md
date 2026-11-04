# B18 の群由来 grading obstruction の独立検証

2026-11-04。原著 packet `AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63`
を、別の literal polynomial tensor と行列微分から検証した。原著 producer は実行しない。

実際の代表は座標順 `(h,x,y,z,w)` で
`P0=e_y∧(x²∂h+hx∂z+(h²+z²)∂w)` である。固定した四項から
`M P+P Mᵀ−(Mx)·∂P` を直接計算して全25 GL 方向を再構成した。
GL 軌道次元は19、`E_yy P0=P0` なのでテンソル scalar を加えても19である。
全 pure GL stabilizer は中心 `I` と
`E_yh,E_yx,E_yz,E_yw,E_wx` の張る6次元空間である。
5次元 tracefree 部分の任意の行列 `N` について、独立により強い恒等式 `N³=0`
を全係数で確認した。対角 stabilizer だけの検査ではない。

有限作用 `gλ=diag(1,1,λ⁻¹,1,1)` とテンソル倍率 `λ` が実際に `P0` を固定する。
受理済みの**全非被約 D**の物理接平面39と軌道19を用い、各物理 weight 空間への
射影の階数を独立に計算した。接平面は `0^32,1^7`、軌道は `0^15,1^4`、従って
横断20次元は `0^17,1^3` である。`u63=P^{xy}_{hz}` の微分は軌道上で零、
全接平面上で非零であり、weight0 の非零横断 cotangent class となる。
この線型函数自体の nilpotence はここでは主張しない。

## 全 cocharacter に及ぶ理由

中心 GL scalar は全二次 bivector 表現に自明に作用する。
pure stabilizer と SL5 の交叉の単位成分は、Lie algebra の全要素が nilpotent matrix
なので unipotent である。したがって pure GL stabilizer の effective torus rank は0。
extended stabilizer `S⊂GL5×Gm` を中心 GL で割ると、テンソル Gm への射の kernel は
連結 torus を持たず、effective torus rank は高々1である。表示した有限作用は非自明で
テンソル射への section でもあるため、rank は正確に1。

任意の cocharacter の像は `S⁰` に入る。連結線型代数群の maximal tori は内部で共役で
あり、上の rank1 torus は effective maximal torus である。`S` は `D` をスキームとして
保ち、軌道接空間も保つので、共役は横断接表現の weight multiset を変えない。
従って任意の非自明 effective cocharacter に17個の weight0 が残る。
effective に自明なら全 weight は0である。これで**群作用に由来し、P0を固定する
どの cocharacter も、全横断最小座標に strictly positive weight を与えられない**。
cotangent の双対作用または pullback の符号規約を変えても weight0 の結論は同じである。

任意の抽象的 formal local-ring action、零 weight parameter を残す相対的方法、
元 D の全 ideal・nil・socle・正確な depth の判定は対象外である。
この packet は新しい source jet、点、GB、原著計算の再実行を用いない。

再現コマンド:

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_positive_grading_obstruction_review_20261104_v63.py --output work/v62-review/B18-positive-grading-independent
```

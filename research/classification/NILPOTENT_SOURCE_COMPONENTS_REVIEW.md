# 冪零Jordan4から得た四成分の独立検証

2026-11-07。主課題 E_multiple_hyperplane。
`MULTIPLE_HYPERPLANE_NILPOTENT_SUPPORT.md` 第6節の四つの全源について、
全周囲の成分性と相互区別を受理する。
下記の再計算は全 source の有限還元を証明した元 producer の算術を import しない。
源の全射性は既存の別担当による全33未知数の独立レビューを入力とする。

## 今回独立に再構成した数値証明

`code/classification_nilpotent_source_components_review.py` は、
別作者による `classification_rank3_mixed_component_review.py` の疎な整数係数ループを再利用する。
元 producer の legacy Jacobi・軌道・階数 helper は一切 import しない。
単項式指数の列挙、全十 affine Jacobi と Euler 外積、全350行、全150列を再生成した。
源は元 JSON の完全な literal 係数を読み、`h²f(x/h)` で別に斉次化し、
全四パラメータを微分してから指定点へ代入した。

| s | rank D | vertical | vertical・軌道・scale | 全四源を加える |
|---|---:|---:|---:|---:|
| 1/3 | 103 | 24 | 47 | 47 |
| 1/7 | 103 | 24 | 46 | 47 |
| 1/5 | 103 | 24 | 47 | 47 |
| 1/15 | 103 | 24 | 47 | 47 |

全て QQ の完全な行還元であり、有限体標本による上界ではない。
各 rank の非零 minor、元の全行列 hash、実四源方向の全係数を保存した。
各代表の全350式と、全55本の実方向に対する D の零化を確認した。
さらに各 **全パラメータ源** の十 raw Jacobi を多項式として再計算し、全て零と確認した。

永続証明書は `reports/CLASSIFICATION_NILPOTENT_SOURCE_COMPONENTS_REVIEW_20261107.json`。
E_multiple_hyperplane の90秒・1536MiB契約で、実測6.779秒。
入力とコードの hash は証明書と `jobs/nilpotent_source_components_review.json` に保存した。

## 成分性と全源への論理

生の150係数から24次元の Euler 鉛直方向と射影scale一方向を除けば、
全接空間の次元は `150−103−24−1=22`。
全 PGL と実四源が与える次元下界も `47−24−1=22`。
従って指定点の局所環は22次元の正則局所環となり、その点を通る全成分は唯一である。
各源の rank四開は既約で、PGL 飽和像の閉包も既約である。
この閉包は上記点を含むため唯一の成分に入り、同じ次元を持つので成分そのものである。
これで全源の帰属が従い、一点の認識から根拠なく全 fiber を推論してはいない。

s=1/7 の実軌道は指定点で21次元であり、源を加えた22次元下界が必要である。
この一点の軌道次元を一般軌道次元と同一視しない。
他の三代表の軌道は22次元なので、それぞれの成分の稠密開軌道を与える。

## 相違の独立読解

全四源の実退化因子は非被約である。
非被約五次式は射影的積写像の像の有限和という閉条件なので、
全 PGL 飽和閉包にもこの条件が成立する。
従って SNC 点を持つ MPS 成分とは一致しない。
同次元の既存非MPS三成分とは一般因子型で異なり、他の次元の成分とは次元で異なる。

最初の二源は一般因子型 `2H+C3`、`2H+H'+Q3`。
第一の C3 は x1 について一次で、係数は非零倍の x4²、定数項は x4 で割れない。
従って primitive な一次式として Gauss の補題で既約である。
第二の Q3 は rank三なので既約である。
この二源は相互にも、最後の二つの二超平面型にも異なる。

最後の二源では multiplicity三の超平面が一意である。
正準 lift `Π=P−E∧N/5` とその超平面の式 ℓ に対して
`Xℓ=Π♯dℓ/ℓ` は射影同型で同変である。
全係数から直接

\[
 \Pi^\sharp dh=-hN/5,\qquad
 \Pi^\sharp dx_4=((25s+1)/20)x_4N
\]

と確認できる。s=1/5 では ℓ=h は `im N*` に入らず、
s=1/15 では ℓ=x4 は `im (N*)³` に入る。
従って二代表は射影同型ではない。
両者はそれぞれ22次元の開軌道なので、同じ既約成分に属することもない。
この区別は、単に同じ重複度の超平面を任意に名前で区別するものではない。

以上から四つの全源は、既存60成分とは別の四つの22次元全成分を与える。
他の modular 型の被覆や全分類は、この検証だけからは従わない。

## 再現

```sh
python3 scripts/classification.py gate research/classification/jobs/nilpotent_source_components_review.json
.venv/bin/python code/classification_nilpotent_source_components_review.py --output /tmp/nilpotent-components-review.json
```

旧成分の全証明書は今回再実行していない。
既存の `KNOWN_COMPONENT_AUDIT.md` と新rank三二重二次成分の証明を、相違の入力に用いる。

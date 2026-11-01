# rank-one A1射影端点の独立監査

2026年11月1日。
[主証明](S61_MINIMAL_RANKONE_ENDPOINT_20261101.md)の指定ambient開のt²=0と、
十六次元端点層の十分一般β上の除去を承認する。
全環自由66、他の旗・binary層、アフィン尽くしは承認対象でない。

## 独立QQ検算

元134式と普遍G6の保存係数から、主コードの局所生成関数を使わず次を再構成する。

- 全九source変数とsを保持した二重直線のε次数0,1、計270恒等式。
- 全分母、s=0の正則性、二接方向の独立性と良い素数32003での単元性。
- 元135式の指定点Jacobian、接核、左核、全二次障害。階数34/35/5/7。
- 同じ固定βの滑らかなbinary三次元族と、二重直線の実際の接平面。
- 元source β、有限endpoint行列のsymplectic性、全βのpullback、最高双ベクトルの変換。
- H9の有理復元、Sp4のsource法線rank2、βおよびjoint rank16。

q0の固定chart値15を維持し、1へ取り替えていない。
独立QQ検算は31.70秒でPASSした。
六次有限体閉包とMacaulay2のprimary分解は主producerを別に再実行し、29.75秒でPASSした。
この部分の独立アルゴリズム再実装は主張しない。

## 有限曲線から端点への被覆

全元A1曲線には `h(d−3ha)≠0` の非空開がある。
その一点で有限Sp4正常化を固定し、同じ固定β曲線を標準Q(0,0,1,1)を通るA1曲線として表す。
旧有限A1証明のbad-image除外と正常形を使えば、十分一般βの曲線にはこの表示を適用できる。
変換は曲線全体に固定して作用し、その同じ曲線の射影端点を移す。

正常化された曲線の最高w微分を直接計算し、主noteの有限endpoint行列が
最高双ベクトルを `−3Q0/h0` へ移すことを確認した。
全source βのpullbackとH9の開の復元も検算した。
従ってこれは一つの任意の部分族だけから端点全体を推論する議論ではない。

## 二次拡大と二重直線

canonical正常化のλは `λ²=1/dh` を満たす。
主producerは分子分母が同じ偶奇を持つLaurent多項式であることを確認し、
共通λの奇偶因子を消してλ²を置き換える。
この操作は有理式への正しい降下である。
その構成はsourceを読み主producerで再現したが、独立の二次拡大実装とはしない。
降下後の最終有理式が全270恒等式を満たすことは独立に確認した。

接minor非零により、完成局所環から `D_A=A[[s,ε]]/(ε²)` への写像は全射である。
この商はβ底A上flatであり、最大ideal filtrationも商filtrationに一致する。
八つの接平面一次式とv6²が実在商の随伴次数idealを与える。
被約方向のv6成分は零、nilpotent方向は593/972である。

## 局所環と相対降下

主計算の `J=B65∩Dlin` と同じbinary族の普遍下界65を用いる。
strict filtered quotientにより追加初期核は `Dlin/J⊂Fp[v]/B65` へ入り、
底上torsion-freeの階数65を下界が使い切るため、追加核は零となる。
二重直線商の核N0は底上torsion-free階数65である。

D_Aのflat性からTor完全列により `N/pN=N0` を得る。
良い還元・同じ特性零binary族・等階数特殊化を組み合わせる既監査のDVR補題が適用できる。
一般因子のHF=(1,7,23,34)、I⁴=0、t∈I²とD上のt=0からt²∈Nとなり、
Nのtorsion-free性で全体のt²=0を得る。
商が非被約であることは、このflat性と核の議論を妨げない。

p反転後の底完成と元q,t原点因子の完成を区別した忠実平坦降下も確認した。
同じflat商を相対H9へ保ち、全Sp4-sourceのβ rank16によるsmooth降下を行う。
既約十六次元Zの非空ambient除去開の補集合は次元15以下であり、
十分一般β上ではZの全点が除去される。

## 再現

```
.venv/bin/python code/s61_minimal_rankone_endpoint_review_20261101.py --output /tmp/rankone-endpoint-review.json
.venv/bin/python code/s61_minimal_rankone_endpoint_20261101.py --output /tmp/rankone-endpoint-replay.json
```

依存する旧有限A1・良い還元・DVR証明と、対象ソースのSHA256を証明書に固定する。

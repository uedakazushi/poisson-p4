# 新しいfirst-jet hypersurfaceの全20次不変式を除外：次数下界22

2026-11-03。既存の18次packetを変更せず、一回の新しいmodular probeと一回のexact QQ minorにより、**全20次Sp4不変式の制限が単射**であることを証明する。新成分の既約first-jet hypersurfaceの方程式次数は22以上となる。方程式そのものや22次の存在は主張しない。

## 完全な20次空間と一回のprobe

全16 weightを用いた整数DPとC2 Weyl交代和を20次まで計算する。各次数の総weight数を `binomial(15+d,d)` と照合し、

\[
\dim\operatorname{Inv}_{20}=36
\]

を得る。全8 Weyl項と全次数の値をprobe reportに保持する。

受理済みのω・moment・θ/H共役性により、長さ10の二元necklace108個のtraceは普遍的な20次不変式である。一回のprobeでは候補不足も検出するため、次も加えた。

* I4と、受理済みの全16次基底20個の積。
* 全8次基底5個と全12次基底9個の積45個。

合計173個を素数32003で、周囲の64点と実dense NEW25 sourceの64点に評価した。実点は既存20個のQQ点を保ち、44個の新しい整数source点を追加する。全source DAGのinverse node、value pivot、Pfaffianが非零となる点だけを使い、既存20点のβ16はQQ記録の還元とも一致する。

結果は108traceだけで周囲・実像ともrank36、全173個でも同じrank36である。したがってambient spanを不足した候補集合と混同していない。probeは60秒・2.2GB制限下の一回で完了し、実際の所要時間は約1秒であった。標本上のkernelを多項式方程式へ読み替えることはしていない。

## 同じ36式・36実点のexact QQ minor

選択36点を**実際のdense S61 DAGでQQを使って再評価**し、primitive β24と全16基底座標を再構成する。value pivotとPfaffianは全てQQ上非零である。

全10表現行列、ωとその逆、trace Gram逆から、FLINT有理行列コードがθ,Hを新たに作る。全108traceを全36点で厳密に評価し、全3,888個のQQ値が最初のmodular probeへ還元されることを検査する。選択した**36×36 QQ行列式は非零**であり、完全な有理数を `exact-evaluation.json` に保存する。このexact native計算は一回である。

この36式は実像への制限で独立であり、同時に周囲でも独立である。完全characterの上界36から、これらは20次不変式全体の基底であり、その制限が単射である。従って実first-jet像上で消える非零20次不変式は存在しない。

## 方程式次数と範囲

既存の証明により、実first-jet像はSp4安定な既約15次元錐hypersurfaceであり、その斉次生成方程式はSp4不変である。中心−Iから奇数次数は全て零。既受理の18次以下の除外と今回の20次により、生成方程式の次数は**22以上**となる。

次数22以上のprobe、巨大なLie制約GB、方程式の同定、全一般fiberの尽くしは実行・主張していない。

## 保存と再現

* `code/q_new_family_degree20_invariant_probe_20261103_v63.py`
* `code/q_new_family_degree20_invariant_exclusion_20261103_v63.py`
* `code/q_new_family_degree20_invariant_exact_20261103_v63.cpp`
* `reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_PROBE_20261103_v63/report.json`
* `reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63/` の全QQ値・exact minor・実QQ点・native入力とログ。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_degree20_invariant_probe_20261103_v63.py \
  --output work/degree20-probe-replay
.venv/bin/python code/q_new_family_degree20_invariant_exclusion_20261103_v63.py \
  --probe work/degree20-probe-replay --output work/degree20-exact-replay
```

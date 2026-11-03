# 新しい first-jet 超曲面の方程式次数は20以上

2026-11-03。新成分 `R4_POSITIVE_WEIGHT6_25` の実 first-jet 像は、既受理の通り W16内の既約 Sp4不変錐超曲面である。本 packet は**全18次不変式の制限が単射**であることを証明し、定義方程式の次数の下界を18から20へ改善する。方程式そのものは構成していない。

## 有限な候補空間と一回の probe

既受理の交代形式ω、moment map、二次行列θ,Hは全 Sp4作用に対して同じ共役変換を受ける。従って長さ9の二元 necklace60個の trace は普遍的な18次不変式である。

既存の完全 Weyl character 計算による `dim Inv18=9` を使う。なお `Inv14` は零ではなく**3次元**である。欠けた候補を見落とさないため、一回の modular probe では既存 Inv14基底3個とI4の積も追加し、計63個を評価した。

素数32003、24個の周囲の β と、既に保存された20個の実 QQ source 点で新たに評価する。実点は全 dense S61 DAGで再評価し、β16と value pivot/Pfaffianを既存のQQ記録の還元と照合する。60個のtraceだけで周囲・実像とも階数9であり、追加3積を含めても同じ階数9となった。probeは60秒・2.2GB制限内に一回だけ実行した。

この段階では sampled kernel を方程式と呼ばない。今回非零の9 minorが得られたので、次にその同じ9式・9点をQQで厳密に評価した。

## 元の実 source からの QQ 証明

選択9点を実 dense S61 DAG上で**QQを使って再評価**し、primitive β24とβ16を再構成する。各点の value pivot、Pfaffianは非零であり、定義域内にある。選択した60列中の9列は、necklaceの0始まり番号

```
6,10,12,18,19,21,23,26,32
```

である。0をθ、1をHとした9語は証明書 `samples.json` に保存する。

別の FLINT 有理行列 evaluator が、全60 traceを全9点で正確に評価する。単に有限体値をliftするのではなく、ω、trace Gram逆行列、全10表現行列とβ16からmoment、θ,Hを有理数で再計算する。全540個のQQ値が最初のprobeの還元に一致する。選択した9×9の**QQ行列式が厳密に非零**であり、その完全な有理数を `exact-evaluation.json` に保存する。

従ってこの9式は実像への制限で線形独立、特に周囲の全不変式空間でも独立である。完全な character による上界9と合わせ、これらは Inv18全体の基底をなし、その制限は単射である。実 first-jet 像上で零になる非零の18次不変式は存在しない。

## 定義方程式への帰結

既約な Sp4安定 hypersurface の高さ1素idealの生成元は、Sp4の指標が自明なので不変式となる。scalar安定性から斉次に取れ、中心−Iの作用により奇数次数はない。既存の完全な16次以下の除外に今回の18次を加え、定義方程式は**次数20以上**となる。

これは方程式の同定、次数20の存在・除外、全一般fiberの尽くしの証明ではない。高次数の探索は自動的に続けていない。

## 保存と再現

* Probe: `code/q_new_family_degree18_invariant_probe_20261103_v63.py`、`reports/Q_NEW_FAMILY_DEGREE18_INVARIANT_PROBE_20261103_v63/report.json`。
* QQ計算: `code/q_new_family_degree18_invariant_exclusion_20261103_v63.py`、`code/q_new_family_degree18_invariant_exact_20261103_v63.cpp`。
* 証拠: `reports/Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63/` の report、全540値とexact minor、実QQ点、native入力とログ。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_degree18_invariant_probe_20261103_v63.py \
  --output work/degree18-probe-replay
.venv/bin/python code/q_new_family_degree18_invariant_exclusion_20261103_v63.py \
  --probe work/degree18-probe-replay --output work/degree18-exact-replay
```

既存の表現・ω・θ/H共役性・Weyl上界と、今回の新しい実点・exact minorを区別して hash で結び付ける。

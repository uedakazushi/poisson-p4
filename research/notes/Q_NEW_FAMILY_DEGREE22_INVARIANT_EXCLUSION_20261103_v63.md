# 新 first-jet hypersurface の全22次不変式を除外：次数下界24

2026-11-03。新たな点は生成せず、既に受理済みの36個の NEW25 QQ 点だけを使う。一回の32003 probe と、それが full rank の場合にだけ実行する一回の exact QQ minor により、**全22次不変式の制限は単射**と証明した。NEW25 の定義方程式の次数は24以上となる。24次の probe、方程式の存在・同定は含まない。

## 21次元の空間と24個の普遍的な積

16ウェイトによる整数DPと C₂ Weyl 交代和を22次まで計算し、

\[
 \dim\operatorname{Inv}_{22}=21
\]

を確認した。各次数の全ウェイト数を binomial(15+d,d) と比較し、全8 Weyl 項を保存する。

候補は既に受理済みの基底だけから作る次の24個である。

* F₄ と Inv₁₈ の9基底の積9個。
* Inv₈ の5基底と Inv₁₄ の3基底の積15個。

これらは全て普遍的な Sp₄ 不変式であり、NEW25 の点に依存して式を選んではいない。

規約を固定する。F₄=μᵀGram⁻¹μ。Inv₈ は θ/H の長さ4 necklace の最初の5 trace、Inv₁₈ は受理済み長さ9 trace の列番号

\[
 [6,10,12,18,19,21,23,26,32]
\]

を使う。Inv₁₄ は長さ7 trace と取り違えず、受理済みの constant Poisson bracket 基底

\[
 \{I_{8,1},I_{8,2}\},\quad
 \{I_{8,1},I_{8,3}\},\quad
 \{I_{8,2},I_{8,3}\},\qquad
 \{f,g\}=\nabla f^{\mathsf T}\Omega^{-1}\nabla g
\]

である。基底番号は0始まりで、従来の全10 bracket 列では [4,5,7] に対応する。全候補の word と bracket pair を probe に保持した。

## 一回の modular probe

`Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63/samples.json` に保存済みの**同じ36個の QQ 点と β₁₆**を使う。点探索、周囲の追加点、別素数は使わない。各係数は32003で還元可能で、既存の value pivot と Pfaffian もその素数で非零である。

Ω、10個の元表現行列、trace Gram から θ,H を作り、Inv₈、Inv₁₈ と Inv₁₄ の勾配を求めた。modular 側の勾配は行列の adjoint contraction で計算する。この式を同じ最初の既存点で従来の全16方向の直接 dθ,dH と比較し、必要な3本の8次勾配の全48係数が一致することを確認した。

得られた36×24積評価行列の rank は **21**。character 上界21と一致する。候補の span が不十分かもしれないという条件を、この full rank によって解消している。probe は一回、約0.107秒で完了した。

## 一回の exact QQ 21×21 minor

probe が選んだ21行と21列を固定し、既存 QQ 点の記録をそのまま使う。新しい source 点も β 点も生成しない。

FLINT 有理行列コードは Ω と全10表現行列から θ,H を新たに作る。Inv₁₄ については modular 側の contraction を繰り返さず、元の

\[
 d\theta=\sum_j(\mathrm{Gram}^{-1}V^{\mathsf T})_{j,k}\rho(A_j),
\]
\[
 dH=d\theta+\Omega^{-1}(dV_kU^{\mathsf T}+UdV_k^{\mathsf T}),
 \qquad U=V\mathrm{Gram}^{-1}
\]

と trace の cyclic product rule で**全16方向を直接微分**する。

全21点で24個の積504値と、F₄・Inv₈五つ・Inv₁₄三つ・Inv₁₈九つの全378因子値を QQ 上で算出した。これら全てが probe の値に還元されることを確認した。選ばれた21×21行列式は非零である。

* 分子4922 bit、分母3186 bit。
* mod 32003 の値は3746。
* 完全な QQ 行列式と全QQ値を `exact-evaluation.json` に保持。

native exact 計算は一回だけで、compile・整合確認を含む所要時間は約1.72秒であった。各計算の上限は60秒・2.2GBである。

## 方程式次数への結論と範囲

21個の普遍的な22次不変式の制限が QQ 上独立で、全22次空間の次元は21だから、制限写像全体が単射である。従って NEW25 の first-jet hypersurface 上で消える非零22次不変式は存在しない。

既存の受理結果により、この像は Sp₄ 安定な既約斉次 hypersurface であり、その生成方程式は Sp₄ 不変である。18次以下と20次の除外を引き継ぎ、中心−Iが W に−1で作用するので23次を含む奇数次は存在しない。以上から

\[
 \boxed{\deg(\text{NEW25 の first-jet hypersurface 方程式})\ge24.}
\]

24次の式の存在、既存の24次判別式との同一性、全不変式環の表示、全 fiber の尽くしは結論しない。24次以上の評価行列へ自動的に進まない。

## 保存と再現

* `code/q_new_family_degree22_invariant_probe_20261103_v63.py`
* `code/q_new_family_degree22_invariant_exclusion_20261103_v63.py`
* `code/q_new_family_degree22_invariant_exact_20261103_v63.cpp`
* `reports/Q_NEW_FAMILY_DEGREE22_INVARIANT_PROBE_20261103_v63/report.json`
* `reports/Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63/` の全 QQ 値、21 minor、選択した既存点、native 入力とログ。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_degree22_invariant_probe_20261103_v63.py \
  --output work/degree22-probe-replay
.venv/bin/python code/q_new_family_degree22_invariant_exclusion_20261103_v63.py \
  --probe work/degree22-probe-replay --output work/degree22-exact-replay
```

全使用コード、基底と既存点 packet の SHA256 を保持する。今回の exact 計算が actual source DAG を再評価したとは記さず、既に actual DAG で受理済みの QQ 点記録を再利用したことを明記する。

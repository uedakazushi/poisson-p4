# 新成分の first-jet 超曲面の方程式次数は18以上

新成分 (R4\_POSITIVE\_WEIGHT6\_25) の primitive first-jet 像は、既に証明した通り (W_{16}) 内の既約な15次元錐である。その既約な定義多項式の次数は、既存の16以上から **18以上** に改善される。方程式そのものや18次での方程式の存在は主張しない。

## 完全な16次不変式空間への制限

既存の独立検証済み symplectic 表現、形式 \(\omega\)、moment map から構成した \(\theta,H\) を使う。両行列が同じ共役作用を受けることは、`Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW` において、全10生成子と anti-action の符号を含めて検証済みである。従って、行列積の trace、および既存の不変式の積は、全ての (W_{16}) 上で普遍的に不変である。

各行列entryは \(b\in W_{16}\) の2次式である。今回は次の60個の16次不変式を評価した。

- \(\theta,H\) の長さ8の二元 necklace 36個の trace。
- 既存の完全な12次基底9個と、4次式 \(\mu^T G^{-1}\mu\) との積。
- 既存の8次基底5個の対称積15個。

16個全ての重みを用いて \(\prod_\lambda(1-tx^\lambda)^{-1}\) の係数を整数で計算し、Weyl の交代和をとると、15、16、17、18次の不変式次元はそれぞれ0、20、0、9である。全重み多重度の和が \(\binom{15+d}{d}\) に一致することも各次数で確認した。

最初に素数32003で40個の実際の新成分の first-jet を評価し、60列の制限行列の rank が20となる20行20列を選んだ。有限体の使用はこの選択だけである。実際の18入力 S61 mapを **QQ上で** 再評価し、20点の primitive first-jet 16座標を生成した。各点の value pivot と Pfaffian は非零であり、定義域内にある。

その後、別の FLINT 有理行列実装で全60個を20点において正確に評価した。選択した20×20行列式はQQ上非零である。全1,200個のQQ値が、最初の有限体評価の剰余とも一致する。証拠にはQQ値全体、正確な行列式、実際の18入力点、primitive jet、native入力とログを保存した。これは補間による消失判定ではなく、非零行列式による線形独立性の証明である。

従って、これらのうち選択した20式は完全な16次不変式空間の基底をなし、新しいfirst-jet像への制限は単射である。この像を消す非零な16次不変式は存在しない。

## 定義方程式への帰結と限界

像は既約な (Sp_4)-不変錐超曲面である。その高さ1素イデアルの既約斉次生成元は、連結半単純群 (Sp_4) の指標が自明であるため、不変式である。既存の検証は14次以下を全て排除しており、15次と17次の不変式空間は零、今回16次への制限が単射なので、最小の可能な次数は18となる。

18次不変式空間は9次元だが、その制限や核は今回計算していない。一般first-jetの五点尽くしや、他成分の分類にも結論を加えない。探索を高次数へ自動的に延長せず、この一回の16次計算で区切る。

## 再現

既存環境を activate した後、別の出力先で実行する。

```sh
. /workspace/cloud-setup/activate.sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/q_new_family_degree16_invariant_probe_20261103_v63.py \
  --output work/degree16-probe-replay
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/q_new_family_degree16_invariant_exclusion_20261103_v63.py \
  --probe work/degree16-probe-replay --output work/degree16-exact-replay
```

公開 receipt は `reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63/report.json` である。独立検証は別の receipt として保存する。

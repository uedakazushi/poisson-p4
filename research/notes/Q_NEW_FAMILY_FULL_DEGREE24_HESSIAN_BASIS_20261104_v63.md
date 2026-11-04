# 全Inv₂₄の67次元基底：同一形式平面でのθ/H修正

2026-11-04。**普遍的なQQ不変式空間Inv₂₄の67本の基底を構成した。** C/D-only検査のrank8停止は保存し、旧16×16 θ/H modelの成分で補った。NEW25に制限した関係式の計算はまだ行っておらず、その方程式の既知次数下界24は変わらない。

## 二つの異なる表現

C,DはW₁₆から得た二つのquadratic adjoint covariantの4×4表現である。一方、旧不変式基底の

\[
 \theta=\rho_{W_{16}}(\operatorname{Gram}^{-1}\mu),\qquad
 H_{16}=\theta+\Omega^{-1}V\operatorname{Gram}^{-1}V^T
\]

は16×16である。Vのj列はΩρ(Aⱼ)β、μⱼ=βᵀΩρ(Aⱼ)β/2。旧producerの式をそのまま用い、H₁₆を4×4Dで置換しない。この区別がC/D-onlyの候補不足を補う。

旧probe最初の三つのambient vector β₀,β₁,β₂による**同じ**形式平面β=β₀+sβ₁+tβ₂、同じ32003上で一度だけ計算した。actual NEW25の点・別平面・source展開は使っていない。前packetの全C/D係数配列をhash付きで読み戻し、その上へθ/Hのwordを加えた。

## 全てのrank gateと選択した生成元

Inv₈の順序を(F₄²,p₁,p₂,p₃,p₄)とする。12次の5本F₄Inv₈と、次の4本を基底とした。

\[
 \begin{aligned}
 P_0&=N_{12},\\
 P_1&=\operatorname{tr}(C^3DCD),\\
 P_2&=\operatorname{tr}(C^3D^3),\\
 P_3&=\operatorname{tr}(\theta^5H_{16}).
 \end{aligned}
\]

追加した旧trace6は一つだけ、word000001。C/D-onlyの8本の独立性を保持してこれが9本目となった。

16次はF₄Inv₁₂の9本とSym²Inv₈の15本のspanがrank19。保存済み全36 C/D trace8を加えてもこのspanを補えず、許可された旧θ/H trace8に進んだ。選択した追加成分は

\[
 T_{16}=\operatorname{tr}(\theta^7H_{16})
\]

のword00000001で、rank20となった。

24次の候補70本は以下である。

1. F₄とp₁,…,p₄だけからなるweighted degree24の35 monomial。
2. Pᵢ×(F₄³,F₄p₁,…,F₄p₄)、20本。
3. PᵢPⱼ (i≤j)、10本。
4. T₁₆×(F₄²,p₁,p₂,p₃,p₄)、5本。

|検査|係数行数|rank|選択minor mod32003|
|---|---:|---:|---:|
|12次|91|9|20149|
|16次product|153|19|別途全係数を保存|
|16次全体|153|20|2363|
|24次70列|325|67|26663|

係数行はsⁱtʲをi=0,…,d、j=0,…,d−iの順とする。24次の選択行は0,…,66。選択列は0,…,69から**62,68,69を除いた67列**である。全70列、選択行列、word定義は保存した。全係数のPython消去でも同じ選択列・row・minorを得た。

## QQ上で全空間を張る理由

全候補はSp₄の普遍的なhomogeneous不変式で、β次数24。元のρ、Ω、Gramと逆行列の全分母は32003で単元である。従って非零のmod32003係数minorは、選択した67本がQQ上で線形独立であることを証明する。一方、独立に確定済みのcharacter計算がdim Inv₂₄=67を与える。この二つから**選択67本は全Inv₂₄のQQ基底**となる。12次9、16次20についても同じ理由で全空間である。

70本には普遍的な線形関係が3本ある。将来source上の関係を求める際、70列の任意のkernel vectorをそのままNEW25の方程式と扱ってはならない。ここで固定した67本の基底を使うか、別に普遍的非零性を確認する必要がある。

## source展開の未実行範囲

既存のF₄,pᵢ,N₁₂に加えて必要な新しい低次数source出力はP₁,P₂,P₃,T₁₆の四つ。source cを代数的に1とした後、homogeneityで復元するc次数はそれぞれ6,6,6,8となる。P₃とT₁₆の選択wordについてpreflightのcoefficient majorantは499bit、665bit。二つのC/D trace6の既存最悪上界と合わせたraw JSON上界は34,527,115 bytesで、32MiB内の完走保証にはならない。実際の係数supportを計算した値ではない。

後続source展開と67本のcoefficient relation検査には、別の明示的なbounded budgetが必要。**このpacketはNEW25上の関係・次数24方程式の存在・24次の排除のいずれも主張しない。** actual orbit invariant fieldの有限index問題も未解決のままである。

## 再現と保存

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_full_degree24_hessian_basis_20261104_v63.py \
  --output work/full24-hessian-replay-report --work work/full24-hessian-replay-native
```

output/workには新しいdirectoryを指定する。nativeは一回、CPU・wall60秒/2GB/certificate32MiB上限。実測native約0.336秒、wrapper約3.89秒。最大support325、scalar polynomial products130394回、保存全体約157KB。

一次保存：`reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63/`。
`formal-input.txt.gz`、`native-certificate.json.gz`、`definitions.json`、`minor-witnesses.json`、`report.json`を含む。producerは同名のPython/C++、C++は凍結済みC/D-only producerのpolynomial/rank helperをincludeし、そのmainは実行しない。全旧source・character・C/D-only証拠と定数入力をhashで固定した。

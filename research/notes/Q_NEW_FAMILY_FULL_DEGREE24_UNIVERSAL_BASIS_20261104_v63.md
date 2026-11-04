# 全24次不変式基底への一つの形式平面検査：12次 gate で停止

2026-11-04。`Q_NEW_FAMILY_FULL_DEGREE24_PREFLIGHT_20261104_v63` の一つの有限な案を実行した。**12次の候補係数rankが8に留まり、必要な9に届かなかった。全Inv₂₄基底は未構成であり、NEW25の方程式・次数24の可否について新しい結論はない。**

## 入力と厳密な計算範囲

凍結済み `Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json` の最初の三つの ambient β vector を、そのまま整数liftして

\[
 \beta=\beta_0+s\beta_1+t\beta_2
\]

とした。固定素数は32003。一つの形式的ambient二平面での**多項式係数**を比較したもので、actual NEW25の点・追加のsource点・点評価行列は使っていない。三vector、元の136個のquadratic monomial順序、C/D全32×136係数を入力に保存した。C/D係数は4倍すると整数であり、4は32003で単元。

既存のSpin5 packetで定義したC,Dをその係数から直接構成し、全14個のtrace6 necklaceを計算した。F₄とdet(tC+uD)の五係数も直接構成し、

\[
 \operatorname{tr}C^2=2\operatorname{tr}D^2,
 \quad \operatorname{tr}CD=-\operatorname{tr}D^2,
 \quad F_4^2=400p_0+1600p_4
\]

を形式多項式として再確認した。N₁₂は

\[
 N_{12}=\operatorname{tr}(CD^2CD^2)
 -\tfrac12\operatorname{tr}(D^2)\operatorname{tr}(C^2D^2)
 +p_4\operatorname{tr}(C^2)
\]

で構成した。Inv₈の順序は(F₄²,p₁,p₂,p₃,p₄)。

## 実際の gate 結果

比較した12次20列は、F₄Inv₈の5列、N₁₂、全14個のC/D trace6である。sⁱtʲ (i=0,…,12、j=0,…,12−i) の91係数を用いた。

|列の範囲|係数rank mod32003|
|---|---:|
|F₄Inv₈|5|
|さらにN₁₂|6|
|全trace6を追加|8|

選択列は0-basedで[0,1,2,3,4,5,9,10]、選択係数行は[0,1,2,3,4,5,6,7]。8×8 minorは **19943 mod32003** で非零。追加の二つのwordは000101と000111である。Python側でも全91×20係数の消去とminorを再計算して一致した。

characterで証明済みのdim Inv₁₂=9に対し、この一つの形式平面でrank8だったため、以後のproduct16 rank19、Inv₁₆ rank20、70列のInv₂₄ rank67というgateには進まなかった。native実装はgate判定前に全36 trace8の153係数も前計算・保存しているが、**16次のrank検査・24次の列構成はいずれも未実行**。別平面、候補増殖、source展開、再試行は行っていない。

非零minorはこれら8本がQQ上でも独立であることを証明する。**一つの平面でのrank8は、普遍的なtrace6空間の次元上界8を証明しない。** 普遍的に9次元でもこの平面への制限が退化した可能性があり、候補そのものが不足している可能性もある。どちらかは未判定。

## 保存と次の有限な問題

実行はnative一回だけ、約0.016秒、wrapper約1.83秒。CPU・wall60秒、address space2GB、certificate32MiBの上限内。全係数certificateは圧縮約22KB。

再現コマンド（output/workは新しいdirectoryを指定）:

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_full_degree24_universal_basis_20261104_v63.py \
  --output work/full24-basis-replay-report --work work/full24-basis-replay-native
```

一次保存は `reports/Q_NEW_FAMILY_FULL_DEGREE24_UNIVERSAL_BASIS_20261104_v63/`。
`formal-input.txt.gz` が三vectorと全C/D係数、`native-certificate.json.gz` が全低次数係数・実際の選択行列、`minor-witnesses.json` が非零minor、`definitions.json` が候補word、`report.json` が入力・証拠のhashとscopeを保持する。`all_passed=true` はこの有限計算と記録の検証が成功したことを表し、`universal_complete_QQ_Inv24_basis=false` と `status=degree12_rank_gate_failed` が数学的な停止点である。

残る有限な問題は、まず**全Inv₁₂の9次元を独立に確保し、既存N₁₂を含む4個のprimitive12を選ぶこと**。旧θ/Hから9次元基底は既知だが、それをC/D trace6のみへ置き換えられるとは今回仮定も証明もしなかった。全Inv₂₄基底が得られるまでは、41部分空間の単射結果を全67次元へ拡張してはならない。

# 全Inv₂₄基底に必要な四つの低次数source不変式：全係数展開

2026-11-04。**C/D word二本と旧16×16 θ/H word二本のsource上の全係数を、native一回で厳密に展開した。** 既存のF₄,pᵢ,N₁₂と合わせれば、`Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63`で固定した全67本をsourceへ制限するための低次数入力が揃う。67列の関係式検査は本packetでは行っていない。

## sourceとT²=Bによる前処理

凍結済み六変数sourceのcを代数的に1とし、δ-constantを一意に記録できる既存t=2 section

\[
 v=-98304X/H,\qquad H=147456X^2-(3072W-1)^2
\]

を使う。元の整数β numerator b=1080BHβは、先の独立検証済み低不変式展開の入力をそのまま読み戻した。新しい点は作らず評価もしない。

元のρから第二Cartan cocharacterを取り、その16個のweight wᵢを対角行列として再検証した。直接B^((wᵢ−1)/2)を使うと、nonzero座標1,7のweight2で半整数冪が出る。そこで**T²=B**を使い

\[
 \bar b_i=T^{w_i-1}b_i(B=T^2)
 =T^{-1}(\operatorname{torus}(T)b)_i
\]

とした。B⁰のβ座標はweight−1で、bᵢに既にB=T²因子があるため全入力は多項式。B⁻¹座標のweight1は不変、weight2はTを一つ持つ。実際にnonzero入力のT次数は0か1であることを全項で確認した。

β次数dの不変式に対し

\[
 I(\bar b)=B^{-d/2}I(b),\qquad
 I(\beta)=\frac{I(\bar b)}{1080^d H^d B^{d/2}}.
\]

これはtorus不変性とhomogeneityの恒等式であり、新しいsource仮定ではない。既知の余分な共通B因子を行列冪の**前**に取り除くため、中間係数を抑えられる。最後に全出力でTが偶数次数だけを持つことを確認し、B=T²へ戻した。

## 整数vector recursionと規格化

C̃=4C(b̄)、D̃=4D(b̄)を全32×136係数から構成した。μ̂=2μ(b̄)=[(−C̃+2D̃)/2]ᵀの整数性とJ⁻¹μ̂の対称性をsource多項式として検証し、ρの元の基底でΘ=2θ(b̄)を作った。Ωとρは整数、G=2Gram⁻¹も整数である。

\[
 Z=(\rho(A_j)\bar b)_j,\quad V=\Omega Z,\quad Z_{k+1}=\Theta Z_k
\]

を二levelだけ保持してk=7まで進め、k=5,7で

\[
 \operatorname{tr}(\theta^kH_{16})(\bar b)
 =\frac{\operatorname{tr}\Theta^{k+1}+
 \operatorname{tr}_{10}(G V^T Z_k)}{2^{k+1}}
\]

を用いた。純trace6/8項には、事前にweight identityと同じ保存済み形式平面の全係数で確認したs₁,s₂の式を使った。H₁₆自身やその高い冪を密に保存していない。

C/Dの二本は短い4×4 prefixからtr(C̃³D̃C̃D̃)、tr(C̃³D̃³)を得て4⁶で割る。旧word二本は上式の2⁶、2⁸で割る。さらに1080ᵈ、Bᵈᐟ²、Hᵈを追跡し、整数contentと追加のB因子だけを取り除いた。H因子は取り除いていない。

全出力でuとXの次数和も偶数であることを検証し、a=u/X、Y=X²としてuⁱXʲ→aⁱY^((i+j)/2)を行った。従って出力はQQ[a,Y,W,B]の整数numeratorと明示分母である。sourceのcを戻す次数はd/2。

## 完成した全係数

各JSONの`terms`は[整数係数文字列,a次数,Y次数,W次数,B次数]。関数値はこの多項式を`denominator_integer * B^denominator_B_power * H^denominator_D_power`で割ったもの。

|出力|不変式|項数|係数最大bit|整数分母|B分母|H分母|復元c次数|
|---|---|---:|---:|---:|---:|---:|---:|
|CD000101|tr(C³DCD)|3449|382|340122240000000|6|12|6|
|CD000111|tr(C³D³)|3618|382|680244480000000|6|12|6|
|theta5H|tr(θ⁵H₁₆)|3948|388|136048896000|6|12|6|
|theta7H|tr(θ⁷H₁₆)|10819|517|881596846080000|8|16|8|

全raw係数保存は2,156,304 bytes、圧縮係数は約0.97MBで32MiB以内。上界が僅かに超えていたpreflightを、実際の項数が大きく下回った。vector recursionの各levelで中間150k項、出力100k項のcapを監視し、全T parityとu/X parityを全項で確認した。nativeは一回だけ、wrapper約10.62秒で完了した。

## 再現・証拠・未完了点

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_hessian_low_invariant_expand_20261104_v63.py \
  --output work/hessian-source-replay-report --work work/hessian-source-replay-native
```

output/workは新規directoryを指定する。hard capはnative CPU・wall60秒、address space2GB、全raw32MiB。完成順に出力を保存し、自動再試行をしない。producerは同名Python/C++。C++は既存整数FLINT helperをincludeし、旧mainは実行しない。

一次保存は`reports/Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63/`。四つの圧縮係数JSON、完全なnative入力、`recipe.json`、native log/report、入力・証拠hashを含む`report.json`を保持する。

この計算は**低次数source全係数の供給**である。全Inv₂₄の67本のsource関係、NEW25の24次方程式の存在または排除、67列のcoefficient rank、actual orbit-field有限indexの解消はまだ主張しない。後続は固定済み67本の普遍基底を用い、70候補の普遍関係3本と混同しない別の有限な検査となる。

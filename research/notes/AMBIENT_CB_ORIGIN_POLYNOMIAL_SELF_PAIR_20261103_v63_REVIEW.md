# 新しい CB 多項式 self-pair：全 source・二枝・非 CM の独立監査

2026-11-03。対象は `Q_b=P0−b(2C[1]+C[9])` という新しい族である。先行するCAの `P_a=P0−a(3C[1]+C[9])` と区別する。固定旗20基底Cは同じであるが、実際のCB sourceは別に検証した。

## 1. 実 tensor と二つの CB source

独立コード `code/ambient_CB_origin_polynomial_self_pair_review_20261103_v63.py` は作者の raw / canonical / action / push helperを import しない。raw CA20基底を物理式から再構成し、P0とQ_bの全150係数を照合する。実際のvector場は

\[
V_b=(x^2+8bhx/5)\partial_h-2bx^2\partial_x/5
 +(hx+3bxz/5)\partial_z+(h^2+z^2-7bxw/5)\partial_w.
\]

`Q_b=∂y∧V_b` であり、25 divergenceと全350 strict cubic Jacobi係数はQQ[b]で零である。Euler∧Qのx³係数は1なので、各parameterの構造の一般階数は正確に2である。

二つの実行列は

```
G1: (h_new,x_new,y_new,z_new,w_new)=(x,x+bh-b²z,y,h,w), det=-b²;
G2: (h_new,x_new,y_new,z_new,w_new)=(x+2bh,x,y,z,w), det=2b.
```

QQ(b)上の両逆行列を確認し、全5×5 tensorの有限pushforward `G_i Q_b(G_i⁻¹X)G_iᵀ` を直接計算する。その全150係数は、標準CB式 `canonical(rawCB(f_i,0,H_i))` に一致する。f1=−bhx、f2=2bx²、二つの有限H cubicは全係数で作者証明書と一致する。

元の旗のmarked lineはxとx+2bhであり、b≠0で相異なる。L=e_yは共通で、W1=span(x,x+bh−b²z)、W2=span(h,x)。これは異なるchart名だけの比較ではなく、同じQ_b上の異なる実旗である。

## 2. full ratio source と全接方向

普遍的な20変数sourceを

\[
\Phi(f_3,H_{16},\mu)=F(f)+\mu G(f)+H
\]

として作り、**μfを残したまま微分する**。三つの mixed μ/f微分がG3の全150係数に一致することも確認する。μ=0の時点でratio列を捨てたり、独立なbinary g3に置換したりしない。

全9旗方向の物理的 Lie 微分を加えた両150×29 Jacobianが保存値と一致する。固定19方向の正確な左逆でquotientを取り、各150×10行列のminorを自分で計算した。

\[
432b^5/15625,\qquad 786432b^{10}/15625.
\]

従って**全b≠0**で両sourceはrank29のimmersionである。全GL25列のQQ(b)階数は19、scalarを追加しても19、dQ/dbを追加すると20となる。保存minorと全係数を照合したので、この族のG-saturationの実像はcone20、射影19次元である。

## 3. 同一Qでの全58接方向

独立計算では、**両方の全29方向**をそれぞれG1⁻¹,G2⁻¹で元のQ_bへpushbackし、全150×58の階数が37であることを確認する。さらに作者の第一frameにおける30列quotientも全係数で再構成し、rank18とminor

\[
82944b^{25}/15625
\]

を確認する。従って全b≠0で二つのsource接空間の和は37、交わりは21次元である。37をwhole imageのembedding dimensionの等号と同一視しない。cleanなscheme intersectionもこの線形計算からは従わない。

## 4. 全proper fiberからの一般bの完成環

既受理の全proper CB旗模型と、P0での**全flag fiber一台長2**を明示依存とする。証明書は `AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW` と `AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW` であり、その既存GBを再実行していない。

propernessと有限なP0 fiberから、P0を含むtarget openで全source mapはfiniteになる。finite moduleのfiber長の上半連続性により、さらに全fiber長≤2の開を取れる。多項式target曲線Q_bの逆像はb=0を含む非空開である。そのb≠0には二つの実旗が既にあり、全fiber長≤2なので、**全fiberは正確に二つの被約点**になる。新しい二sourceの行列にbの逆数が現れても、このproper argumentはtarget曲線全体とP0 fiberを使っているので問題ない。

この主張はb=0を除くある非空開とそのG-saturationで成立する。全b≠0の全旗分類は、このpacketでは行っていない。

sourceはsmooth29で、実source mapのrank29により、ambient formal ring Tから各source completion B_iへの写像は全射となる。`B_i≅C[[29]]`、kernelをI_iとする。二つの接空間の和の次元37>29だからI1≠I2である。

全finite sourceのreduced imageは実CBであり、その局所環から全source代数への単射はflatなcompletionで保存される。全fiberの二点を全て用いたため

\[
\widehat{\mathcal O}_{CB}=T/(I_1\cap I_2),
\qquad T/I_i\simeq\mathbb C[[29]]
\]

となり、第三のformal branchを落としていない。global normalizationの未検証なbirationalityを先に仮定する議論ではない。射影ではsmooth28の二つのformal branchである。

## 5. meeting schemeとdepthの範囲

`C'=T/(I1+I2)` の接空間は二つの接空間の交わりなので、edimは21である。G×bの実paired-sourceが共通像20次元を与え、20≤dim C'≤21となるが、schemeの正確な型・clean性・conductorの方程式はまだ決まらない。

完全列

\[
0\to T/(I_1\cap I_2)\to B_1\oplus B_2\to C'\to0
\]

でdepth(B1⊕B2)=29、depth C'≤21だから、depth lemmaより個別reduced CBのdepthは22以下、次元29より小さい。従ってこのG-saturated族の一般点で、CBはsingular・nonnormal・non-CMである。射影のdepth上限は21である。

既存の原点quadratic base上の表示 `A=R⊕It⊂R[t]/(t²−ν)` に移す場合、この二つの被約fiberを持つpunctured self-pairの点ではν≠0となる。étaleに二sheetへ分裂させると `A=R×_(R/I)R` であり、meeting schemeがR/Iになる。この事実を原点conductorの新しい下界へ使うには、別のspecialization・multiplicity論証が必要である。

**原点conductorのmultiplicity4、全odd ideal、CM性、B18原点の完成型、actual CA/CB intersection Dは、この独立reviewから結論しない。**

## 再現と保存

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/ambient_CB_origin_polynomial_self_pair_review_20261103_v63.py \
  --output work/cb-origin-self-pair-independent-review
```

`reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW/report.json` に、全source・minor・全58接方向と論証の適用範囲を保存する。作者の新しい全b≠0 companionがある場合は、その全旗被覆を別に監査する。

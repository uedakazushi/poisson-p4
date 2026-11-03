# Q_(8/9)の明示T38フレームと二つの固定valuation検査：独立監査

次の三つの凍結packetを、以下の範囲で受理する。

- `Q_T38_QA8OVER9_COMPACT_FRAME_20261103_v63`：四次体上の実際のT38 source curveと全150係数のQ_(8/9)極限。
- `Q_T38_QA8OVER9_EDGE_WEIGHT_TEST_20261103_v63`：固定frame・level−1・定数edge amplitudeというansatzの限定的な否定。
- `Q_T38_QA8OVER9_SEED_WEIGHT_TEST_20261103_v63`：同じframe・固定period normalization・level5の定数seedについて、同じ物理的Q_(8/9)を保つidealが零seedのmaximal idealであること。

P+、P2、正weight6軸のT38所属や、T38全閉包からの除外は、いずれも結論しない。独立検算は0.60秒で終了し、大規模GBは使っていない。

## 1. 独立した四次体演算とadmissibility

主producerのSymPy algebraic fieldやその行列・tensor関数をimportせず、各元を四つの有理数のvectorとして実装した。

\[
 K=\mathbf Q[\theta]/(\theta^4-50\theta^2+305),\qquad
 \theta^4=50\theta^2-305.
\]

積は上の式で次数を下げ、逆元とnormはQQ上の4×4乗算行列からGaussian消去で求める。5-Eisenstein条件を確認したため、これは積代数ではなく次数4の体である。実平方根の一つを数値的に選んで等式を判断したものではない。

rho,z,x,y、alpha,beta,s,muの全式を再計算し、rho²−rho−1=0、`2rho(z²+1)+7z=0`、a=8/9を確認した。五固有値の全10差について保存逆元とQQ normを別演算で再計算し、さらにgamma12とPfaffianの逆元・normも確認した。s、s+3、muも非零単位である。従って全subdiagonal link、固有frame、punctured periodの分母は正当に使える。

## 2. T38 periodと元の全sourceへのbinding

gammaの全六成分を式から再構成し、保存値と一致することを確認した。gamma12=pとして

\[
 (\gamma_{34},-\gamma_{24},\gamma_{23},\gamma_{14},-\gamma_{13},\gamma_{12})
     =p(0,1/2,0,0,-1/2,1),\quad
 \operatorname{Pf}(\gamma)=p^2/4\ne0.
\]

これは記録されたperiod conventionの **Pfaffian dual** であり、gamma自体をperiod basisに代入したものではない。rhoとVandermondeによるPfaffianの式も別に確認した。固有値のlabel permutationは用いていない。

さらに `MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json` の **全homogeneous bivector** を、p0を掛けてQQ多項式として比較した。toric三項、四つのlinear edge、三つのquadratic cross correctionの全十項が一致し、他の項を省略していない。toric seedは正確に

\[
 p_0x_1x_2\partial_{12}+(p_0/2)x_1x_3\partial_{13}
       -(p_0/2)x_2x_4\partial_{24}
\]

で、記録されたopenはp0²/4≠0である。この元sourceのT38への帰属は既受理のMPS source証明を依存として保持する。任意のT01退化曲線をT38と呼び替えたものではない。

## 3. 実際の150曲線・固有frame・物理式

N_tの全25 entry、四つのtracefree冪W1,…,W4、および

\[
 \Pi_t=\mu(\alpha t^2W_1\wedge W_2+\beta tW_1\wedge W_3
                         +\rho W_1\wedge W_4+W_2\wedge W_3)
\]

を独立なK[t,t^-1] sparse演算で構成した。主側のwedge関数も使用していない。全150係数が保存されたK[t]の曲線と一致し、負のt冪が無いことを確認した。全150個のt=0係数は物理canonical Q_(8/9)と一致した。

同じ150係数から五つの物理変数を直接微分・乗算し、**全25 divergence係数と全350 strict Jacobi係数** がK[t]で厳密に零となることを確認した。projective Jacobiだけの検査ではない。

S0をrecurrenceから再構成し、`S_ij=S0_ij t^(j-i)` とその有限nilpotent逆を計算した。両逆積、N_t S=S t diag(d)、四つ全てのtracefree冪の対角化を確認した。さらに三項log seedのEuler projection

\[
 P_{can}=P_{raw}-(Euler\wedge\operatorname{div}P_{raw})/5
\]

を独立に作り、**このcanonical tensorを実際にSでpushforwardし、mu t^5を掛けた全150係数がPi_tに等しい** ことを確認した。従ってt≠0では旧T38 sourceのedge=0部分と実際の可逆なprojective線形変換で同定される。source period p0=mu p t^5はそこで単位である。t=0のS自体の正則性を仮定せず、物理曲線の正則性を直接確認している。

## 4. 全150×4 edge行列とlevel−1の限定的否定

四つの元edge tensorを実S0でpushforwardし、保存された150×4行列Eの全600個の四次体係数と比較した。四本のdivergenceは零である。三つのcross tensorも直接pushforwardしてdivergence零を確認したため、canonical projectionによる未記録のedge/cross項は無い。

SとS^-1の全entryが同じ指数則 `t^(j-i)` を持つことを確認した。各pushforward項の指数を足すと

\[
 t^{a-i+b-j+u-p+v-q}=t^{w_j-W},\quad
 w=(-5,-2,-4,0),\quad W=i+j-p-q
\]

となる。従ってこれは全150 rowに成立する指数則であり、有限個のrowからの推測ではない。

W≥0の実rowから選ぶ4×4 minorを独立な四次体Gaussian消去で計算し、保存された値・逆元・QQ normと一致する非零単位を得た。各edgeの最高Wは実際に7であった。

固定ansatz `tau_j=c_j t^(-1-w_j)` では、同じtarget rowのtoric、linear、cross部分のt冪はそれぞれ5−W、−1−W、−7−Wとなる。三つは異なるので、W≥0でのlinear poleを他の二部分は相殺できない。Eのrank4から全定数amplitude c_j=0が必要となる。

受理するのはこの **固定frame・固定period normalization・定数aligned amplitude** の必要条件だけである。非aligned valuation、t依存frameの摂動、別period/固有値seed、多段階のedge調整は扱わない。

## 5. 同level5 seedの全8×7係数証明

`tau_j=c_j t^(5-w_j)` では、period p0=mu p t^5を含めた三cross correctionもtoric・linearと同じlevel5になる。全sourceの各rowはt^(5−W)倍となる。regularで同じ **物理的** Q_(8/9)に到る条件は、W≥5の全8 row

`90,105,120,121,135,136,137,140`

について追加係数が零となることである。

この八式の七monomial `(c0,c1,c2,c3,c0c3,c1c2,c2c3)` に対する **実8×7係数行列** を、上で独立に作ったEと三つのcross pushforwardから再構成し、全56 entryを照合した。最初の七式の7×7行列式・逆元・QQ normを再計算し、保存逆行列の両側積も単位行列であることを確認した。

従って七monomialそれぞれが七式の定数K係数線形結合であり、特にc0,…,c3がそのidealに入る。逆に全八式は(c0,…,c3)に入るので、全八式のidealは **正確に(c0,c1,c2,c3)** である。根基の一致やGB標本からの推論ではない。

これは固定frameで零以外のconstant seedから同じQ_(8/9)を出発点にする経路を排除するだけで、T38全閉包、P+、P2、他のMPS成分の所属・非所属を決めない。

独立codeは `code/q_T38_Qa8over9_review_20261103_v63.py` と `code/q_T38_quartic_arithmetic_review_20261103_v63.py`。reportと小certificateは `reports/Q_T38_QA8OVER9_20261103_v63_REVIEW/` に保存した。三つの主packetの全input/evidence hashを照合し、producerのfield/tensor計算を呼ばずに検算した。

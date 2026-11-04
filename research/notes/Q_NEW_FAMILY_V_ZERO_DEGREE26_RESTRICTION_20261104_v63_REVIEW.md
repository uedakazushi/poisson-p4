# Dense v0切片の全Inv26制限：独立係数レビュー

2026-11-04。**完全Inv26の制限が単射、NEW定義式の次数下界28を受理する。**
28次方程式の存在、source28計算、新しい点・prime・候補追加は主張しない。
authorのFLINT exact-QQ workerを再実行せず、独立Python整数辞書でGF32003係数を再構築する。

## Actual sourceと完全基底

受理済みv0 sectionは同じactual B_NEW15へのSp4飽和がdenseである。
その16式へxi0=1,r0=R⁴,B0=T²を入れ、全Cartan weightsから
`diag(T/R,R/T,T²R²,T^-2R^-2)` の有限pullbackを独立に計算する。
Rの指数をH=R²,L=R³Tへ変換すると、保存された九非零成分/十二項と全係数が一致する。
これは新しい点評価ではなく、formal四変数多項式である。

各非零beta_iのxi0指数e_iは `2e_i=1-w1_i-w2_i` を満たす。
degree d・Cartan weight0の不変式ではxi0総指数がd/2になるので、26次では正確にxi0¹³を戻せる。
xi0=1へ置いたことが同次不変式のvanishingを失わせない。
補助R,Tの有限被覆はdenseで、Sp4 pullbackは不変式を変えない。

完全45基底は別の独立レビューが全16-variable Weyl character上界45と係数下界45を一致させている。
今回基底を追加・再選択せず、全57候補のsemantic順と選択45 indexを固定して使用する。

## Ambient微分を先に保つ

元136個のbeta_i beta_jに対するC,D全32二次tensorから値と512微分を独立に作る。
beta_iが切片上で零でも、そのambient偏微分は捨てない。
全16方向を保持してから十二項betaを代入する。Euler式も全32で検査する。

五trace4、三trace6についてcyclic derivativeと全global quadratic微分を縮約し、
全128 gradientを保存QQ多項式の還元と一致させる。各Euler式も全係数で確認する。
四つのsource座標に関するchain-rule微分をPoisson bracketへ代用していない。
bracketは元16×16 Omega逆行列で縮約する。

さらにOmega/rho/trace Gramからtheta全256 quadratic係数を再構成し、
元の四次moment polynomialの全ambient二階微分から本来の16×16 H16を作る。
native inputの全512×136正規化係数と全Omega逆係数を独立にbindする。
H16を4×4随伴Dと取り違えない。
tr(theta⁵ H16)も独立多項式積で計算し、Inv12の九基底を正しく使用する。

## 全45列と非零係数minor

low19値、trace8値、ambient128勾配、全45制限多項式について、保存された
`QQ content × primitive integer coefficients` のGF32003還元と全係数一致を確認する。
全分母が32003の単元であることも検査する。
packed exponentはbase64で、pipeline最大ambient degree26、各beta変数の各source指数≤1なので
各指数は26以下。整数加算にcarryは起きない。

保存QQ全22,674項のL偶性と2H+3Lの4合同を確認し、
`H^a L^b→r0^((2a+3b)/4) B0^(b/2)` へ下降する。
この指数写像は許された格子で単射なので、異なるH,L項を勝手に同一視していない。

選択された同じ45個の係数monomialについて、全2,025成分を独立に再構築する。
別の整数modular eliminationでdet=7069≠0。
これは点評価matrixでなく、formal source多項式の係数matrixである。
従ってQQ上でもそのminorは非零で、完全Inv26の制限は単射となる。
exact QQ全係数を独立に再展開したという主張ではなく、p-unitなglobal定義からの
独立modular係数構築が非零QQ minorを証明している。

## 下界28の範囲

先行する24次までの完全除外と合わせ、26次不変式もNEW上で消えない。
中心-Iは全16 weightsで-1として作用するため奇数次数不変式は零。
actual NEW像は既受理のhomogeneous Sp4安定な既約hypersurfaceであり、Sp4はcharacterを持たない。
従ってその同次定義式は不変式であり、degreeは少なくとも28となる。
degree28で実際に消える式があることや、その正確なdegreeを決めたものではない。

再現:

```
python code/q_new_family_v_zero_degree26_restriction_review_20261104_v63.py \
  --output /tmp/v0-Inv26-independent-review
```

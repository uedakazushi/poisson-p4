# 楕円 Cremona 変換の primitive 一次ジェット恒等式：独立点検

2026-11-01。elliptic_jets 担当の普遍恒等式を、保存済み Jacobi 行列・余因子・M2 出力を
読み込まずに、五つの二次式から独立に再計算した。
コード：`code/elliptic_cremona_review_20261101.py`。
証明書：`reports/ELLIPTIC_CREMONA_REVIEW_20261101.json`。

## 1. 全パラメータで確認した式

添字を mod 5 とし

\[
 Q_i(x)=\lambda x_i^2+\lambda^2x_{i-1}x_{i+1}-x_{i-2}x_{i+2}
\]

とする。f=[Q₀:…:Q₄] は二次 Cremona 写像。独立コードは以下を有理係数の多項式として
確認した。有限点の標本ではない。

1. 既存の q(5,1),q(5,2) 族に対して、射影 Poisson 写像の関係 f_*P₁=−P₂ が成立する。
   斉次 lift の radial 項を許容する正しい確認として、残差行列と Q の wedge の10成分を
   全て検証した。
2. K=λ³P₁、F=Pf(K)、N=F(−K⁻¹) と置くと、det(DQ)|_{x₀=1}=2F。
   余因子は二次式から再計算し、DQ adj(DQ)=2F I₅ も確認した。
3. Γ_hom=A_y(DQ)⁻¹ Hess(Q) の接続項を primitive 射影したものは
   (3/2)β₁。ここで β₁ は −dΩ₁ の primitive 部分であり、64 tensor 成分を
   共通分母 F² を払って検証した。

## 2. 斉次 Hessian とアフィン Hessian の差

アフィン座標で f_i=Q_i/Q₀、ℓ=d log Q₀ とする。Euler の斉次恒等式を使うと

\[
 ((Df)^{-1}D^2f)^c_{ab}
 =\bigl(A_y(DQ)^{-1}\operatorname{Hess}(Q)\bigr)^c_{ab}
  -\ell_a\delta_b^c-\ell_b\delta_a^c.
\]

最後の二項は projective perspective の接続項であり、primitive 射影で消える。
従って斉次式を用いた上の検証は、実際のアフィン一次ジェットの変換を検証している。

f_*P₁=−P₂ から f*Ω₂=−Ω₁。これを微分すると、tensor としての pullback は

\[
 f^*\beta_2=-\beta_1+
 \operatorname{prim}(-\Gamma^t\Omega_1-\Omega_1\Gamma)
 =\tfrac12\beta_1.
\]

ここで非線形な Poisson 同値一般が primitive 一次ジェットを保つとは主張していない。
この特定の二次写像の Hessian に対する普遍恒等式が不可欠である。

## 3. 正規化と次数

source の frame F₁ は cP₁(p) を J に、target の F₂ は c′P₂(f(p)) を J に正規化する。
c′=−4c とし

\[
 S_0=-2F_2DfF_1^{-1}
\]

と置く。この F₂ は **c′ を含めた** frame である。
f_*P₁=−P₂ により S₀JS₀ᵗ=J が成立する。
上記の 1/2 の pullback 係数と三つの tensor 添字の −2 による変換が一致するため、
二つの正規化された primitive 一次ジェットは S₀ で結ばれる。

従って (λ,p,c,S) を (λ,f(p),−4c,S S₀⁻¹) に送る写像は、両表示空間から W16 への
写像を可換にする。f の有理逆が別途証明されれば、これは W16 上の双有理同型である。
f の birationality は古典的 quadro-cubic 変換であり、一次資料では Polishchuk,
[Feigin–Odesskii brackets, syzygies, and Cremona transformations](https://arxiv.org/html/2204.09152v2),
Theorem B が対応する一般の構成を与える。具体的逆式の証明書は generic_fiber 担当が
別に検証する。この依存を本チェックだけで省略しない。

その逆式を合わせると、一般表示ファイバー長 N₃₉=N₄₀。
既に独立に証明した ρ₃₉=ρ₄₀=3000 と d₄₀=1 を用いて

\[
 d_{39}=d_{40}=1
\]

を得る。T01 の次数3と合わせ、既知三成分上の一般延長は3+1+1=5個になる。
未知の支配成分や全 S61 の有限性を排除したという意味ではない。

## 再現

```sh
.venv/bin/python code/elliptic_cremona_review_20261101.py \
  --output /tmp/elliptic-cremona-review.json
```

Python の多項式環によるこの検証は、主実装の M2 恒等式検証とも独立に一致した。

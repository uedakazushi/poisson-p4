# 連結 Spin4 の joint biform 不変関数体：独立検証

2026-11-04。対象は複素数体上の
`Sym²(A)⊗B ⊕ A⊗Sym²(B)` と `SL2(A)×SL2(B)` の不変**関数体**。
既受理の物理 Spin4 切片との表現同定を入力とし、著者 producer を再実行せず
全係数作用と完全性を検証する。

一般 U の二つの quadratic p,q は resultant 非零で、Jacobian quadratic の
discriminant はその16倍である。独立コードでは著者の Sylvester 行列式と
違い、一変数 subresultant 計算で全六係数恒等式を確認する。従って Jacobian
を `4c²x0x1` に SL2(A) で正規化でき、両端条件から p,q の mixed 係数は零。
ここには `b*detM=a*(bg−df)+d*(af−be)` 等の正しい plus 符号の恒等式を使う。
detM=c² より SL2(B) の変換で `U=c(x0²y0+x1²y1)` を得る。正規化の有限な
root 選択は次の**全 normalizer**で回収し、見落とした有限指数を仮定しない。

U0 line を保存する A は span(x0²,x1²) を保存する。その mixed 係数が零で
detA=1なので A は対角または反対角のみ。誘導二平面の determinant は各々
1,−1で、補償する B を完全に解くと U0 の倍率 ζ は各々 ζ²=1,−1。
ζ=1 の exact stabilizer は T のみ。四つの line-normalizer coset は
`T,σT,σ²T,σ³T` で尽くされ、σ⁴=1、σ²∉T、σtσ⁻¹=t⁻¹。
従って cross-section の N 商は一般 G 軌道を**完全に**記述する。

全 V 係数の重みと σ permutation を元の二変数ずつの polynomial 置換から
別に算出する。weight1 の v1 と五つの torus 不変座標から全 V を有理的に
回収するので、この段階の invariant field の同一性には有限指数が残らない。
σ² は c の符号だけを残す。z=c² として二つの paired block を取り直すと、
残る involution は
`(z,r,P,Q',S,T'')→(-z,-r,Q',P,T'',S)` となる。

著者の六つの不変式 I0,…,I5 を全係数で確認し、その Jacobian は独立な六
変数上で `8/z²`。z²=I0 を付加すれば全六変数が明示的な逆式で戻り、残る
involution は非自明な z→−z のみ。よって固定体は厳密に C(I0,…,I5)、六つは
代数的独立である。単独 biform の定理を直和へ流用した議論ではない。

これで**連結 Spin4 の joint biform quotient 関数体**の有理性を受理する。
元の W16/Sp4 の別の normal-line component は U,V を交換し、その追加の
作用は今回計算していない。従って W16/Sp4 の有理性、全多項式不変式環や
global quotient、NEW25 の方程式は結論に含めない。新しい点・高次数 probe・
GB は使わない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_firstjet_spin4_biform_field_review_20261104_v63.py --output work/spin4-joint-field-independent
```

60秒・2.2GB上限。著者の修正済み凍結コード・note・report のハッシュを束縛する。

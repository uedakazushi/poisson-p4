# Z₁,Z₂係数のfull content：891項のprimitive候補

2026-11-04。既知境界を除いた6,919項のQについて、**(Z₁,Z₂)係数15群のfull contentをexact FLINT gcd/cofactor計算で求めた**。全係数の等号

\[
 \boxed{Q=C(r,R)\,G(r,R,Z_0,Z_1,Z_2)}
\]

を確認し、商の全係数gcdが1であるというnative最大性検査も成功した。Cは70項・r次数9、Gは891項・r次数10。独立な最大性certificateを別gateとし、**本packetでは一般source r体次数10を確定したとは記録しない**。NEWの定義方程式・actual target体indexの解決でもない。

## full contentを使う理由

独立に確定済みのraw determinant identityと、sに関する有理写像の次数3,5を用いる。底をK=QQ(R,Z₀)とする。二つの正規化有理関数f(s),g(s)について

\[
 [K(r,s):K(r,f)]=3,\qquad [K(r,s):K(r,g)]=5.
\]

共通部分の次数は3と5を両方割るので、K(r,f,g)=K(r,s)。従ってplane imageはbirationalにparametrizeされ、その既約implicit曲線はprimitive resultantに一乗で現れる。全てのtargetに共通する射影base-point因子は、K[r]上で(Z₁,Z₂)全係数のcontentに入る。

境界J、Fn、∞を個別に全て場合分けしなくても、この**full content**を除く方法が使える。r,R,Z₀は実source体では独立なので、非零contentはその函数体の単元である。ただしこれはgeneric function-fieldの操作であり、元の全chartでcontent零の幾何を捨てる・propernessを得るという主張ではない。

## 正確な入力と一回の計算

変数順(s,r,R,Z₀,Z₁,Z₂)のQを全6,919係数で読み、Z₁ⁱZ₂ʲの15群へ分けた。各係数群は整数3変数(r,R,Z₀)の多項式。最小の群(0,3)は124項、(5,0)は133項で、どちらもZ₀を含まず、項数が違うので比例しない。これらを最初のgcd入力とした。残りは124～1123項で、全群のr次数は19。

Qのdecompositionを別の6変数ctxへ再構成し、元の全6,919係数に一致することもnativeで確認した。係数の選択は**数値のZ₀=0などへのspecializationではない**。

FLINT3変数`fmpz_mpoly_gcd_cofactors`を順次用い、各int returnをstatusとして検査した。gcd×cofactorが両入力へ戻ることを毎回確認し、現在のgcdで次群が割り切れる場合は完全除算だけを行った。全15群の最終cofactorも保存した。

最初の二群がZ₀-freeなので、contentはfield unitを除きZ₀-freeでなければならない。最終Cも実際にZ₀次数0。Gauss lemmaにより、この3変数のprimitive gcdはQQ(R,Z₀)[r]のcontentと単元を除き一致する。R,Z₀だけの共通因子は底のfield unitなので、本methodに支障を与えない。

最後に商係数15群についてもnative gcd/cofactorを確認し、gcd=1となった。全gcd API callsは20、全statusは1。一般factorization、resultant、PRS、点・素数probe、source26は使っていない。

## 保存したcontent・全商・全体identity

|多項式|項数|変数順|多重次数|最大係数bit|
|---|---:|---|---|---:|
|C|70|(r,R,Z₀)|(9,8,0)|26|
|G|891|(s,r,R,Z₀,Z₁,Z₂)|(0,10,5,5,5,3)|57|

canonical係数SHA256は

```text
C: 7519dbba4d4dacd7ff59944a150bfe9495d8a9c5d7a275d3612518d04a4004b2
G: c31b7e3596aacc7ec66d4668670904594f4e02cac57bf77bba378e36533d3337
```

各fileのheaderを除き、整数係数と各指数をspaceで区切ったnewline列のhashである。全15群でcᵢⱼ=C·qᵢⱼ、さらに**literal 6変数の乗算Q=C·G**を検証した。全商、API status、係数hashを保存した。Gの最初の二つの係数群の商は、それぞれ7項、8項のrだけの多項式となっており、独立なQQ(R)[r] Euclid/Bézoutに利用できる。他の群のZ₀係数も全て残している。

nativeは約0.032秒、compile・入力hash・保存を含むshared約1.52秒。60秒/2GB/32MiB以内で、全保存約0.115MB。

## 最大性gate後に得られる帰結と、未解決点

nativeのfull gcd結果は、別のfield Euclid/Bézoutまたは別algorithmのcertificateで最大性を独立確認する。元のraw resultantの独立監査と3/5のbirational性は既に受理済みである。最大性が確定すれば、GはGauss lemmaによりK[r,Z₁,Z₂]で既約で、実際にrを含むためK(Z₁,Z₂)[r]でも既約となる。sはK(r,Z₁,Z₂)から回収可能だから、その時点で

\[
 [\mathbb Q(R,Z_0,r,s):\mathbb Q(R,Z_0,Z_1,Z_2)]=\deg_rG
\]

という**既知sourceのshape体の正確な次数**が得られる。

これをactual NEWの全不変函数体と既知source体の指数1や、NEWのN₁₂比を含むminimal relationと混同しない。Gはrを含むspectral base上のモデルであり、NEWの定義方程式ではない。N比によるr消去、実際のtarget index、定義式・正確な次数は依然未解決。新しい処理は自動開始していない。

## 再現

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_two_shape_content_gcd_20261104_v63.py \
  --output work/two-shape-content-replay-report \
  --work work/two-shape-content-replay-native
```

output/workは新しいdirectoryを指定する。producerは同名Python/C++。
一次保存は`reports/Q_NEW_FAMILY_TWO_SHAPE_CONTENT_GCD_20261104_v63/`。
係数15群入力、method preflight、C、全15 cofactor、G、商係数gcd1、literal全体identityと入力・証拠hashを含む。raw resultantと前の境界除算packetは変更していない。

# 斉次 module の良い還元から逆包含を証明する補題

2026-11-03 / v62。以下は exact QQ Gröbner 候補の検証に使う一般補題である。現在の非平方 J≠0 matrix に対する必要な大きい計算は未完了であり、この補題を適用済みとはしない。

## 補題

V=Z_(p)、S=V[x0,…,xr] とし、target S^n の row degree は全て0とする。斉次 column f1,…,fm が生成する graded module を I とする。G は QQ[x0,…,xr]^n の斉次 column の有限集合で、次を満たすと仮定する。

1. G は指定した同じ module monomial order に対する exact QQ Gröbner basis であり、I_QQ⊂(G) である。
2. G の各先頭係数を1に正規化した後、その全係数は V に属する。
3. 元の斉次 column を直接還元した I_Fp と (Gbar) は等しい。これは両方向の exact module membership で検査する。

このとき I_QQ=(G) である。入力 column ごとの次数は異なってよい。row shifts を用いる場合も、同じ shifts と有限な各次数の Macaulay map を保てば同じ議論が成り立つ。

## 証明

G は monic で全係数が p-integral なので、Buchberger の各 S-pair の零への除算は V 上で実行できる。従って Gbar も Gröbner basis であり、QQ と Fp で同じ leading module を持つ。特に各次数 j で

    HF_j(S_QQ^n/(G)) = HF_j(S_Fp^n/(Gbar))

が成り立つ。

一方、I の各次数 j の Macaulay map は V 上の有限行列である。その QQ 階数は Fp 階数以上なので、cokernel に対して

    HF_j(S_QQ^n/I_QQ) ≤ HF_j(S_Fp^n/I_Fp)

を得る。仮定3により右辺は HF_j(S_QQ^n/(G)) と等しい。仮定1の包含は逆向きの Hilbert-function 不等式を与える。従って全次数で等号となり、包含の graded quotient (G)/I_QQ は全ての次数で零である。これが結論である。I の V-flatness は仮定していない。

## 非斉次入力への使い方

非斉次 column それぞれを、その column の固定した全次数で新変数 w により斉次化する。target shifts=0 を保ち、上の補題を斉次 module に適用する。結論を得た後に w=1 とすれば、元の非斉次 module に戻る。localized relation zU−1 も固定次数で z*U^h−w^(1+deg U) にする。

全 column の斉次性、QQ/Fp の同じ monomial order と shifts、候補係数の p-integrality、先頭 module の保存、元の**斉次** input の mod-p GB と候補の両方向の equality を記録する。入力を dehomogenize してから mod-p equality を検査する方法では、この次数ごとの証明は使えない。

## 今回の計算の射程

既存 Cloud の Singular 4.4.1 library を再利用する薄い [実行系](../handoff/20261103_v62/SINGULAR_RUNTIME_20261103_v62.md) を追加した。`modstd.lib` の `finalTest_std_comm` は、入力の候補への包含と候補自身の Gröbner 性を検査する。研究用の有限証明書では、これに加えて上の逆包含の根拠を明示する。

J=0 の12×44と J≠0 の14×46の raw 非斉次 unit-basis 出力は、現在は候補としてのみ保存した。corrected 斉次計算は各600秒まで実行したが、前者は明示的な6GB virtual-memory cap、後者は時間制限に達し、いずれも候補 GB を返さなかった。従ってこの実行から全開の product 階数41や周囲冪零性は主張しない。J=0 の欠損有限性は、独立した四つの exact QQ minor と unit gcd による [別の証明](S61_BINARY_NONSQUARE_PRODUCT_J_ZERO_20261103_v62.md) である。

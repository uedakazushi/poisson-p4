# 新25次元被約成分のQQ/C有理性

2026-11-04。新しい既約成分 `R4_POSITIVE_WEIGHT6_25` の**被約な射影成分全体**が
QQ上で有理25次元であることを、二つの独立監査済みの証明から導く。
第一jet像B_NEW15だけの有理性を全成分へ取り違えた推論ではない。
下の二gateがFINALで通った場合に限り、専用reportがこの結論を採用する。

第一gateは `Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63_REVIEW`。
これはfull tensorの有限projective変換を用い、実S61成分Y16の関数体について

    K=QQ(B_NEW15) ⊂ QQ(Y16) ⊂ K(z)

を同じK上で構成する。sourceのkernelはfull tensorの対称性とはされていない。
raw pairのframeを合わせた後のfull tensorはなおzに依存する。
既知有限tauのfull tensor同値を使ってE5(z)上へ下降することが本質である。
次元差が一なので拡大は非定数。LürothによりQQ(Y16)=K(z0)であり、
Kの有理15次元性と合わせ、Y16はQQ上で有理16次元となる。
zそのものがLüroth生成元であること、sourceの一変数写像が次数1であることは不要。

第二gateは `S61_RATIONAL_COMPONENT_TRANSFER_20261104_v63_REVIEW`。
これは固定評価点での非退化値を六つの有理係数Cとして保持し、平方根なしの
Darboux行列と旧S61の四つの有限perspective正規化を使って

    affine-cone open ≅ U6×A4×Y16

を両向きに構成する。projective scaleを先にC03=1で固定すると

    X25 ⇢ U5×A4×Y16

が双有理となる。評価点を動かす余分な四変数を数え足していない。
連結PGL5は当該被約既約成分を保つので、全ambientスキームの積表示がこの
実成分に制限することも確認済みである。従って

    QQ(X25) ≅ QQ(Y16)(t1,...,t9) ≅ QQ(t1,...,t25).

両証明はQQ上の実有理写像と幾何的既約性を保持しており、C上でも同様に成立する。
これは既存の新成分の幾何を追加で決定したもので、成分を再度数えない。
全成分下界58、rank4下界52、総上界未定は変わらない。

有理性を主張するのはX25の被約台である。既知generic full Poisson完成環の長さ5と
CM type2は維持され、全Jacobiスキームが被約・滑らかになるとは言わない。
NEW firstjetの実定義多項式・正確な次数、全boundaryの軌道、一般betaの五点尽くし、
全Poisson成分の上界と完全分類は引き続き未解決である。

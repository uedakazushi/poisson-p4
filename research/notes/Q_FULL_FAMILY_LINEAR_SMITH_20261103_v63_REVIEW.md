# 全 Q_a 一次 Jacobi・軌道 Smith 証明書の独立監査

2026-11-03。新規実行の report と圧縮証明書は保存資料と時間以外一致した。
旧108列 helper を使わず、明記された実 Qa と canonical126 個の実二次 tensor
方向から物理的 Jacobi の全350係数を独立に再構築し、保存350×126行列と一致した。
divergence 行列の rank24、C の divergence 零と free126 行の単位行列も確認した。
全25 GL 方向は独立の単項式係数作用から再構築し、全150×25係数と一致した。
J O=0 と weighted diagonal 方向=5Qa も厳密に成立する。

全重み block の保存 S,B,T を QQ[a] 上で掛け直し、diagonal と一致した。
S,T の determinant は非零定数であり、任意の a に特殊化しても可逆である。
したがって標本点での rank 計算に依存せず、非単元 invariant は厳密に
weight6 の a−2/9、weight7 の a、weight8 の a の三つだけである。
Jacobi rank は a=0 で73、a=2/9 で74、それ以外で75となる。
軌道の17 invariant と最初の非正重み群方向7 invariant は全て単元である。
よって affine 軌道 rank17、projective 軌道 rank16 は全 a で一定である。
元の群 frame の a(9a+4) は、この内在的軌道 rank の跳びを示さない。
a=0 での Jacobi rank の実際の低下は別に残る。

全 canonical affine tangent は126−rank、projective tangent は125−rank、
群横断 tangent は109−rankである。これは一次変形の完全な結論である。
全環の平坦性、高次 ideal の枝の尽くし、Smith frame からの高次有限模型が
既に構築済みであることは主張しない。特に generic tangent rank が一定の
a=−8/9 でも、別途証明済み T40 vertical 成分による非平坦性は矛盾なく残る。

監査結果: `reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63_REVIEW.json`。

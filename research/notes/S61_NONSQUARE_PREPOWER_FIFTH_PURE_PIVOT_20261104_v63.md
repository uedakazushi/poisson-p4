# heavy の既知欠損を除いた pure 問題と第五 pivot

2026-11-04。新しい parameter 多項式の係数評価、minor/rank 探索、点・素数は使わない。
受理済み actual 6×13 と、その (0,9) 成分の非零証明を入力とする。

## 全 FULLU の heavy 被覆

`p=g1*Theta/750000`、`(g1)+(Theta)=R` は既存の全係数恒等式である。
D(p) では既存 sparse pair の determinant は `-6 B p²/5`、g1=0 では
別の既知単元 pair があり、いずれも raw-H rank は `4+rank(B7×3)` となる。
Theta 上では既存の exact cokernel が `T/(P80*P100)` で、選択 B-minor は
全域単元である。従って raw-H rank<7 の幾何学的集合は、実 Bbad13 と二本の
Theta 曲線の和集合に一致する。これは heavy Fitting ideal の scheme 等式を
主張しない。Bbad13 と両曲線上では、それぞれ別の独立 source 証明により full
cubic57 が全点で成立済みである。

従って未解決の full-cubic 欠損は H-rank7 開集合だけにある。その上で prePOWER
の heavy5 は全射であり、15 source の heavy kernel は局所自由10、pure targets
は17である。Q9 は独立受理済み a29 恒等式によって全域 unimodular であるため、
その商は局所自由9対16の pure 問題になる。ただしこの projected source の
full rank は actual full cubic57 の**十分条件**であり、失敗から実欠損を結論しない。
raw-H の Cdelta 補正も既存 source operator に残す。

## 指定された第五 pivot と同じ cokernel

先行 direct operator の residual を E6×13 とする。pivot は v=E[0,9] で、

    v = rho_u * Saux / (delta2 * rhoLCM * D4)

という全係数 source 恒等式を持つ。Saux は受理済み7022項の非零多項式である。
領域は厳密に FULLU∩D(R434*Delta5*Saux)。h22 は反転しない。p=0、Delta5=0、
R434=0、Saux=0 の側を解決したとはしない。

列9を先頭へ動かし E=[[v,r],[c,D]] と書く。parameter のみの可逆演算

    L=[[1/v,0],[-c/v,I5]], U=[[1,-r/v],[0,I12]]

により `L E U=diag(1,T)`、`T=D-c*r/v` を得る。LとUの逆、元15/22演算子への
埋込み、残る5 source の元25 sourceへの輸送を literal typed operator として保存する。
全8 Q9 tail は残る12列の先頭8列に保持される。残る source は current15
`[8,9,11,12,14]`、元25では `[18,19,21,22,24]`、列は元32
`[1,2,3,4,6,7,8,9,14,17,18,19]` である。

自由な universal E の Laurent polynomial 環で両逆積と全78成分の block 等式を
確認する。新しい actual 5×12 coefficient を展開しない。rank5 も global57 も
主張しない。指定領域で I6(E)=I5(T) となり、これは当該十分条件の未解決 locus を
正確に指定するだけである。

## 次の scalar 目標を小行列式に戻す

各行の heavy 消去を、既存 pair2 と B3 Cramer 係数、既知 rho の LCM だけで
cleared pure 行 Bhat_i にする。pair-only 行は B3 を使わない。P4 の後の numerator は

    S_i,j = D4*Bhat_i,j - Bhat_i,piv4*Adj(P4)*P4_j.

既存 Saux=S_0,15 を第五 pivot とする二重 Schur の積差には、普遍的に

    Saux*S_i,j - S_i,15*S_0,j
      = D4 * det([P4; Bhat_0; Bhat_i], columns=piv4+[15,j])

がある。従って7022項を二回掛けてから D4 で除す手順を採用する必要はない。
次の source-only 支持評価では直接6×6 bordered determinant を subset recurrence で
評価し、全60候補の支持上界と、other pair-only 行/四つの B-coupled 行を区別する。
これは係数評価ではなく source box の上界であり、最小候補でも live/height/denominator
guard を別途通すまで計算しない。全60行列を一括展開する指示でもない。

有用な追加 simplification として二つの Bhat 行の heavy clearing を wedge に戻せる。
pair-only の二行なら、pair2 の Jacobi complementary-minor identity により、この
6×6式は元の P4 四行＋heavy pair二行＋二候補行、pure六列＋heavy二列の8×8 source
minorへ戻せる。これも parameter minor の未評価 recipe であり、未知因子を単位にしない。

この packet は小さい source-bound reduction と資源前見積であり、有限 bad locus、
全 FULLU57、quartic239、generic-beta exhaustion の証明ではない。

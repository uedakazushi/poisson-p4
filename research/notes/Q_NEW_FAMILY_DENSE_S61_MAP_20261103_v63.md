# 新しい第52 rank4成分の全61座標のdense S61 map

2026-11-03。新成分R4_POSITIVE_WEIGHT6_25について、既知51族と同じ形式の
**全u1,…,u61を出力する具体的な有理arithmetic circuit**を保存した。
sourceは族n0,n1,n2、全P scalar、一般評価点四変数、Sp4のbig cell十変数の合計18変数である。
実際のS61像次元は16、primitive beta像は15、像上の一般fiberは1次元。
成分のdefining primeや一次像hypersurfaceの展開方程式はまだ求めていない。

## 全sourceを保つ有限構成

独立検証済みの全canonical150族を使い、まず一般評価点(1:base1:…:base4)を固定-one chartへ移す。
非退化な値Cの一般6係数に対する既存の有理symplectic frameを掛け、Sp4 big cell全十変数も掛ける。
次に既存のEuler quotientとperspective gaugeを使ってS61の全61 free座標を得る。
既知MPS用に導入・監査済みの汎用normalizerを再使用し、setupを再実行していない。
sourceのPoisson恒等式は全パラメータで既に証明済みであり、
各射影座標変換とS61正常化定理からgenericに元134式を満たす。
一点検算を一般恒等式の外挿には使わない。

構成は6079nodeの有限arithmetic circuitとなった。
入力名・全node・全61output・値pivot/Pfaffian node・Sp4分母・正常化行列を保存する。
分母の非消失はn=(1,1,1)、base=0、scale=1、Sp4恒等点で検査し、
値pivot=1/36、Pfaffian=−16/81を得る。

## dense性の正確な証人

同じQQ点で全18入力のforward ADを正確なFractionで計算する。
全61×18の微分rankは16で、その16次有理minorを保存する。
新しい全projective成分は既証25次元、value6＋perspective4を正常化したS61切断は16次元なので、
像は上下から16次元となり、全新S61成分にdenseである。
primitive一次部分16行のrank15のminorも別に保存するが、そのgeneric上界15は
全7変数関数体による既存独立証明を用いる。

さらに同じsourceを正常化行列で直接pushforwardしてmodel.normalizeし、
circuitの全61値が一致することと、元S61式の全チェックを確認する。
これらは実装の検算であり、generic恒等式の根拠は前節の有限構成である。

## hypersurface認識への一つの厳密な除外

この像上の記録QQ点で、intrinsic Pfaffian elliptic-quintic curve C_betaの
Fisher–Sadek不変量c4,c6,Delta,jを正確に計算した。
Deltaは非零である。従って新しい15次元一次像は、
既証のdegree24 canonical-curve discriminant hypersurfaceとは**一致しない**。
また像は既約なので、その一般点でもC_betaは滑らかである。
一次像がhypersurfaceであることを、曲線が特異であることと同一視しない。
他の既知15次元一次像との比較や、低次数Sp4不変式による認識は次の課題である。

producer `code/q_new_family_dense_s61_map_20261103_v63.py`、
receipt `reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/`。
`map.json.gz`が全source、`map-record.json`が既存catalogへ追加できる独立record、
`exact-sample.json`が全QQ値・微分・minor・曲線不変量を持つ。
既存historical known51 catalogは書き換えていない。
一般W16の五点尽くしや、全新成分のscheme的な被約性は主張しない。

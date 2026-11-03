# R19 の Log1112 包含・Poisson 有理同値の独立監査

2026-11-03。対象 producer の新規 QQ(A,B,F) replay と保存結果は時間以外一致した。
新 producer から独立に、旧 `verify_rank2_components.py` の実際の
`log_bivector` を四式 (h,x,z,K2) と参数 (5/4,1/4,3/4) に直接適用し、
全150係数が −300 倍の canonical P_Q に一致することを確認した。
旧28次元 Log1112 が full Poisson 空間の実際の成分であることは、既存の
ambient tangent/family affine29 証明書を hash 付きで再利用し、大規模計算は再実行していない。
既監査の R\\N の全標準形被覆と閉 PGL 飽和により、全閉 R は既存 Log1112 に含まれる。

K2 の Darboux 恒等式、二つの同次次数零 Casimirs、函数体の双方向の有理逆式を
独立に検証した。w を K2/(120x) に換える写像の全五変数 Jacobian は determinant1、
Euler を保ち、raw tensor の actual pushforward は厳密に P_0 に等しい。
従って P4 上でも Poisson 有理同値を与える。Q が零でなければ二次 Cremona、
Q=0 では共通因子を除いて恒等写像である。
この結論は全 R\\N に適用し、R∩N の全境界へ拡張していない。

h=1 の函数体は Casimir 値 r,t と自由変数 u,Y の純有理函数体になる。
rank2 foliation の一般葉は r,t の一般 fiber であり、有理曲面である。
rF≠0 では商 P3 の曲線の同次化は次数6の四式で、第二式 rT^6 と第四式の
T=0 における値 −FU^6/(10r) により basepoint がない。
比 z/h=U/T が逆を与えるため像の次数も6である。
元の葉閉包は vertex C∂y を持つ次数6の有理 cone となる。
曲面の正規性や R∩N 上の全葉の分類を主張しない。

監査結果: `reports/AMBIENT_R19_LOG_AND_BIRATIONAL_20261103_v63_REVIEW.json`。
新しい full Poisson 成分の追加はない。

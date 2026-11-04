# raw resultantの既知境界因子：6,919項・r次数19への厳密な除算

2026-11-04。保存済みの47,114項の整数多項式Pについて、一回の小さい計算で

\[
 \boxed{P=7378697629483820646400\,
 r^{22}Z_0(r-1)^{16}T^{15}Q},\qquad T=(R-2)r-R
\]

という**全係数の厳密な因子恒等式**を確認した。整数定数は25·2⁶⁸。Qは6,919項、r次数19である。これは保存された整数多項式に関する除算結果であり、Qのprimitiveな一般fiber方程式・一般次数19や、NEW25の定義式を主張しない。

入力Pが8×8 Sylvester determinantと一致するという大きいlibrary出力の独立identity auditは、別の未完了gateとして保持する。今回新しいresultantを計算していない。

## 除いた因子と局所化

まず全項の指数から見えるr²²Z₀を除き、整数係数のgcdを正のcontentとして除いた。次に**r−1とTだけ**で多項式除算を繰り返した。成功回数は16、15。各段階で

\[
 P_{\mathrm{old}}=f\,P_{\mathrm{new}}+\mathrm{remainder}
\]

を係数乗算で検証し、成功段階の余りが零であることを確認した。最後には全因子をまとめてQへ掛け、元のPの全47,114係数に一致することを改めて検証した。Qの整数contentは1である。

元の二shape chartではr,r−1,T,Z₀は既に単元であり、上の整数定数もQQで単元。従ってP=0とQ=0は**その局所化上で等価**で、境界を新しく無条件に捨てたわけではない。J、Fn、Hn、Nq、s、ξなど他の既知境界の影響は、この計算では調べていない。

## 最大性の証明書

最終Qをもう一度r−1、Tでそれぞれ割り、両方について非零余りを保存した。

|除数|division quotient項数|remainder項数|remainder多重次数(s,r,R,Z₀,Z₁,Z₂)|
|---|---:|---:|---|
|r−1|6478|3|(0,0,4,5,0,0)|
|T|5747|769|(0,19,0,5,5,3)|

両方のquotientとremainderを全係数で保存し、q·f+rem=Qを確認した。従って除いた16、15という冪は**この二つの指定因子に関して最大**である。一般のfactorizationを行っていない。Tは先頭係数1の既知multivariate divisorとして扱い、分母を導入していない。

全ての成功中間quotientは、最終Qに残りの因子chainを掛ければ再構成できる。各段階の項数・多重次数・canonical SHA256と正確な除算identityを`native-report.json`に保持し、同じ大きな係数表を31回重複保存する必要を避けた。

## 最終多項式Q

変数順は前と同じ(s,r,R,Z₀,Z₁,Z₂)。

\[
 \#\mathrm{terms}=6919,\quad
 \deg Q=(0,19,13,5,5,3),\quad
 \max\mathrm{coefficient\ bits}=84.
\]

canonical係数行の総bytesは218569、SHA256は

```text
9a73dabdbe8868707290c085211c32b587ddde21d0be7700fdd7595aabd4d685
```

である。`quotient.txt`先頭の`6 6919`行を除き、整数係数と六指数をspaceで区切ってnewlineを付けたhash。wrapperでも全保存行を読み直して一致した。

r次数19は**現在の未分解r射影多項式の次数**である。他の境界因子・重複・generic minimal polynomial、sの回収、N₁₂比との整合を未確認なので、一般fiber次数やactual target体indexと等置しない。NEW25の定義次数下界26も変わらない。

## 再現と保存

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_two_shape_boundary_strip_20261104_v63.py \
  --output work/two-shape-strip-replay-report \
  --work work/two-shape-strip-replay-native
```

新しいoutput/workを指定する。native一回約0.763秒、compile・hash照合を含むshared実行約2.05秒で、60秒/2GB/全artifact32MiBの範囲内。全保存約0.65MB。OpenSSLは既存libraryを係数hashだけに使い、networkを使わない。

producerは同名Python/C++、一次保存は
`reports/Q_NEW_FAMILY_TWO_SHAPE_BOUNDARY_STRIP_20261104_v63/`。
Q、二つの非整除quotient/remainder、全因子chain、全体identity、step hash、元Pのhash、入力・証拠hashを保持する。raw resultant packetは変更していない。一般因数分解、新しいresultant、別素数・点、source26、r–N消去、subresultantは未実行である。

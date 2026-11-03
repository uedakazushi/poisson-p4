# CB の一般 V における全 actual odd 計算と完成型の独立監査

2026-11-03。`AMBIENT_CB_V_FAMILY_FRAME_20261103_v63` と `AMBIENT_CB_V_ACTUAL_ODD_LINEAR_20261103_v63`、および全完成環への適用鎖を独立に監査した。係数体は QQ(V) とし、V は数値に固定していない。これは一般の V に関する定理であり、全 V の分類ではない。

## 1. 全物理 source と別実装による jet

独立コードは `code/ambient_CB_V_actual_odd_review_20261103_v63.py`。作者の tensor / Lie / jet helper を import しない。raw CB 式を自分で作り、divergence による Euler 補正と h,x の交換を行い、全150係数の F3、G3、H16を再構成する。source は必ず `F(f)−qG(f)+H` として使い、二次項 `−δq·δf_i·G_i` を保持する。

旧全 critical ideal の七つの monic 消去式に、保存された polynomial V 曲線を代入し、全 generator の消滅と V=0 で実際の P0 に一致することを確認する。全150係数の点・20 fiber 微分・19 GL 方向を再構成する。旧 fiber20＋旗9に対する source 座標変換の全29列も再構成し、行列式

\[
V^{12}(V^2-2)/3072
\]

を厳密に照合する。全 GL25列は GL19列の正確な span に入り、その階数は19である。GL19＋補 fiber10 の全 target Jacobian は rank28で、保存された全29成分の kernel を消す。

Lie 作用は `I+δE_ab` による**物理的な pushforward** と逆行列 `I−δE_ab` を δ²=0 上で使い、その δ 係数から再構成する。次に ordered exponentials の合成を、純 ε³と、29個の別々の δ 方向の全 δε係数について伝播する。この方法は作者の `L_min L_max` による Hessian の列公式とは別の計算である。初期 fiber の二次項とその全 mixed 微分も伝播する。全150×29の Hκ,j、150個ずつの W=Fκκ/2 と Z=Fκκκ/6が保存値と一致し、Hκκ=2Wも確認する。

## 2. actual base と tracefree odd 係数

28個の独立な実 target 関数 Y と、その補空間 Jacobian `D` を使い、選択 Y-minor を正確に反転する。`rem(M)=M−D(D_Y)⁻¹M_Y` は28接方向の除去である。κと28補方向の mixed Hessian を H と書けば、Yを固定する二次の補正は `−(D_Y)⁻¹W_Y` であり、純三次は

\[
Z_c=\operatorname{rem}(Z-H(D_Y)^{-1}W_Y)
\]

となる。raw な純 κ 三次の Z だけを使わない。branch の二次係数 w≠0を選び、`W_c=rem(W)` による branch 成分を全 output から引く。これは同じ Y、同じ branch 値を持つ二つの source 点の差を取る際の trace 補正である。二次 branch の linear Y 部分も残して actual R29座標へ戻す。

こうして得た **全150×29 actual odd firstjet** は保存行列と厳密に一致する。九次 minor は

\[
\frac{157837977V^{28}(V^2-2)^9}
 {125440000000000000000000000000000000000}
\]

で非零である。選択9列による全列の factorization も確認し、上界・下界の両方から rank9を得る。全 source critical cotangent の rankも9となる。V=1での同じ QQ minor の32003還元は非零だが、これは rational function の検算であり、V=1 が後述の全 proper finite 開に入るという主張には使わない。

GL19方向に `dP/dV` を加えた実 target tangent の階数は20である。従って19次元の一つの GL 軌道で modulus を消していない。GL `E_yy` の保存列は全 V で P 自身に等しく、scalar方向も明示的に確認する。

## 3. 一般 V の全 proper fiber

全 proper CB 旗模型、P0での全旗 fiber一台長2、全 source critical ideal については、既受理の独立証明書 `AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW`、`AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW`、`AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW` を明示的に使う。これらの既存 GB は再実行しない。

proper morphism の全 P0 fiber が有限なので、P0を含む target 開で写像は有限になる。有限 algebra の fiber 次元の上半連続性により、さらに全 fiber長≤2の開を取れる。この開の V 曲線への逆像は V=0を含む非空開であり、generic Vを含む。generic source 微分の corank1により選択点の local fiber長は2以上であるから、全 fiber長≤2と合わせ、全 fiberは一台長2に確定する。ここで V=2等の任意の標本がその開に入るとは仮定しない。

28実 target の可逆微分と非零 branch 二次係数により、29実 target 関数の parameter fiber は正確に長2となる。完備 Nakayama と正則 source の CM 性から、source B は正則29次元の actual底 R 上で free rank2である。trace を引き、`B=R[t]/(t²−ν)` と書ける。全 finite source fiber を含むため、proper imageとして定義された actual CB の完成像も、この B の中の `A=R⊕I t` である。

QQ(V) による計算は、基底拡大後の Vを座標とする点で全29 source方向を残して行う。元の曲線の generic prime における28次元局所環と取り違えない。有限個の非零式を spread out することで、複素数上の V 曲線の非空開に結論が適用される。

## 4. I=J と normalization の順序

全 source critical quotient は `B/(t,I)=R/(ν,I)` である。旧の正確な critical ideal は、V≠0の枝で滑らかな20次元となる。新たな全 critical cotangent rank9と GL19＋modulus の実20方向もこれに一致する。従って `J=(I,ν)` は正則底 R29の正則九生成 ideal である。

actual I の全 firstjet rank9は J の九つの cotangent class を張る。従って `J=I+m_RJ`、Nakayamaにより `I=J` となる。B が正則であるため ν は R の正則 parameter の一つである。よって R の座標を選び、`I=(s1,...,s8,ν)` とできる。

ここまでで I≠0も分かる。非零 c∈I を取れば `t=(ct)/c∈Frac(A)` なので Frac(A)=Frac(B)となる。**この段階で初めて** finite normal B を actual normalization と同定する。全 global proper 模型の未確認な birationality を前提にしていない。

## 5. 受理する一般 V の完成型と範囲

coneの全完成像は

\[
A=\mathbb C[[g_{20},s_8,t^2,s_1t,\ldots,s_8t,t^3]]
\subset B=\mathbb C[[g_{20},s_8,t]].
\]

これはCB自身の計算から得られる。R-moduleとして `A=R⊕I t`、I は正則列9本の idealなので、coneの dim29、edim38、depth21を得る。odd生成元の平方は m_R³に入り、m_RA は m_A の reductionになるため、R上の階数2から multiplicity2となる。scalar方向を滑らかな因子から一つ除く projective版は dim28、edim37、depth20、multiplicity2である。

正則列の Koszul 関係36本と pair-product関係45本が全 kernelを与える。pair productで次数を1以下に戻し、残る linear relationが正則列の第一Koszul syzygyで尽くされるためで、新しいGBは必要ない。conductorは B内 `(t²,s8)B`、A内 `(t²,s8,t³,s8t)A`、`B/A= C[[g20]]t` であり、この formal germの singular / nonnormal / non-CM supportは滑らかな conductor supportである。`t²,t³∈A` により、この有限正規化は elementary subintegral、従って universal homeomorphismである。

全 V、V=0、frameの V²=2、B18原点の全 odd ideal、actual intersection Dの型を結論しない。frame の極は座標証明の除外であり、その点を幾何学的例外と認定しない。

再現:

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/ambient_CB_V_actual_odd_review_20261103_v63.py \
  --output work/cb-V-odd-independent-review
```

保存結果は `reports/AMBIENT_CB_V_ACTUAL_ODD_LINEAR_20261103_v63_REVIEW/report.json`。既存の長い幾何証明は上の三つの独立 receiptを依存として保持し、新しい全 source / jet / rank の検算と区別する。

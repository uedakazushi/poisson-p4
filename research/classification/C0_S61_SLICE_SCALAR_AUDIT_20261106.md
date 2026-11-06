# S61の全rank4被覆・一般切片・scalar判定の依存監査

2026-11-06 / v64-classification-plan。担当C0。対象は複素数体上の射影Poissonスキームであり、PGL5商を取る前の被約成分分類である。本ノートは分類の完了を主張しない。

## 確認した矢印と証拠の強さ

| 矢印 | 今回の確認 | 原証拠と残る条件 |
|---|---|---|
| rank4の全点→群移動したS61 | 下記の幾何学的証明を読解・再構成 | `S61_MODEL_AUDIT_20261030.md`、v02本文§3。55消去・20Bianchiの全係数再生成は今回行っていない |
| 保存S61一次16→primitive beta16 | 元150係数から正確な線形写像を独立再構成。det=36/625、両側逆積を確認 | 本ノートの明示16式。二次45座標のgauge変更も必要である点を保持 |
| 標準六切片の優越性 | 元24テンソルの有限pullbackを微分し、全160作用列とdet=8428695240を再現 | `primitive_curve_dominant_slice_20261101.py`。Fisher不変量の計算は不要で再実行していない |
| 既存疎五切片の利用 | Sp4×Gm×源5の優越性を独立QQ Jacobianで確認。全134式の重みとlambda除去恒等式も確認 | 疎六切片と`S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63`。Sp4単独の五切片を優越としない |
| 全134式→normal42/J0 | finite5・commutator・global quinticの証明を読解。環の有限性・閉埋込みの論理を点検 | 既存独立exact証明書に依存。今回巨大126逆行列や全normal係数は再計算していない |
| Bの空性⇔delta_E所属 | regular性の証明と所属判定の両方向を点検 | 下記の仮定表。delta_E所属そのものは未証明 |

C0の用途はGの入力と論理を確定することであり、証拠の所在を確認したことを全計算の再実行と同一視しない。

## 1. rank4の点を取りこぼさない理由

Pを一般階数4の任意のPoisson構造とする。rank4の点をp=[e0]へ移し、接空間の線形座標を変えると、その値をJ=∂12+∂34にできる。これは点ごとの標準化であり、一般betaへの仮定を使わない。非退化評価の開Uの評価写像はGL4/CSp4上のassociated bundleとなる。CSp4は連結なので、各fiberの有限個の既約成分を置換しない。PGL5も連結なので全スキームの各既約成分を保つ。従ってS61の全点を既知族の和で覆えば、群作用によってrank4の全点を覆える。

定数項Jを固定したfiberには、affine次数1、2、3にそれぞれ24、60、36係数がある。四つのperspective方向を除く操作は接空間だけの商でなく有限の代数的正規化である。aを四成分、yをaffine座標、f=aᵗyとすると一次部分は

    L' = L - 2 f J + y∧(Ja),
    tau(L') = tau(L) - 5 aᵗy,   tau(L)=L12+L34.

したがってa=grad(tau(L))/5でprimitive gaugeを全域で一意に選べる。有限変換は元の二次・三次部分も変える。保存gaugeからprimitive一次座標を得た後も、保存二次座標45個を無断でprimitive二次座標へ読み替えてはならない。

保存模型の55消去は非零有理定数の三角消去であり、20Bianchiはideal内の恒等式である、という原exact証明を用いる。これらが全恒等式として成立する条件で、形式近傍に限定されずUは基底5方向・perspective4方向とS61で記述される。成分の次元差は9となる。旧本文の既知40成分という歴史的候補数は現行の58下界や尽くしと混同しない。

## 2. 元S61からprimitive betaへの明示的な橋

uは元の`parameters[k]=u(k+1)`。重み1の番号は

    4,5,12,13,14,25,26,27,34,35,36,37,48,49,50,51.

beta=-J L' Jを、`weight1_structure.constraints()`のkernelの標準SymPy nullspace基底B24×16でbeta=B bと表す。独立再構成した橋は次の16式である。

    b0=u34;                 b1=u36;
    b2=-u25;                b3=u37+4u48/5;
    b4=-u27;                b5=u35+u50/5;
    b6=-u12;                b7=-u13;
    b8=-u26+u51/5;          b9=u4;
    b10=-u14+4u49/5;        b11=u5;
    b12=3u48/5;             b13=3u49/5;
    b14=2u50/5;             b15=2u51/5.

上記順序でdet=36/625。この行列と両側逆行列の積、betaの閉性・原始性の全8条件をQQで確認した。再構成はadapterの`primitive_beta` helperを使わず、保存150係数の自由列と55復元式の一次単項式から直接行った。Euler dehomogenizationは

    Paff_ij = Pij|x0=1 - y_i P0j|x0=1 + y_j P0i|x0=1

で行い、一次部分を取り出した。これは保存S61座標との橋であり、任意の一次betaが一意にPoissonへ積分できるという主張ではない。

標準六切片はv=(1,...,16)と標準基底e11,...,e15が張る空間である。十個のSp4方向と六切片方向の16角接行列式は8428695240。作用列はA=I+t aについて`Aᵗ beta(Ay) A`のt一次項を直接取り、十生成元×十六基底の全160列を照合した。

## 3. 既存generic5回路を使える正確な意味

疎六切片は`section6.json`の

    beta6(lambda,H,L,s0,q0,epsilon)
      = lambda*b(H,L,s0,q0) + epsilon*e4

である。lambda=H=L=1、s0=q0=epsilon=0において、元24テンソルから再構成した十群方向と六源方向の16角Jacobianは

    137438953472/996451875 = 2^33*16/996451875 != 0.

有限個のQQ値による**非零Jacobianの証拠**なので、標数零の優越性を証明できる。零の標本から大域零を推論しているわけではない。保存された全記号式は

    detJac|epsilon=0 = (2^33/996451875)*lambda^14*H^3*L^2*Gamma,
    Gamma=1+5q0+8s0+15H^2+115H^2q0.

今回独立計算した非零値はこの式とも一致した。全記号Bareiss式は既存の別実装reviewに依存し、今回は繰り返していない。

lambdaを単元としepsilon=lambda*t0とすれば

    beta6(lambda,H,L,s0,q0,lambda*t0)
      = lambda*beta5(H,L,s0,q0,t0)

は全16係数の恒等式である。元134式は重み4が70本・重み5が64本で、weight1/weight2に対する物理Gm作用は(u1,u2)→(lambda*u1,lambda^2*u2)。全134式の各単項式で重みを再確認した。ゆえにgeneric5模型はSp4**とGm**を合わせて一般betaを支配する。同じGm作用で45未知数全体をlambda^2倍し、全元134式と逆写像を運ばなければならない。Sp4×源5だけでは次元15であり優越しない。

群と源の関数体を含む忠実な拡大へ全問題を基底変換した後、同変性で群変数を除ける。特定beta点への代入だけをこの体拡大の代わりにはできない。beta-only selectorの非零は切片への制限で別途確認する必要があり、未知fiber座標の逆元を追加した場合はその補集合を保存する。

## 4. scalar経路の正確な仮定と結論

K=C(W16)、R=K[n0,...,n41]、A=全134式の一般beta代数とする。以下の既存証明を依存として使う。

1. 正規座標変換は全45未知数の可逆変換であり、全134式をF5+G129に可逆行変換する。既知三平面n=0の交叉は有限エタール次数5のE5である。
2. global splitterによりA=E5×B、(n)A=(1-e)A。126個全てのcross積は`z_k*n_j=P_kj(n)`、P(0)=0、次数≤2へ解け、他の三行Hも保持する。選択126minorの分母はbetaだけに依存する非零多項式である。
3. AはR上5生成以下の有限加群。BのR射影は閉埋込みでoriginと交わらず、R→Bが全射である。R上有限であることはK上有限を意味しない。
4. normal-only idealはJ0=m∩I_B=ker(R→A)。2838生成元の次数≤4。原点でのcommutator Jacobianが可逆なので付加originは被約である。
5. global quintic `C=R[h]/(g(h)-H(n))`、`A=C/(r0,...,r133)` は元134式と同型。placeholderのS1,S2,Hは自由変数ではなくexact normal circuitを代入する。

未知幾何点pではR→Bの全射性からA⊗κ(p)=κ(p)。従ってmonic quinticの上で全remainderはideal(h-h0)を生成する。r_i=(h-h0)s_i、D=κ(p)[h]/(g/(h-h0))とすれば(s_i)D=D。transcendental aを使うmoment r(a)=Σa^i r_iのbottom4 minorはNorm_D(Σa^i s_i)となり、各pで非零多項式である。

したがってdelta_E(a)のa係数idealはB上で単位idealであり、delta_EはB[a]・B(a)上の非零因子となる。原点ではdelta_E=0。よって

    B=0  ⇔  delta_E belongs J0 extended to K(a)[n].

十分性はmembershipからdelta_EがB(a)を零化し、regular性からB(a)=0、忠実平坦な体拡大からB=0と戻る。必要性はB=0ならJ0=mであることによる。これは根基所属・unit・有限体空性を混同しないexact criterionである。

Euler版targetの次数28とcofactor次数24以上の障害は既存の独立証明に依存する。今回actual delta_E所属、未知残余空性、一般五点尽くしは証明していない。Gには直接根基所属またはこのscalar証明書という追加の数学が必要である。

## 5. 再現と入力束縛

今回の小計算はC0のjob契約を`classification.py gate`で確認してから実行した。標準橋0.85秒、疎五切片0.35秒、各55秒・1.5GB制限。Fisher不変量、全126逆行列、B18/T2などのsecondary値は計算していない。

今回の小計算は次の一つの再利用可能なrunnerと永続証明書に統合した。

```sh
.venv/bin/python code/classification_slice_bridge.py --output work/new-slice-bridge.json
```

保存証明書は`reports/CLASSIFICATION_SLICE_BRIDGE_20261106.json`で、明示16行列、24×16埋込み、標準16角接行列、疎六切片の16角接行列、元入力hash、runner hashを含む。両側逆積・全160有限pullback微分列・元134重みの一致をrunner内で確認する。このrunnerは全55消去や134式自体の再生成ではない。

数学的な再現手順としては、2節の16式からCを構成し、`weight1_structure.constraints()`のnullspace Bと`infinitesimal_action`で接行列を作ればよい。独立性のため作用列を有限pullback微分で照合した。疎六切片の非零接行列は、`section6.json`のbetaに十作用と六源微分を並べ、3節の有理点を代入してQQ行列式を取る。復元55式と元134式の入力hashは下表の通り。

| 入力 | SHA256 |
|---|---|
| `code/s61_model.py` | `347fd285c5342d6e6444bf6d049a49342650a22b29857f87bbd5f4b689400fc1` |
| `code/weight1_structure.py` | `204578a6299862e22494586e94cf7e69633ec2f9812a82e595457a6e75d8c22b` |
| `legacy/v02/data/symplectic_slice_certificate.json` | `f9b1841de278d3fb663caaf6d3d64196fde60f8bc1831bf2bf1baa28812013a5` |
| `legacy/v02/data/symplectic_slice_equations.json` | `f66156dbbeb1e6c56bef852dd99fae0d534a2b631363d140b372a1773deeeabb` |
| `legacy/v02/data/symplectic_slice_reduced_equations.json` | `e8e38d8c63f2e6c4db2a905420b6bdebd283f5acdc5814a04bf346283f4f3df5` |

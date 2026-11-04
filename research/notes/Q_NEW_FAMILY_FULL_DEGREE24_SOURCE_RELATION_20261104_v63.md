# 全24次不変式のsource制限は単射：NEW25定義次数は26以上

2026-11-04。**dim Inv₂₄=67の全空間についてNEW25への制限がQQ上で単射であることを証明した。従って24次の定義方程式は存在せず、既知の下界24と奇数次の消滅を合わせ、NEW25の定義方程式の次数は26以上となる。** 方程式そのものと正確な次数は未解決である。

これは先の41部分空間の検査とは異なる。独立に受理済みの全67本の普遍基底と、独立に全係数が検証されたsource mapを使い、一度だけ固定素数32003の**多項式係数rank**を計算した。sourceの新しい点、点評価行列、別の素数、高次数への探索は行っていない。

## 実行前に固定した二つのgate

1. `Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63_REVIEW`：同じ形式ambient平面上の係数minorとcharacter上界により、選択した67本が普遍的QQ空間Inv₂₄の基底であることを独立確認した。旧70候補から0-basedの62,68,69を除いた列である。
2. `Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW`：actual1032-node source、元の16×16 quartic Hessian、全4本21834係数の直接計算から、新しいsource係数を独立確認した。

両receiptの`all_passed=true`と全入力hashを確認してからnativeを一回だけ開始した。既存のF₄,pᵢ,N₁₂の係数と、四つの新しい係数を使う。物理的に保存した11 mapsのうちp₀はF₄²=400p₀+1600p₄から従属なので、この検査の生成元は10個である。

\[
 (F_4,p_1,p_2,p_3,p_4,
 N_{12},\operatorname{tr}(C^3DCD),\operatorname{tr}(C^3D^3),
 \operatorname{tr}(\theta^5H_{16}),\operatorname{tr}(\theta^7H_{16})).
\]

β次数は(4,8,8,8,8,12,12,12,12,16)。全70個のmonomialと選択67個の**実際の10生成元における指数**を`generators.json`へ固定した。70列の普遍関係3本をNEW25の関係と誤認する余地はない。

## 係数証明書

sourceの変数はa=u/X、Y=X²、W、Bであり、H=147456Y−(3072W−1)²。各生成元は整数numeratorを明示的な整数・B冪・H冪で割って保存されている。24次の選択67 monomialは全て、B¹²H²⁴を共通分母として、さらに列ごとの正の整数dⱼを追跡すれば

\[
 P_j=d_j B^{12}H^{24}M_j\big|_{\mathrm{source}}
 \in\mathbb Z[a,Y,W,B]
\]

となる。全dⱼは32003で単元。

この67整数多項式を全係数で比較した。選択67個のmonomial係数行に対する67×67整数行列の行列式は

\[
 \boxed{9857\pmod{32003}\ne0}.
\]

従ってPⱼはQQ上で線形独立であり、非零の共通分母と列定数を戻した67個のsource有理関数もQQ上で独立である。これは有限体上でrankがたまたま上がるという推測ではなく、**整数係数行列の非零minorによる特性零の厳密な証明**である。QQ kernelを推測したり、QQへのmodular relation liftを使ったりしていない。

nativeは全67列を構成し、全50,277個の係数monomialの和集合を記録した。選択係数行と全67×67整数行列は`native-certificate.json.gz`に保存した。選択行は[0,0,0,0]から始まり、最後は[6,3,0,3]。各成分は(a,Y,W,B)の指数である。Pythonも保存された全整数行列からmod32003のminor9857を独立に再計算した。全67列の個別項数はfrontierとcertificateに残した。

## 全Inv₂₄の排除と次数下界

選択67本は既に全Inv₂₄のQQ基底である。従ってこの係数minorは、**全24次不変式空間のsourceへの制限が単射**であることを示す。source sectionはactual NEW25の構成から得られるため、NEW25全体で消える不変式はこのsourceでも消えなければならない。従ってNEW25で消える非零24次不変式は存在しない。

primitive first-jetの像NEW25は受理済みのSp₄安定な既約homogeneous hypersurfaceである。Sp₄には非自明なcharacterがないので、その既約定義多項式は不変式として取れる。中心元−IはW₁₆に−1として作用するから、不変式の奇数次成分は零。従って25次も存在しない。既に確定した下界24と合わせ、定義次数は**26以上**である。

sourceのδ-constant/t=2 sectionに関する正のrelation検出criterionも独立に整理済みだが、今回のnegative結果はsource restrictionの単射から直接従う。actual NEW25 orbit invariant fieldとsource定数体の有限indexが1であるという仮定は使わない。

## 実行範囲と再現

native一回、GF32003の係数rank一回だけ。実測native約29.77秒、wrapper約31.50秒。最大列・中間support49,015項、polynomial products187回。hard capは60秒・2GB、列100k項、中間150k項、certificate32MiB。全て上限内で完了した。

native入力の圧縮約1.11MB、選択整数行列を含むcertificate圧縮約48KB。完成列ごとの項数と実行段階を逐次`frontier.json`へ保存し、全67列完了とrank67を記録した。

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_full_degree24_source_relation_20261104_v63.py \
  --output work/full24-source-replay-report --work work/full24-source-replay-native
```

output/workは新しいdirectoryを指定する。一次保存は
`reports/Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63/`。
`generators.json`、全入力、全選択係数行列、native log、逐次frontier、入力・gate・証拠hashを含む`report.json`を保持する。producerは同名Python/C++。

`exact_all_integer_polynomial_coefficients_identity_verified=false`はrank67のためrelation候補・QQ kernel経路を実行しなかったことを表す。計算失敗ではなく、`all_Inv24_QQ_restrictions_injective=true`が結論である。

定義方程式、正確な次数、multisectionの次数・genus、actual quotient-fieldの有限indexは依然未解決。26次以上の不変式検査、別の素数・点、自動再試行は開始していない。

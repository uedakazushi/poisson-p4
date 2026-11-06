# 次の課題 — 全被約成分の尽くし

2026-11-07 / v64-classification-first。classification/PLAN.jsonが現行キュー。

```sh
.venv/bin/python scripts/classification.py check
.venv/bin/python scripts/classification.py next
```

GとE_quadric3を並行する。Gの新しい入口はclassification/DOMINANT_COMPONENT_GEOMETRY.md。
未知beta優越成分は射影次元25〜27。正規被約退化因子のT39以外の分岐は、
可解unimodular零点を持ち、その退化因子重複度は3以上となる。
全高次延長を含む非可解一次部は形式線形化により排除済み。
次のG-normalは、この正規分岐のS61からprimitive betaへの像が次元15以下となること。
像上界は未証明で、非被約・被約非正規退化因子も別に残す。

全134式のnormal rank4枝とF5、有限分離忠実chart、直接splitterは利用可能だが、
それらだけからB=0とはしない。独立reviewの五点+孤立第六点反例が限界を示す。
五回未完了の直接GB、低次trace探索、旧compact方向graphの再包装を反復しない。

E_quadricは完了。実退化因子H+2Q、rank Q>=4の全支持を既知閉包で覆った。
接触配置の全点はT18へ入るisotropic源と、実変形でMPSへ入るanisotropic一軌道に分かれる。
rank3横断配置も、全20次元の閉逆形式を保った変形でbeta零飽和閉包へ入る。
E_quadric3では残るrank3二配置を処理する：
Q=u²+v²+z²の全19次元閉形式と、Q=w+u²の全25次元閉形式。
後者の六追加クラスβ=u dr−r du, r∈Sym²(v,z,w)を保持し、Pfaffian=kQ²,k≠0を全て課す。
新計算はactive E_quadric3に対する有限job契約を作ってから実行する。

E_zeroとE_poly、A_rank2は完了。現在の既知58成分に対するAの依存監査も完了。
旧証明書を入力にする範囲と今回の独立再計算を区別し、新候補だけを必要に応じて追加監査する。
全成分下界58、rank4下界52は維持し、上界・一般五点・全分類を未証明と表示する。
各checkpointは完了宣言ではない。mainへの保存後も、許可の再確認なく次の有効な数学へ進む。

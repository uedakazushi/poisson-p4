# CA self-pairの実CB lift：全非線形source微分の独立監査

2026-11-03。指定した一つのCB liftのimmersionを検証し、whole CB fiberやtargetの滑らかさとは区別する。

著者の微分ルーチンを使用せず、実5×5テンソル

    canonical(rawCB(f,mu*f,H))

を19fiber変数とratio muについて先に組み立て、全20方向で微分した。
fの三係数とmuの積を保持し、全三つのmixed二階微分も全150行で照合した。
actual f=3a x²、mu=0でのratio列はcanonical(rawCB(0,3a x²,0))であり、自由な三つのg列で置き換えない。
既存の全CB flag parabolicの補空間9方向を、実5×5 Lie式で独立に微分した。

全150×29行列が保存証拠と一致し、QQ(a)でrank29。
fixed19を除いた150×10行列と19次inverseも全係数で一致した。
選択10次minorは正確に53747712 a^10/15625であり、全a≠0でこの指定liftはimmersionである。
fresh著者certificateのbyte一致と現在の依存hashも検証した。

これだけではwhole CB proper fiberの点数・length、個別CBのnormality、smoothness、completion型は分からない。
CAの二旗をCBへ転用せず、Nの一旗の結論とも区別する。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CA_self_pair_CB_rank_review_20261103_v63.py --output reports/AMBIENT_CA_SELF_PAIR_CB_SOURCE_RANK_20261103_v63_REVIEW.json --replay work/v62-review/CA-self-pair-CB-source-rank-replay

GBなし。公開著者証拠は変更しない。

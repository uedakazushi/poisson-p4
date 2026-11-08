# HA_1・t=2/15：全軌道閉包は滑らかな23成分候補

2026-11-08、`/root/quartic_review`、PRIMARY E。
whole-source純証明 `PLANE_MOVING_SIMPLE_TWO_WHOLE_ORBITS.md` §4と、
一つのcanonical fullambient算術を接続する。**新成分の台帳採録は別作者の全算術／全prior128監査後**。
現時点でcomponent countを変更しない。

## 完全なactual particular

座標(a,b,c,z,w)、Π=Q+∂c∧Wで

    Qab=w², Qaw=−9z²/14, Qbz=z², Qbw=2zw,
    Qaz=Qzw=0,
    W=(a²−27bw/7,−2ab,−2az−3w²/5,2aw/5+27bz/5).

W0だけをactual particularと読み替えず、旧affine系が強制した
−27bw/7及び−3w²/5を両方含む。
整数scale70で元全十括弧を入力し、全350 original projective Jacobiを零と確認した。
全Euler-Pfaffianはscale70の平方を除き

    D=−(3/5)w(w²+3az)²−(27/14)b w²z²−(27/5)abz³−(9/14)c z⁴.

## 全周囲350×150と実群作用

既存の別作者Fraction疎多項式engine
`classification_degenerate_quadric_q4_nil31_I_review.py`を不変のまま再利用した。
新wrapper又はshared code appendは行わず、公開input JSONにexact実行programを保存した。
元の全10pair×15monomialの150方向から350 derivativeを再構成し、
全25Euler線形列、全25GL5作用列、tensor scalarも定義から構成した。
rankはそれぞれ

    rank D350x150=102,
    rank Euler150x25=24,
    rank(Euler,GL5,scalar)150x51=48.

各行列の非零最大minor、その正確なdet、完全kernel、全free座標identityを保存し、
元の全rowへのkernel代入零でrank上界も検証した。
全350×51についてD·M=0、actual quintic恒等式も確認した。
有限体、数値rank、部分matrix上界又はnormal jetは用いていない。

raw Jacobi schemeの接次元は150−102=48。
実Euler additive／GL5／非零scalar作用の48次元軌道が通るので、局所次元も48。
したがってこの点はregularで、唯一のfull component上にある。
全24Euler directionsと非零scalar一方向を除いた射影局所次元・actual PGL5軌道次元はともに23。
actual軌道は局所閉かつこの成分で稠密openなので、その閉包が23次元全成分である。
whole-source純証明のglobal inverseを合わせれば、元HA残六源はその同じ軌道に入る。
この数学的結論の**独立算術review**は本稿作成時には未完了。

## 実行identityと一回目のmetadata欠陥

最初の実行では上記三rankと全kernel検査を終えた後、JSON構成時に
job辞書`j`をPfaffian loop indexが上書きしてTypeErrorとなり、reportは生成されなかった。
V1 input／jobの全文は不変に保持し、実stdout／stderr／exit1を
`CLASSIFICATION_PLANE_MOVING_SIMPLE_HA_AMBIENT_V1_FAILURE_20261108.json`へ保存した。
V1にruntime又は保存されなかった行列のreceiptを捏造しない。

**NEW V2**はjob辞書名をcontractへ直し、input／job path identityを変えただけで、
同じcanonical点・同じ全行列／kernel演算を再実行した。
V2 gate PASS、2.418325秒、最大18,476KiB。
保存reportは `CLASSIFICATION_PLANE_MOVING_SIMPLE_HA_CANONICAL_AMBIENT_20261108.json`。
全行列のsparse entry、全kernel、exact executed program及び全入力hashを含む。
追加点・profile・高次モデルは実行していない。

## 残る相違監査

Dはprimitive c-linearなのでC上絶対既約。
唯一の被約特異平面はz=w=0、一般normal cubicは
−27a z²(aw+bz)/5でroot multiplicity21である。
五偏微分は独立（順にcz³,w⁴,bz³,az³,z⁴の係数）なのでessential5。
この二点は既essential≤3成分、actualD可約閉族、二／三特異平面を持つ既99／型7成分、
triple normal leadingを持つ新weak-wall二成分との相違監査に使える。
ただし固定prior128の全rowとそれぞれの閉条件／一般元証明の照合は次の別義務であり、
本noteだけで新129採録済みとは呼ばない。

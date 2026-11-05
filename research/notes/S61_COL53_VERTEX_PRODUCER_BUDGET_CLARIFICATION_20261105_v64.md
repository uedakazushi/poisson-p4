# col53 producer の decoded VALUE 所有量：メタデータ補足

これは再算術を伴わない履歴・所有量の補足である。
報告 SHA256 は `b80ff917eac6eb9be2e5ea0673c83903026b126610900bfc678279e6a6207952`。
元の SOURCE gate、十係数 VALUE、独立 peer の全バイト列は変更していない。

元 gate の84,916は、JSON前の tuple頂点座標を generic counter が数えなかった場合の
20,856＋24,060＋40,000であった。JSON後は40,100個の指数ラベルが list になるため、
decoded graph は同じ conservative counter で60,956 slotsである。
従って84,916を decoded VALUE の literalな全体計数とは読まず、この補足で範囲を訂正する。
独立 reviewer の118,477は別レイアウトの値であり、元 producer の訂正値ではない。

producerについては、decoded graph60,956、最大1,203×10個の QQ 値の分子分母24,060、
既存一時領域40,000を保持すると、全体の正値上界は **125,016** となる。
元の hard150,000以内であり、capraiseはない。

一時領域は最大入力1,762項から因果的に確認した。
元 primitive polynomial は一項当たり numeric string三＋五指数、係数辞書は一項当たり
五指数 key＋整数一で数える。旧 p/coeff を保持したまま新 decode が完成する最大重なりは
38,780、旧 coeff と新辞書 comprehension の重なりは35,248以下である。
新 binding document 全体の最大は6,517で、旧 p/coeff と合わせても前者より小さい。
残り1,220を root/order/result・loop/parse/scalarの一時メタデータに割り当て、
既存40,000の予備内に収まる。

この補足では既存 graph と literal document のサイズだけを読んだ。
頂点係数、全成分、面、parameter point、prime、rank の算術を再実行していない。
十係数の零とその有限 Laurent 順序の証明は、元 VALUE と独立 VALUE peer のままである。
再現コードは `code/s61_col53_vertex_producer_budget_clarification_20261105_v64.py`。

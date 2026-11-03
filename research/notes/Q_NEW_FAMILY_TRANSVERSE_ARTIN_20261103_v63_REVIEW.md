# 新25次元成分の閉transverse fiber長さ5の独立監査

2026-11-03。対象は実族のQQ点n=(1,1,1)で26個の線形座標を固定した**閉fiber**である。
全Poisson完備局所環の一般多重度5・平坦性・26形式変数との直積を、この一点からは結論しない。

## 全126方向・実26切断と98形式解

reviewerは保存実150族から全GL25+scalar方向と三族方向を独立に作る。
rankは群だけ24、族も含め26。26選択座標への非零minorを確認する。
この26方向は全て群方向ではなく、実族の二つの追加方向も必要である。
全canonical126方向の残る100方向で、物理的定数二次積表から全350微分を独立に作り直す。
保存350×100と一致し、与えられた98×98逆行列の厳密積を確認する。
二つの独立kernel方向からrank98の上界も確認し、二変数の実形式presentationを得る。
大きな98逆行列そのものを再構成する計算は繰り返さない。

## 全QQ二変数係数・全350式

reviewerの別C++/FLINT engineで、P+R1+R2+R3+R4の物理的Jacobiを斉次に畳み込む。
全350式、二変数の全QQ係数、次数0–4とdivergence25×5を保存残余と比較する。
方向標本や有限体の検査ではない。
二次一行A2は非退化binary quadraticで、そのdiscriminantも報告の有理数と一致する。

残る349生成元をorder≥3へ変える定数row操作を確認し、
S1*A2と全349三次を並べた351×4行列を独立に再構成する。
rank4なので全三次monomialが真のinitial idealに入る。
保存347行kernelについて零積と独立rank347を確認し、全四次carryも再計算する。
四次rank5は整合するが、Artin性は三次で既に決まる。

実閉fiberの局所環でh3=0ならm³=m⁴=m*m³であり、Nakayamaによりm³=0。
次数0,1,2は厳密に1,2,2なので全長5である。
全高次項は消えるためidealは(A2,m³)となる。
非零discriminantのbinary quadraticはC上でuvへ移るので、
実閉fiberはC[u,v]/(uv,u³,v³)と同型である。
独自に作った小さなQQ M2入力でもdim0、degree5、m³包含を確認する。

## 完備局所環の有限性と、一般構造への限界

選択26 source座標から実族の形式mapを取り、26切断座標への射影は微分同型である。
形式逆写像により、その実Poisson branchは26変数形式環Rのgraphになる。
全完備局所環AにはR→A→Rのsectionがあり、合成は恒等である。
閉fiber A/m_R Aの長さ5から、形式有限性と完備NakayamaによってAはR上有限、5元以下で生成される。
これは「26群変数の形式直積」の証明ではない。

generic Jacobian rank98が別に証明されても、単独では一般長さ5の証明にはならない。
一般二次obstruction rankの上界等による下界と、この有限生成上界との比較が必要である。
その追加の一般関数体結果はこの点packetから推論していない。

review code `code/q_new_family_transverse_artin_review_20261103_v63.py`、
receipt `reports/Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63_REVIEW/`。

# 普遍G6の実first-w-colonとaxis idempotentの合成

2026-11-03。既存の普遍G6証明書をgeneric61変数のまま再生し、raw homogeneous idealの実際のpartial colonとして位置づける。新axis splitterと合成し、未知アフィンの全schemeを保つ中間idealを指定する。全w飽和や尽くしは主張しない。

## 1. 実次数と全generic identity

元70本weight4式をF_i、保存125項linear syzygyのmultiplierをK_i(q)とする。元61変数のQQ polynomial dictionaryで

    G6(beta,q)=sum_i K_i(q)*F_i(beta,q)

を再構築し、保存9257項と全係数一致することを検査した。実際のq次数は**G6が2、K_iが1**である。generic betaの特殊化による外挿はない。

全F_iとG6をq次数2へhomogenizeすると、q³の恒等消去から

    w*G6^h = sum_i K_i(q)*F_i^h

が全61変数+wのQQ多項式identityとなる。両辺を別に構築して全係数を比較した。両辺のq,w次数は3である。

従ってI134=(F_i^h)とすれば

    I134 subset I135=I134+(G6^h) subset I134:w,
    I134:w^infinity = I135:w^infinity.

これはradicalやsupportだけでなく実idealの包含・等号である。D(w)上では両idealが同じであり、未知点とnilpotentsを保持する。

## 2. 孤立axis因子とのactual composition

先のe_axisは元134のchart algebra内のidempotentで、w*e_axis=0を満たす。任意のquotientでもidempotentはdescendするので、I135へG6を加えた後もe_axisの像を使い

    J=I135+(e_axis)

とできる。D(w)ではe_axis=0なので、再び元の全アフィンschemeと一致し、w-saturated closureも変えない。

保存されたgeneric minimal-axis sourceは全135 leading equationsを満たすので、証明済みのaxis componentではG6^h=0である。producerは旧generic135証明書をpinし、実QQ sampleでもG6最高式を確認する。もしDelta saturationが他の孤立因子も選び、その一部でG6が非零なら、それらはI135により消えるだけである。残る因子のidempotentの像を使うcompositionはなお正しい。選ばれた全finite etale algebraでG6が零だと無根拠に仮定しない。

## 3. global quinticへの有限拡張

同じ凍結T01 betaでG6を元generic9257項から正確に評価・中心化し、保存790項と比較する。新global quinticのformal45係数環へ代入し、monic quinticで割る全QQ identityを再検査した。追加remainderはh次数≤4である。

5×670表示へこのremainderのh^j倍(j=0..4)五列を加えれば、同じ実Aを与える5×675表示となる。元134 affine idealの帰結なので新Poisson条件を追加していない。rank探索や高次数Macaulay再試行は行わない。

## 範囲と保存

`code/s61_affine_g6_partial_colon_20261103_v63.py` と `reports/S61_AFFINE_G6_PARTIAL_COLON_20261103_v63/`。全generic cofactor identity、実次数、quinticの追加remainder、axisとのsheaf compositionを保存する。これらは具体的なw-torsion処理の最初の二段階であり、残りraw無限遠の全除去、direction-chart B_iのQQ nilpotence、一般五点尽くしは未完了である。

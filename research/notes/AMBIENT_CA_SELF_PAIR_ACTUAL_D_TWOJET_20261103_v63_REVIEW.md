# 実self-pair上のactual D：全QQ(a)二次差とCB包含の独立監査

2026-11-03。対象は非空parameter openとそのG-saturationである。a=1だけの標本からの外挿は行わない。

著者のraw Lie/Hessian/係数変換ルーチンをimportせず、raw CA20を物理的にcanonical化し、P(a)、Q(a)と両全flag9を再構成した。
著者の有理kernel方向は共通非零因子587068342272 a¹⁸/78125で割り、小さいQQ[a]方向へ正規化した。
両GL生成子Aᵢ³=0を確認し、exp(tAᵢ)=I+tAᵢ+t²Aᵢ²/2と実逆行列による**全有限5×5 pushforward**を直接展開した。
全150の共通一次係数vと真二次差Lが、著者のQQ(a) certificateと正規化因子込みで一致した。
二次Lie公式を再利用しただけのチェックではない。

CBの全sourceはrawCB(f,mu*f,H)の実ratio列g=fを保持して再構成した。
固定19列、ratio1列、全9flag列を実CB座標から同一P(a)へ有限pushbackした。
QQ(a)上の独立rankは

    CA1+CB 30、CA2+CB 38、CA1+CA2 37、+L 38、さらに+CB 38

である。実GL×a familyはrank20、vを加えると21、CBにvを加えると30。
これは各包含と非包含の有理関数体上の恒等式であり、適切な非空algebraic openで成立する。

個別CAのactual ideal元Fを両実curveへ代入すると、二次差は同じHessian項を打ち消す。
従ってdF(L)=0であり、whole CA tangentはJ1,J2だけでなくLも含む。
上の厳密rankにより全CB tangentはwhole CA tangentへ入る。
二branch tangentの和37をwhole tangentと同一視してはいない。

Nはここでsmooth28であり、第一CA branchとCBの両方へ入る。
接交叉次元28とNの包含により、第一branchとCBの全local intersectionはregularなNである。
同様に、第二branchとCBは実smooth20 paired-image Sを含み、接交叉次元20なので、全local intersectionはSである。

reduced全projective交叉=N27∪R19、R∩N=B18、および**全B18⊂Sing(N)**はそれぞれ既存の独立監査済み証明を明示的に使用する。
本族上のNはsmoothなので点はB18から外れ、従ってRからも外れる。
よってactual Dのreduced supportはここでsmooth N28だけである。
actual D tangentはCBの全29次元接空間なので、dim28/edim29であり非被約。
射影ではdim27/edim28である。

このpacketでは全D ring、nil module、associated prime、a=0の型は結論しない。
これらは二CA branchの全81kernelとCBへの実制限による別のideal equalityを要する。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CA_self_pair_actual_D_review_20261103_v63.py --output reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW.json --replay <new-replay-directory>

fresh著者source方向・全150係数とrankを照合し、別の有限物理計算とQQ(a)行列rankで確認する。GBなし、120秒・2.2GB。

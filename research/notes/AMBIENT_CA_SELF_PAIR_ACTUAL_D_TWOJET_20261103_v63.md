# CA自己旗族でのactual CA/CB交叉：真二次差による非被約性

2026-11-03。actual多項式族P_aを対象とする。以下はQQ(a)での計算、従って非零aのある非空開とそのG飽和上の結論である。全特殊点やB18原点のactual Dをここから分類しない。

## 二つのCA branchとCBを同じphysical点へ揃える

第一CA source、第二CA sourceをg^-1で戻したもの、明示的CB sourceをその座標変換の逆で戻したものについて、全150行のsource微分J1,J2,JBを再構成する。全sourceのrankは29。QQ(a)上で

    rank(J1,J2)=37, rank(J1,JB)=30, rank(J2,JB)=38。

Nは一般点でsmooth28で、第一CA branchとCBの両者に入る。従って第一CA branchとCBのactual完成交叉はNそのものである。
GLとaによる20次元の滑らかな局所像Sは第二CA branchとCBの両者に入る。その交叉の接次元が20なので、actual交叉はsmooth20のSとなる。
これらは三つの既約成分の接空間から全Poisson schemeを推論する議論ではない。

## 共有接方向の真二次差

GL+aの20方向を越える共有CA接方向vを取り、J1 u1=J2 u2=vを全150行で解く。
各sourceで、実曲線 exp(t A_i)_*(P_i+t delta_i)を使う。
その二次係数は Lie(A_i)delta_i+Lie(A_i)^2(P_i)/2である。第二sourceを同じphysical点へ戻し、その二次係数の差をLとする。
計算では全150係数を保ち、QQ(a)上で

    rank(J1,J2,L)=38, rank(J1,J2,L,JB)=38。

a=1では、別の有限5×5座標代入による独立計算も全150係数で同一のLを与える。
疎な線形函数 ell=(2/5)e60+(1/10)e63+e68 は二接空間を殺し、ell(L)=−5/12である（この数値は固定点の対応する正規化に限る）。

## actual Dの接空間と非被約性

CAのactual ideal元Fを両実source曲線へ代入する。
その一次項ell_FはJ1,J2を殺し、同じ一次接方向vに沿う二つの二次式を引くとell_F(L)=0となる。
従って**whole CAの接空間はspan(J1,J2,L)を含む**。これは二branch接空間の和37をwhole tangentと同一視することとは異なる。
上のQQ(a)恒等式によりJBはその空間に入り、actual D=CA∩CBの接空間はCBの29次元接空間全体となる。

一方、全体の既証reduced交叉はN27∪R19で、R∩N=B18である。
この族のNは正則で、全B18はNの特異点だから、対象点はRから外れる。
従ってD_redはここでsmooth28のNである。
dim D=28、edim D=29なので、**actual Dは非被約**である。
射影版ではdim27、edim28となる。

## 次の全環判定（本packetでは未主張）

二CA branchのmeetingがdouble smooth Sであることと、その全fiber-product presentationを別に受理すれば、Dの全環には候補

    C[[s20,u8,f]]/(f²,f*u1,...,f*u8)

がある。これを証明するには、二branchの全kernelとCBへの制限を使った両ideal包含が必要である。
本packetのreportは、この全環、nil moduleの正確な長さ、全associated primesを未証明として保存する。

再現コード：`code/ambient_CA_self_pair_actual_D_twojet_20261103_v63.py`。元のQQ(a)全行列と二次差、正確なrank記録を保存する。

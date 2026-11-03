# CA実self-pair族の実N・CB source：旧定理との整合

2026-11-03。CAのordinary self-pairを、次元の類似からR19へ同定しない。
同じ実150係数でNとCBのsourceを明示すると、この族は全aで**実N**に含まれる。
従ってCAがR\Nでsmoothという既知定理と衝突しない。

## Nの実19係数source

CA critical packetとN critical packetは同じh,xを交換した物理frameを使う。
二つの保存済み行列を直接照合すると

    N150×19 = C150×20[:, (0,1,3,4,...,19)]

である。Nのfixed-fiberはC[2]を欠く。
先行実多項式族の第一CA source

    P(a)=P0−a*(3C[1]+C[9])

はC[2]係数が恒等的に零なので、保存した19個の多項式係数から全150行で実N sourceとして再構成される。
これは接空間や閉包の次元からの推論ではない。
N proper source自体の既知同定を用いる、具体的なsource incidenceである。

一方、第二CA source Q=g_*PのC[2]係数は−a²である。
従ってa≠0では**第二CA旗を同じN fixed-fiberの旗として扱えない**。
CAの二つの旗からNの二つの旗を推論するのは誤りである。
この区別により、B18でのN subintegral normalizationと、N self-flag jetの先行不整合はそのまま保持される。

## 実CB sourceの全150係数恒等式

N⊂CBという既知全交叉定理だけに依存せず、a≠0で実CB sourceも与える。
物理座標(h,x,y,z,w)に対し、変換Gを

    h_new=−3a*h−x, x_new=x, y_new=y, z_new=z, w_new=w

とする。det G=−3aで、QQ[a,1/a]上の実逆行列を持つ。
CBの標準raw式

    rawCB(f,0,H) = f∂x∧∂y
       +∂y∧[(-f_x*z/2−H_w)∂z+(-f_x*w/2+H_z)∂w]

へ

    f=3a*x²,
    H=−a*x*w*z+x*w*(h+x)/(3a)+(h+x)²*z/(9a²)+z³/3

を代入し、canonical化すると

    G_*P(a) = canonical(rawCB(f,0,H))

が全150係数でQQ(a)上の恒等式になる。
実5×5行列とその逆座標による有限代入で確認し、25 divergence・350 cubic Jacobiも零である。
これはPaの一つのQQ標本や有限jetからの外挿ではない。
a=0は既知の共通B18点であり、proper閉像の閉性とも整合する。

## 幾何学上の範囲

この族は実Nの中にあるため、R19へのmembershipを示す必要はなく、ここでは示していない。
CAの二つの旗は実際に異なるが、そのうち一方だけが上のN sourceとして表示される。
個別CAの二枝と、Nの一枝・subintegral modelは異なるsource情報である。
R19のacross-CA/CB pairを個別CA conductorと混同する旧推論は復活させない。

このpacketはN/CBの実包含だけを追加し、CA原点の全odd ideal、全conductor、全枝数やactualDの型を分類しない。
Paの族の19次元PGL飽和と、既知R∩Nの18次元という比較による一般的な非R判定は別の幾何学的帰結であり、
このsource identityの実行結果には含めない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CA_origin_self_pair_N_CB_incidence_20261103_v63.py --output <new-output-directory>

codeと全150係数:
`code/ambient_CA_origin_self_pair_N_CB_incidence_20261103_v63.py`、
`reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63/`。
実行約1秒、制限60秒・2.2GB、GBなし。

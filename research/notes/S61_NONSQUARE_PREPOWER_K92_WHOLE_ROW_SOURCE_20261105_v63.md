# K92：係数層を増やさず全Schur行の恒等式を調べるsource

保存済みの実Tbar550次prefixはTrow2=59375 t^68 Trow1を満たす。元の重みの差alpha0−alpha2=68を戻すと、重みを付ける前のS=D−CA^{-1}BでS2=59375 S0が候補となる。raw C,D行自体は比例しておらず、A cleanupを省略しない。

Aは全多項式としてblock lower2+2である。P,Qを対角二角、A21を左下二角とし、deltaP,Qを行列式とする。dc=C2−59375 C0、dd=D2−59375 D0に対し、v=dc_right adjQ、u=deltaQ dc_left−v A21、w=u adjPとすると、全七成分の分子はdeltaP deltaQ dd−w Btop−deltaP v Bbottomである。全sourceをLで整数化すると共通分母はL^5であり、未知の多項式逆を計算しない。

SOURCE段階では指数bit集合と係数絶対値上界だけを伝播する。後続jobは七成分を一つずつ計算し、最初の非零で止める。全七成分零なら、この固定11minorについて全多項式の零が証明される。非零成分一つは候補行恒等式を棄却するだけで、minor非零や全域57を意味しない。

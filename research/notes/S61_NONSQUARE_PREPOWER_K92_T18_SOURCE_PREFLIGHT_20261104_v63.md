# 次の固定G92/target18：source recipeと必要な比例gate

2026-11-04。新しい係数計算を行わない。残余二行のcurrent15=9（G92 lineage）と
original32 target18を一つだけ固定する。既存source見積もりでは、未検査h9二列のうち
18は17より小さい。今回17の新しいsource値は計算していない。

raw pair比例は今までの七列[22,11,13,25,15,19,14]だけで証明済みであり、18へ延長しない。
固定列18についてκ一式と三つのraw交叉式の支持・高さ・lifetimeを別々に見積もる。
未補正literal-pair式を使うには三raw交叉零という強いgateが必要で、κ零だけでは足りない。

G92のmiddle-heavy三列を消す実sourceは

    K92=Delta407 G92/c8−lambda505 G85/c6
        −a lambda419 G129/c12−lambda428 G99/c9.

G129のrhoだけaを含まないため、その係数にaが必要である。全ciは非零QQ定数。
もし18のraw交叉三式が零なら、pure八列上でpair四vectorが比例し、K115、K120、K92の
残るpair-heavy消去の寄与をliteral G20方向へ吸収できる。

さらにA=G115/c10、B=(G120−G99)/c11、D=G92/c8、C=G85/c6、
E=a G129/c12、F=G99/c9と置く。source十一行[P4,G20,A,B,D,C,E,F]と
pure八列＋middle-heavy三列[28,30,29]の行列式det11を定義し、最後三列だけ既知単元L=aFC²で割る。
sourceを最後三行で消すと、K115/g、K120/g、K92/Deltaとなり、最後のH3行列式はDeltaである。
従ってpure八minor D8=g二乗 det11というpolynomial恒等式が得られる。
ここでgやDeltaを数値的に割って係数を作る必要はない。

条件付きで次の残余成分はc8 g det11/(Delta a二乗 N7)となる。
det11の係数・非零・rankは未計算で、whole-worker資源gateもまだ通していない。
この小さいraw18 gateが先に必要というsource frontierを保存する。
全tail2効果・七remaining tail・旧chartの全補集合を保持する。

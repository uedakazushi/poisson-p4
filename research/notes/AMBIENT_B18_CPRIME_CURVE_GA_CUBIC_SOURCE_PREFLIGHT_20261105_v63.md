# 実 conductor 曲線列による残余三次一汎関数

既知 cubic closure は 27 行中 26 行であり、一つの残余汎関数がある。本 packet は固定 minor からその汎関数を取り出すが、以下の実曲線列との pairing はまだ評価しない。

honest C0 の既証の座標 telescoping 恒等式は f63=gA+Σz_i e_i であり、gA は七 zero 座標を零とした pointed 57-unit inverse 上の関数である。gA は r,s だけに依存し、次数 3 未満の係数は零、三次係数は (3/8,2,13/6,2/3) と既存 cache で確定している。

比較行列 V の最後の八列は −U_i である。元の C0 frame が曲線 gA を使えば −gA tauA が文字通り最後の列となる。raw f63 を使う同値 frame でも、−f63 tauA から Σe_i(−z_i tauA) を引く unipotent 列変換で同じ曲線列を得る。従ってこれは実 Cprime の関係であり、raw f63 の非 plane 二次項を無視する仮定ではない。

s=r+ell/5 により gA3=125 r³/24+5 r²ell/3+r ell²/6+2 ell³/375 である。γ=(tauA+11 tauB/12,tauA tauB,tauB) に対する係数は (−gA3,0,11gA3/12) となる。この列の次数 3 未満は全九変数で零なので、その三次式自体が真の initial relation であり、追加の V1 cancellation や source3・gauge3 は不要である。

固定残余汎関数とこの保存済み八項を pairing するだけで、非零なら既知 cubic rank は 27 になり、Cprime/ell³ 長さの上界 9 と既存下界 9 が一致して e=2 を証明できる。零ならこの実列では残余を消せないだけで、e>2 や真の cubic quotient の非零は推論しない。値計算は別の bounded GO と独立監査の対象とする。

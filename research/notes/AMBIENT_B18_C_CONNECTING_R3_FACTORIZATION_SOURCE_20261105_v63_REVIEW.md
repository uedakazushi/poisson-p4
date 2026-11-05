# 元の Told 上の r³ と接続写像：独立 SOURCE 検証

独立検証は PASS した。報告 SHA256 は
`b8e4dc64f6133a0b98b1c61175eb55d840686543b46d267eac6604c1a9ebae00`。
実 Told の HF `(5,3,2,1)`、長さ11、m⁴Told=0、実七次元商 M7 と
接続写像の受理済み source を固定し、新しい head の値は計算していない。

I Told の実四基底は `z5*h4,z5*r*h4,r²*h4,r³*h4` である。
実三次 initial source `qC=z5*r²*h4` の持上げは qC と m⁴ の項の和であり、
m⁴Told=0 によって qC は実加群内で零である。他の三基底の r² 倍は次数4以上なので、
r²(I Told)=0 が実際に成立する。単に initial cover 内での零ではない。

M7 の実関係から `rh0,rh1,rh2−rh3/4` は I Told に属する。
したがって r³h0=r³h1=0、r³h2=r³h3/4 となる。
G38 の実持上げ `r²h3+(11767/3652)z5rh4+O(m³)=0` を r 倍すると、
残りは m⁴=0、もう一項は qC=0 なので r³h3=0。
残る r³h4 は実 degree3 基底 gdeep で非零、かつ m がこれを消す。
よって全五 head の r³ 像は gdeep 座標で `(0,0,0,0,1)` である。

genuine G=C[ell] の cyclic top functional は
`chi=(0,0,0,0,d)`, `d=−36/1681` である。
全五 head が同じ Rbar 上で Told を生成するため、実 Rbar 線形写像として
\[
r^3:Told\longrightarrow k\,gdeep
= (gdeep/d)\,\chi\,\delta
\]
が全加群上で成立する。ここで表示の右辺は写像の合成を表し、
未知の高次 head 係数を零に置いていない。

元の正確列の kernel は正直に Kbar=K/ell K なので、r³Kbar=0 を得る。
元の K へ戻した結論は **r³K⊂ell K** である。ell を約して r³K=0 としてはいけない。
Kbar の全 extension、元 H0 の DVR 指数、nil 三乗や正確な nil 指数は決めていない。

元 primary の code・note・全入力・受理済み peer を hash 固定し、
独自の source reader と Fraction で五 head の係数を照合した。
所有数上界13,703、実行0.0062秒、2 GB・55秒の guard 内で終了した。
再現コードは `code/ambient_B18_C_connecting_r3_factorization_source_review_20261105_v63.py`。

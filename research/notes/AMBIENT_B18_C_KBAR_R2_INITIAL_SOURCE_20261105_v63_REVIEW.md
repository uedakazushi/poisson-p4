# Kbar の r² の最低次射影：独立 SOURCE 検証

独立検証は PASS した。報告 SHA256 は
`662785338116da297cfe47ea8b8e117fee0b4ee2834d808458600f59f95384c7`。
新しいsource jetや完全作用行列は求めていない。

元の37個のlinear initial generatorが、既存linear artifactと文字通り一致することを確認した。
producerのGroebner divisionとは異なり、三つのstandard slotを除いた
固定37×37の転置係数系を、下側row優先の有理数消去で解いた。
その結果、gr¹で
`rh1=alpha1*z5h4`, `rh2=(rh3)/4+alpha2*z5h4`、
`alpha1=−1681/1992`, `alpha2=75645/14608` を得た。
これはinitial/associated-gradedの式であり、filtered Toldの高次補正を捨てた実恒等式ではない。

独立なsource cofactorを作って全40成分へmultiplybackした。
さらにproducerの二つのlinear certificateと七つのdegree2 certificateも、
元40 generatorとそのmonomial cofactorだけから逆順でmultiplybackした。
全て入力に一致した。

G38の実initial関係 `r²h3+(11767/3652)q1=0` を使い、
`q1=z5rh4`, `q2=r²h4` 座標で七基底のr²像をmod old m³へ射影すると、
\[
\begin{pmatrix}
0&0&0&0&0&273412940423/2952943290&-6957653957/295294329\\
0&0&0&0&0&0&0
\end{pmatrix}
\]
となる。受理済みexact Told HF²=2により、q1/q2は単なる既知coverのstandardではなく、
真のgr²基底である。よってこの射影のrankは正確に1である。

実際のr²像には、old m³の非零基底gdeep方向の補正があり得る。
それらは全て未計算のまま保持し、上の行列を完全なr²作用行列とは呼ばない。
一方、非零のgr²射影からr²Kbar≠0が従う。
既に独立証明したr³Kbar=0と合わせると、rのnilpotent blockに少なくとも一つ長さ3のものがある。
完全なr-Jordan分割や、異なるDVR変数ellの指数・元H0・nil三乗は決めない。

独立SOURCE所有上界20,092、before-operation最大115 bit、実行0.0363秒、
55秒・2 GB・RAW32 MiB内で終了した。
再現コードは `code/ambient_B18_C_Kbar_r2_initial_source_review_20261105_v63.py`。

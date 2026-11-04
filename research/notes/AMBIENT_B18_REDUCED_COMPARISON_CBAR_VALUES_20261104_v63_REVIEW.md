# actual reduced Cbar 長さ4：独立監査

actual Cbarは既受理の3 sector `k[r]/r³ ⊕ k[r]/r³ ⊕ k[r]/r²`。
全z作用は第1・第3で0、第2で `zi=-ci r²`、非零ciは0,4,5のみ。
B11の全定義関係と全363 basis積が各sectorで消え／積互換となることを独立に確認した。
従ってB11で計算した係数関数をsectorへ送ることはexactであり、欠けている高次係数を
任意に0と置く操作ではない。

自分の先行独立B11/C44 ring実装と保存されたactualκ9だけを再利用し、
`n1=d(a-d+x-2t)` と `n2=d(d+x)(a+d-x)/6` をC44で計算する。
全2×4 basis倍の8像を、その後にR0-linearにCbarへ送る。
Cbar自体を環やHensel計算の底とはしない。

8像の全一次・二次R0倍360列を独立に再構成する。元Cbarでm³=0を作用表から確認し、
固定4列／4pivotによる像の正確なspanと全8作用に対する不変性を検査した。
新しいpivot探索は行わない。正規形は全作用と整合し、商は長さ4・HF(3,1)。

商基底は `(τA,τB,rτB,τAτB)`、全zi作用0、`rτA=−11/12 rτB`。
基底を変えると正確に `k⊕k⊕k[[r]]/(r²)` となる。
ここでrはellで切った後のCA-F9座標であり、DVR parameter ellとは異なる。

whole E4のnilidealは既受理の二生成n1,n2であり、free4係数環から8個のR0-module像を
使えば全像になる。`Cprime=C/image Nil(E4)` のell-basechangeはright-exactなので、
この像の余核が正確なCprime/ellである。kernel/Torを落としていない。
generic DVR rank、torsion指数、元nilideal生成元やcubeはこの有限商から推論しない。

再現：`.venv/bin/python code/ambient_B18_reduced_comparison_Cbar_review_20261104_v63.py`

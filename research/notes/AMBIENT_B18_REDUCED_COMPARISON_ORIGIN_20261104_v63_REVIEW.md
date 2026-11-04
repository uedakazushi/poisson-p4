# reduced比較の原点rank0：独立監査

actual proper/physical Jacobianのtrace二列を読み、独立検証済みB11/C44保存κ9の
同じ一次係数とも照合した。`a=τA−τB,d=−2τA−τB/2,x=t=τA`。
五pivotはtrace一次部0で、保存された五rational core rootの分子も全て次数2以上。
従ってproper u5は次数2以上である。

原点係数環は `k[τA,τB]/(τA²,τB²)` で最大ideal三乗0。
h一次部は `2τA−τB/2` なのでdhの係数は `(-4,0,1/4)`、従ってdh=0。
d*u5とdcr,dcl,drlは三次以上で0となる。全9 native関係と全36 basis倍が0。
実C0のβ8も原点0なのでVrad原点は全零、従ってcokernelのμ_R0は3である。

これはright-exactなcokernel fiberの計算であり、kernelやH1のbasechange、
元nilidealの生成数・cubeを推測していない。

再現：`.venv/bin/python code/ambient_B18_reduced_comparison_Cbar_review_20261104_v63.py`

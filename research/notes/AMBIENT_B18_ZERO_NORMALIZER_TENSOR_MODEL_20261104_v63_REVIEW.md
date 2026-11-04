# Actual normalizer tensor模型と導体積のTor：独立監査

対象は `ZERO_NORMALIZER_TENSOR_MODEL` と `CONDUCTOR_PRODUCT_CURVE_BASE` の凍結packetである。著者builderを実行せず、新しいsource solve・jet・GB・旧r/g cofactorのpullbackを行わない。

## 同じP0と、全free4/48入力

物理的150係数のzero-weight40方向について、25本すべてのdivergence線形式を独自に作り直す。そのrankは4で、表示された四式と同じrow spaceを持つ。残り四座標を復元するminorは2であり、36座標は自由である。common cut15を除いたP0はregular21変数、同じphysical R0九座標と、zero odd八座標・zero graph四座標で与えられる。後四座標を落とさない。

pruned DAGの元nodeへの対応を、定数・符号・和・積・単位逆元の局所的な普遍環恒等式で全件照合する。非零分母は原点値でも確認する。これは巨大式のhash評価ではなく、各nodeの意味を帰納的に確認する手順である。依存閉包にはc/v/a/m/n/i以外が無く、旧r/gは一つも入らない。

旧expanded HONEST/C0 archiveには、前段レビューが確認済みの未使用duplicate variable aliasがある。physical xの名前を最後の同名aliasへ結ぶことはせず、元のpure sourceのnode IDへ固定する。

pure CA48、CB48、CA tail-prefix57の計153単位式と変数順を元sourceに固定する。三つの既受理のunit blockが三角形に並ぶため、同じphysical R0上の係数函数を全次数で指定する。free4のbasisは `1,tauA,tauB,tauA*tauB`。両tau乗法行列の可換性と各平方νを照合し、全12 target差と全四basis倍、計48列を独立に確認する。basisがfree4を張るので、これら四倍は理想全体のR0-module表示である。C0 free9の自然写像も同じ係数根で一致する。

## Tensor比較における正しい零化ideal

各 `A_i=R0+I_i*taui⊂B_i=R0+R0*taui` で、f∈I_iによる乗法はP0-linear map `B_i→A_i` を与え、inclusionとの両向きの合成がf倍になる。任意のcompleted tensor base change後にも、この二つのmap恒等式は保持される。従って第一のnormalization stepのkernel/cokernelはI_Aに、第二はI_Bに零化される。二段の合成については、kernelとcokernelの双方が **I_A I_B** に零化される。flatnessやTor消滅を仮定しない。

非負gradingのdegree-zero商は `D0=A_A completed_tensor_P0 A_B`、全normalizer比較対象は `E4=B_A completed_tensor_P0 B_B` である。finite integral normalizationのSpec射は全射で、この性質はbase change後も保たれる。したがって両者のR0-supportは同じだが、kernelが零であることや深さは従わない。

**この積はT0全体のannihilatorではない。** また、各conductor plane上で保存されたN56一次補正を独立に代入し、全56式を確認すると、v_Nの一次式は両側とも零になる。提案されていた一次式とg3のcoprimalityによる有限性の十分条件は通らず、比較kernel/cokernelの有限性はこのgateからは分からない。

## CI8二曲線のKoszul計算

同じR0の二つのactual zero odd idealは各々CI8であり、和は実際に

```
J=(z1,...,z7,r²,rs,s²), R0/J=A3.
```

Jのminimal generator数は10である。16本のcombined生成元からJ/mJのbasisを選び、残り六本をそのexactな線形結合として引けば、local GL16変換で六つの零を得る。さらに任意のminimal10組の間の変換はresidue fieldで可逆なので、Koszul16は `(z7,r²,rs,s²,0^6)` のKoszulと同型である。これは正確なideal等号による存在証明であり、巨大な16×16変換行列を計算したとの主張ではない。

regularなz7を除くとS=QQ[[r,s]]上のK(r²,rs,s²)となる。一次cycleは自由basis

```
u=(-s,r,0), v=(0,-s,r)
```

を持つ。実際、Ar²+Brs+Cs²=0なら、rとsの互いに素性からC=rD、A=-sE、B=rE-sDとなる。二次boundaryは `r*u, s*u+r*v, s*v` だからH1はomega_A3に一致する。uをdual basis eps_s、vを-eps_rへ送ると全三関係が一致し、長さ3・socle1となる。

H2も直接確認できる。`rA+sB=0, rB+sC=0` の解は `(A,B,C)=(-s²,rs,-r²)E` で、K3の三次boundaryそのものである。H3はdomain上の非零列のkernelなので零である。従ってH0=A3、H1=omega_A3のみが残る。

六つの零生成元のexterior factorを戻すと、非canonicalなA3-module同型

```
Tor_i(R0/I_A,R0/I_B)
 = A3^binom(6,i) direct_sum omega_A3^binom(6,i-1)
```

を得る。長さはi=0,...,8について `3,21,63,105,105,63,21,3,0`。
特に `(I_A∩I_B)/(I_A I_B)=A3^6⊕omega_A3` は長さ21、socle dimension13、annihilatorは正確にJである。

## Product ring自身と、任意の比較moduleを区別する

二つのCI curveはCM1であり、`R0/(I_A∩I_B)` はその直和へ単射に入るから有限台torsionを持たない。従ってscheme-theoretic unionもCM1である。各idealがradicalとは仮定しない。

`0→(I_A∩I_B)/(I_A I_B)→R0/(I_A I_B)→R0/(I_A∩I_B)→0` から、product ring自身のH0は上の長さ21moduleに正確に等しい。しかし、このring上の任意のmoduleのH0や長さに21という上限は与えない。parameterの任意高冪で割れば長さは増やせる。

比較kernelの有限性・非零性・socle、E4やD0の深さ、T0の全supportは未決定のまま保持する。専用二review reportに、全依存hashと実測時間・メモリを記録する。

# Actual top12 graph circuits：正則telescopingの独立監査

2026-11-04。**六membership lift、同じCA trace基底の六relative graph pair、
十二top pairの有限代数回路を受理する。** T3所属・socle・depth・nilradicalは未判定。
先のblocked preflightのbytesは変更しない。後続のfull186-generator packetも今回の対象外。

## Actual座標gate

前stageの独立source/DAG/Hensel48監査を再利用し、そのhashを固定する。
新checkerは元41,315 nodeを独立Fraction ADで微分し、実CBの48式について
`z_X=-H_z^-1 H_X` を別のexact Gaussian実装で求める。
保存された全48×48 Jacobian、48×9微分、実10×9 odd Jacobianを全係数一致させる。
source変数をphysical Xと誤認していない。

baseは実physical rows `27,108,109,112,113,117,118,119,61`。
odd rowsは `18,19,28,29,60,67,72,63,46,95`、その全origin値は零。
保存された選択index `0,1,2,3,4,5,6,8,9` についてdet=1/96を独立に確認する。
row63の一次部分は零である。
従って選択九本F_jはR0の正則座標であり、十本のidealは最大ideal全体。
単に九変数上に十式あることからrankを仮定していない。

## Pointed inverseと正則な除算

未知順(z48,X9)、式順(H48,F9-Y)のJacobianはSchur formulaにより
`detH_z*det(dF/dX)=(-331776/625)*(1/96)=-3456/625`。
よって一意なpointed algebraic-power-series inverseがある。

全176,910 nodeをread-onlyで検査し、元DAG prefixが完全一致することを確認する。
必要ancestor部分を、各levelの57変数bindingへ**独立に構造置換**することで、
九instanceの全513式とcubic出力を保存nodeへbindする。author/sourceの再実行は無い。

level iのtargetは(Y1,...,Yi,0,...,0)。Yi=0に制限すると、level iとi-1のtargetは一致する。
同じ原点を通るunit57系の一意性により、二つのcubic出力も一致する。
level0ではactual cubic b(0)=0。
従って各差は正則座標Yiで割り切れる。これはYiが単元という主張ではない。
`R0/(Yi)` への制限kernelが(Yi)であり、R0はdomainなので商は一意な正則power series。
分子とYiが代数的なので、その正則商もrational function field上代数的である。

保存された全54個の`r` nodeについて、分子、座標、level/alpha、原点値を照合する。
原点値も別に計算したdb/dYから得る。
通常のinverse nodeは全て非零origin-unitのみを反転し、`r`を無条件の非単元inverseとして扱わない。
方程式Yi*h=numeratorだけに任意のbranchを許す読み方もしない。

各Yi*hは二つの連続出力の差なのでexactにtelescopingする。
full level9へY=F(X)を代入すると、元X,z自身が同じunit57系を満たす。
一意性から最後の出力は元actual b(X)。従って六membership全てが全次数で成立する。
除いたindex7の係数は零にできる。finite jetをzeroにする推論ではない。

## Physical graphとtrace alignment

共通minimal targetはbase10とodd10の同じphysical座標である。
zero oddのgraph liftにはw³、positive oddのgraph liftにはw²を掛けるので、
どちらも正しくweight3。hは同じphysical零base R0上のscalar函数である。
CB自身へ代入すると、その六raw cubic outputを回収する。

CAへ代入する際は、CA自身の `tau_A²=nu_A` を用いる。
CBのtau_Bをtau_Aと同一視せず、同じraw physical qのeven shiftの差を全て保つ。
checkerは六つの(a,b)について全even/odd式をread-only DAGで再構築し、全output IDを一致させる。
さらに `j_A*tau_A*(a+b*tau_A)=(nu_A*j_A*b,j_A*a)` から十二pairを全て照合する。

graph liftの選択差はCB kernelへ入るため、full CB kernelとlower graph relationsを残す限り
交叉idealは変わらない。それらを捨て、十二pairだけからT3全関係や所属を結論してはいけない。

今回の受理はactual係数回路の構成までであり、四加群socle判定、depth19/20、全nilradicalを
決定していない。新しいsource、点、GB、full186生成元jobは実行していない。

再現:

```
python code/ambient_B18_top_graph_telescoping_review_20261104_v63.py \
  --output /tmp/B18-top12-independent-review
```

# P*のCPLUS親排除：whole generic rank上界の正確な旧証明入力

2026-11-08。純入力監査のみ、旧算術を今回再実行したとは呼ばない。
`UNIMODULAR_POLE_GLOBAL_MARKED_JET_AND_T36_T38_EXCLUSION.md` §2のCPLUS排除を、
一つの固定点のD100からgeneric上界を外挿する主張と混同しないための補足。

旧 `legacy/v08/code/cplus/certify_cplus_component.py` は、全(c,beta)族Fの全350×150微分Dと
全group/scale＋二source微分の行列Aについて、DA=0をpolynomial恒等式として確認する。
`universal_normal.json` はbidegree(2,1)の全150係数Ruを保持し、同codeは

```
D(c,beta) Ru(c,beta)=0    （全350係数恒等式）,
rank[A(-1,1) | Ru(-1,1)]=50
```

を確認する（同code44–51行）。後者は50rankの非空openを与える。
従ってsourceのdense openでkernelDは少なくとも50、rankD≤100。
この弱い上界は全101角minor消滅というclosed条件なので、
同fullsourceの全PGL飽和閉包CPLUS_EXTにも延長する。
固定点D100は一般階数の下界に使う旧入力であり、本排除の上界には不要。

旧certificateは`field=QQ(c,beta),generic_Jacobi_rank=100`、universal tangent identity、
sample_actual_plus_normal_rank50を明示する。
これらのファイルは今回改変も再実行もしない。旧全source証明の入力を明記している。
P*のfullD rank101は新同一点証明書であり、このclosed rank≤100と矛盾する。
従ってwhole CPLUS_EXTはP*の親ではない。
Ga9 cusp24の一つのD100など、別族へこの全source上界を流用しない。

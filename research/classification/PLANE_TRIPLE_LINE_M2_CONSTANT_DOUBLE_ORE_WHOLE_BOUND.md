# R=3L・m二constant double残二源：全P2九を保つOreの短い全源証明

2026-11-08、root純手証明、独立確認前。
八fixed P0全幾何還元は別作者受理済み。六sem源は全bad u-row零で空。
本稿は残NIL0/NIL1のfullP1 kernelに**P2を全九係数自由に戻した**必要supersourceを用いる。
homogeneous K2だけに限定した旧NIL0 a-factor推論は不採用で、別erratumに明記する。
新CAS、operator、minor、quadratic jet又はGBは実行しない。

## 1. 全kernelからcを除くEuler gauge

座標は(a,b,c,u,v)。両nil P0は

```
{b,c}=3a², {c,u}=a*u, {b,v}=a*u,
{c,v}=epsilon*b*u+a*v, epsilon=0(NIL0),1(NIL1),
other brackets=0.
```

全K1は不変八source reportのcomplete kernel basisを用いる。
各basis rowはcodeのfull30 basis順であり、coefficientsをt0,t1,...とする。
全c依存はP_bcのc*u/c*vの二項だけである。
P_bcのc係数をB(u,v)とし、linear fieldT=B(u,v)*partial_bを取る。
P→P+E∧Tは同じprojective tensorで、P_bcの−cBが全c依存を消す。
他の追加成分はP_ab+=aB、P_bu-=uB、P_bv-=vBだけでcを生じない。
fullP2九は初めからc非依存なのでこの計算に影響しない。

従って全sourceは文字通り

`P'=Q(a,b,u,v)+partial_c wedge W(a,b,u,v)`

となる。全W係数と元sourceとの相関を保つ必要包含であり、Wを固定していない。
fullP2のab成分は任意F2(u,v)をQ_abへ加える。ac/bc成分は全てWへ入る。
特にP2のinhomogeneous particularを失っていない。

全projective Jacobiをc次数で比較すると

`[Q,Q]=0, E_y wedge[Q,W]=0, y=(a,b,u,v)`

が必要である。strict条件は**Qについて**だけの結論で、元PのstrictJacobi零を追加仮定しない。

## 2. constant factorの全点補題

Q=partial_b wedge V、Vがb非依存ならPfQ零で、actualFはc非依存、
全体signを除いてF=det(E_y,partial_b,V,W)である。
module式E_y wedge[Q,W]=0へpartial_bをwedgeすると

`det(E_y,partial_b,V,partial_b W)=0`。

E_yのb微分はpartial_b、Vのb微分零なので、これはpartial_b F零に等しい。
ゆえにFはb,c両方に依存せず、essential高々三。
非零既約homogeneous quinticがessential二以下になることは無いので、
対象は既全essential-three八23成分へ全pointwiseに被覆される。
Wがhomogeneous kernel、particular又は全十八parameterに属するかはこの証明に不要。

## 3. NIL0：全必要Qがconstant factor

reportの十二K1係数をt0,...,t11とする。
B=t2*u+t5*vである。全P2九を戻したEuler gauge後の全Qは

```
Q_ab=a*(L*u+M*v)+F2(u,v),
Q_bu=-t2*u²+(3*t7-t5)*u*v,
Q_bv=a*u+t8*u²+(t10-t2)*u*v-t5*v²,
Q_au=Q_av=Q_uv=0,
L=(t2+t4+3*t6+t10)/2, M=(t5-t7)/2.
```

全Qはpartial_b wedge Vで、全coefficientsはbに依存しない。
従って§2から、**NIL0の全fullsolutionはessential三の八23成分へ入る**。
NIL0全空とは主張しない。旧a-factorの誤推論を使わず、全P2九を保持した全点E認識である。

## 4. NIL1：一つのstrictQ二係数だけで全sourceを二枝へ

reportの十一K1係数をt0,...,t10とする。

```
A=t3/3+t4-2*t6-t9/3, T=t6, J=t9, R=t7, m=t5,
B=(t3/3+t4-2*t6-t9/3)*u+(23*m/3)*v,
L=2*A+3*T+J.
```

全P2九を自由に戻したQは

```
Q_ab=a*(L*u+(10*m/3)*v)+F2(u,v),
Q_av=T*u²,
Q_bu=-A*u²-(14*m/3)*u*v,
Q_bv=a*u+R*u²+(J-A)*u*v-(23*m/3)*v²,
Q_au=Q_uv=0.
```

F2は任意binary quadraticであり、仮想的に零としていない。
strict Jacobiのabu/abv係数を元定義から直接計算すると

```
J_abu=-(14/3)*T*m*u³,
J_abv=-T*(A+3*T)*u³-(28/3)*T*m*u²*v.
```

F2はu,vが互いにcommuteしa,uもcommuteするため、この二式に寄与しない。
従って全reduced complex solutionでT*m=0、T*(A+3*T)=0。
T零ではQ_avも零、全Qはpartial_b wedge Vでb非依存となり、§2から全essential三被覆。
T非零ではm零。元全P1の唯一bad u-row v²Cはm*v²*partial_cだから、
m零は元full u-rowがu倍であることを意味する。
Euler補正T=B*partial_bもu-rowにu因子を足すだけで、この条件を変えない。
全P2九にはnormal columnsが無いので、任意のaffine particularを含めactualFはuで割れ、既約性に反する。

これでNIL1の全sourceも、空又はessential三の八23成分へ全点被覆された。
T,m,Aの零境界及びF2の全係数を落としていない。

## 5. 結論と検証状態

六sem全空と、上の二nil全点被覆を合わせれば、**R=3L・m二constant-double全sourceは空又は既八23成分へ入る**。
これは単なるG次元上界より強い全点帰属であるが、独立確認を受けるまではproducer純証明の状態を保つ。
証明に使うK1全basisの完全性は別作者exact replayへ依存する。
八reportの誤metadataやNIL0 a-factor flagを使わず、Qへの全包含・全P2九・全projective方程式のc次数を使う。
他のR、SF00/01/02又はm三・全G/Eは本稿で閉じない。

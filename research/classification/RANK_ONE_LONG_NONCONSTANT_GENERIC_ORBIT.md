# rank一long・非定数g成分のgeneric軌道22と商三次元

2026-11-08。主成分証明及びそのimmutable入力を変更せず、全源のactual stabilizerを補足する。統合担当の全parameter手導出をここで係数から独立に確認した。新しいambient点、minor、profile又は計算jobは使用しない。

## 1. 全Sに二つの独立affine Poisson場がある

`RANK_ONE_LONG_NONCONSTANT_COMPONENT_ENTRANCE.md` の全sourceはX=∂x+x∂yをmodular場に持ち、[X,P]=0。もう一つ、g≠0上で閉一形式dlog gに対応するHamiltonian場をcolumn-sharp規約で

```
Yg=P#(dlog g)
  =d∂x+[d*x+K_w*g_z−K_z*g_w]∂y
      +g_w∂z−g_z∂w,     d=div U
```

と置く。全係数の導出はU(g)=d*gとH=gKから直ちに得られる。例えばy-componentは

```
[Pyz*g_z+Pyw*g_w]/g
=d*x+[(gK)_w*g_z−(gK)_z*g_w]/g
=d*x+K_w*g_z−K_z*g_w.
```

gはaffine、Kは二次以下なのでYgの係数はaffine一次以下。従ってg零集合にもregularに延長する。g非零上でHamiltonianでありL_Yg P=0なので、その多項式恒等式は全affine空間で成り立つ。

gの線形部分は非零なので(g_w,−g_z)≠0。一方Xのq-componentは零。従ってX,Ygは全Sで独立なaffine vector fieldsで、P4上にも独立なprojective linear fieldsとして延長する。h-rowが零の二場の線形結合がEuler fieldになることもない。

char零のalgebraic stabilizerはそのLie代数と同じdimensionを持つため、全Sでprojective stabilizerは少なくとも二次元。PGL5軌道の次元は全域で22以下である。これは一点の軌道rankから全域の上界を推測したものではない。

## 2. 保存済み一点の下界とgeneric値

既実行した唯一の点のfull Euler／GL／scalar行列Mはrank47、Euler24、scalar一なので、その実PGL5軌道は22次元。既nonzero47-minorによりrank≥47はsourceの非空開条件である。前節の全source上界と合わせ、25次元成分Zのgeneric PGL5軌道次元は正確に22。

従ってgeneric reduced Rosenlicht quotientのdimensionは25−22=3である。特別点でも全軌道が22だとする主張ではなく、quotientのcompactification又は全軌道分類を追加していない。

## 3. generic実因子の幾何

generic D=h²*L*QでQはsmooth rank五quadric、h=0は唯一double hyperplane。g=zの表示では

```
Q=h*y−x²/2+K2(h,z,w),
p=[h:x:y:z:w]=[0:0:1:0:0].
```

Q(p)=0でdQ(p)=dhなので、h=0はQのpにおける接hyperplane。L=g_hom=0もpを通る。gの線形部分が非零なのでLはhと異なり、generic LはQの非接hyperplaneである。例えばK=z²+w²,g=zではQ|L=hy−x²/2+w²がrank四となり、この開条件が空でない。

接点p、double-H、simple-L及びQの特殊衝突は同じ全source閉包へ保持される。このgeneric記述はそれらの境界に新たな非零条件を課さない。

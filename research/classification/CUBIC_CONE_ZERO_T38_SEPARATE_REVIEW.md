# Twisted cubic-cone rank零全源の T38 帰属：独立手証明確認

2026-11-08、PRIMARY `E_reduced_nonnormal`。
**判定 PASS：被約既約 actual quintic を持つ指定全 rank零 cubic-cone source は T38 に入る。**
全 kernel・全軌道の旧独立証明書を再利用し、新しい数値計算、rank、minor、M1、CASを使わない。
今回新しく確認したのは全六 vertex source の literal T38 同定、
全 nilpotent fiber への vector-bundle 閉性の橋、Q2 全 fiber の実一次因子である。

## 1. 全入口と残る Jacobi 条件

`CUBIC_CONE_ABELIAN_CONORMAL_SEPARATE_REVIEW.md` の受理済み全入口を精査した。
nonplane rank零・被約既約 quintic なら generic multiplicity二と Frobenius first-jet から
conormal は abelian。その必要条件の complete30×18 source は旧別作者が全行・全kernelを確認済み。
cone vertex の smoothness は要しない。

座標 `(x0,x1,x2,x3,c)` と三 ideal quadric q0,q1,q2 は新 producer note と同順序。
全十八源は、六 base 変数の P と、全十二変数の `P_i,c=W_i=Σqg*t_(4g+i)` に分かれる。
旧 `∂c∧V` convention なら V=-Wである。

この分離が remaining Jacobi を落とさないことも手で確認した。
coefficients が c に依存しないので `[∂c∧V,∂c∧V]=0`。
projective Jacobi の c係数は strict `[P,P]=0` を強制し、残りは
`Euler_base∧[P,V]=0` という V に線形な全条件だけである。
旧全350×12 matrix はこの全条件を保持し、別の vertex二次条件は無い。

全 reduced base support L3∪Q2、actual SL2×scalar の全三非零軌道とその逆は
旧 ORBITS_V2 と独立 receipt に依存する。
初回の誤 semisimple ラベルを用いない。
base零は全 `∂c∧V` が decomposable、rank二以下として別に保持する。

## 2. 全 semisimple 六 vertex 源と全十 entry

semisimple 代表の base は

```
P02=q0, P03=2q1, P13=q2;
W=(k0*q0+2k4*q1+k3*q2, k4*q2, k2*q0/2,
   k1*q0+k2*q1+k5*q2).
```

保存された全12×6 kernel の行4g+iを直接読んで、この全四 entry を独立に復元した。
全六係数を保持している。

actual shear `c'=c+ell*x` は W_i を `W_i+ΣPij*ell_j` に送る。
四つの変化は

```
k0'=k0+ell2, k4'=k4+ell3,
k2'=k2-2ell0, k5'=k5-ell1,
```

であり、k1,k3は不変。ell=(k2/2,k5,-k0,-k4) で同時に消える。
逆は opposite shear、det一、唯一の除算は定数二。
k1/k3零を含め source parameter の反転は無い。

この shear 後、元 coordinate pair 順の全十 entry は正確に

```
(P01,P02,P03,P0c,P12,P13,P1c,P23,P2c,P3c)
=(0,q0,2q1,k3*q2,0,q2,0,0,0,k1*q0).
```

真の completed T38 の全式へ
`p0=2,t2=-1,t3=1,t0=-k3,t1=-k1` を戻すと、三つの higher corrections も含め

```
P_y0,y1=-k3*(y2*y4-y3²),
P_y0,y2=-k1*(y1*y3-y4²),
P_y1,y2=2*(y1*y2-y3*y4),
P_y1,y3=y1*y3-y4²,
P_y2,y4=-y2*y4+y3²,
```

他 entry は零。period Pfaffian open は正確に一。
permutation `(y0,y1,y2,y3,y4)=(c,x0,x3,x2,x1)` で、上の全十 entry と literal に一致する。
raw `P_i,c=-P_c,i` の符号も確認した。
Euler 補正、normal symbol、あるいは zero-magnetic 代表の一致ではない。
よって全 semisimple six-dimensional vertex fiber と全 actual orbit逆の像が T38 に入る。

## 3. 全 nilpotent fiber：constant-rank の適用条件

L3\{0} は既約 smooth な三空間の puncture。
semisimple と非零 nilpotent の二全 actual orbit で、
旧 complete350×12 vertex matrix の rank はともに厳密に六。
旧独立証明書は全kernelの代入零・独立性と非零 minor を確認済みなので、
これは部分matrixだけの上界又は一標本からの一般性推論ではない。
actual coordinate/source scalar transformation は全 linear fiber の同型を与える。

したがって全 L3\{0} 上で constant rank六である。
係数が base に regular な matrix の kernel は、非零六次 minor の開被覆で解ける rank六 vector bundle。
この一般的な局所線形代数を用いており、未証明 flatness を仮定していない。
base が既約なので total space は既約、discriminant非零の逆像は dense open。
その全像が T38 affine cone に入り、tensor map は regular、同 cone は閉であるため、
全 total space の像も T38 に入る。
これにより各 nilpotent fiber の全六係数を一つも残さず戻せる。

base零で kernel rank が跳ぶ十二 fiber には、この bundle を外挿しない。
同補集合は前述の rank二として全対象から除かれている。
zero section は affine cone で保持し、対象非零 tensor を最後に射影化する。

## 4. Q2 全三 fiber の actual 一次因子

Q2 非零代表では base `P12=q0,P13=3q1/2,P23=3q2`、row0は零。
旧 complete12×3 kernel を保存行順に読むと、raw vertex は

```
W=(0,-k2*q0/3+k1*q1/2,2k0*q0/3+k1*q2,k0*q1+k2*q2).
```

従って full 五変数 tensor の row0 は全 k0,k1,k2 上で零。
signed Euler Pfaffian の i≠0 項は row0を含むため零、残りは
`actualD=x0*Pf(P restricted to1,2,3,c)`。
非零なら物理次数一×四という proper factor、零なら rank四 quintic 対象ではない。
全 source係数零、nil/base boundary、actual SL2/scalar逆も因子性を保つ。
したがって Q2 の指定被約既約 quintic 枝は全空。

以上で L3 semisimple、全 nilpotent、Q2、base零と全 zero complements を尽くす。
新成分は零。可約 D の全点分類、cone上のrank二制限、他曲面は別義務。

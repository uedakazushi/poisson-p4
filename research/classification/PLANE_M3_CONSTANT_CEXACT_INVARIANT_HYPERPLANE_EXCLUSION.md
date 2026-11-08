# C-exact 三全源：不変な normal hyperplane による全既約五次排除

2026-11-08。PRIMARY `E_reduced_nonnormal`。producer 純証明完成。
全源必要K1の完全性・全complex壁は既別作者受理済みであり、新M1計算は不要。
本結論の別作者読解は統合担当の工程である。

## 命題

rank二plane、actual Dのnormal位数三、generic normal cubicがdouble+simple、
両primitive normal rootsがconstantで、R=2L+Mの非logC-exact connectionを持つ全三源
`DOUBLE_LINE_E0_EXACT_C/E1_EXACT_C/E2_EXACT_C` は、
全parameter τ・全例外を含め **被約既約quintic actual Dを持たない**。
これらをG上界22/24/23だけの未帰属源から、全点の空性へ強める。
reducible/nonreduced Dの別E、他primitive degrees、全分類の完了は主張しない。

## 1. 全first-symbolと全higher必要source

座標(a,b,c,z,w)、e=0,1,2とする。三つの全P0は

```
Pab=Paz=Pbz=0, Pac=−a², Pbc=2ab,
Paw=−aw, Pbw=2bw,
Pcz=(4−3e)az, Pcw=3τaw−cw, Pzw=(3e−4)zw.
```

c/aのirregular exact connectionを落とさず、K=(3e−4)zwを保持する。
normal z,wの両rowはこのP0でそれぞれz,w倍である。

保存済みfull350×30必要operatorの全kernel K1とfullweight二source P2九を用い、
`P=P0+P1(K1*p)+P2(r)` とする。P1のnormal columnsは正確に

```
P1_i,z = z*w*B_i+w²*C_i,
P1_i,w = z²*D_i+z*w*E_i  (i=a,b,c).
```

元source順ではD_iがrow18,19,20、C_iがrow27,28,29である。
P1のbase/base部分は全zL+wMを保持し、P2九はbase/baseのquadratic(z,w)だけなので
normal rowsを一つも加えない。P2を零又はhomogeneous kernelへ制限する議論ではない。

## 2. 最弱の一般補題

任意のskew quadratic matrix Pについて、あるcoordinate vの全rowがv倍なら、
actual Euler-Pfaffian

`D=Σ_i (−1)^i x_i Pf(P_omit_i)`

はvで割れる。i≠vの各Pfaffianにはv-rowから一つのentryを使い、v因子がある。
i=vの項は外側のmultiplier x_v=vが因子を与える。Poisson条件さえこの補題には不要。

従ってD_i全零なら全Pのw-rowがw倍でDにw因子、
C_i全零なら全Pのz-rowがz倍でDにz因子となる。
必要K1の線形像上の恒等式だけで十分である。

## 3. 全三QQ(τ)源と全八complex壁

independent constant14 source reviewは、generic完全核と全例外多項式の全根を元matrixへ戻している。
三源における全例外は次表で尽くされる。

| e | 保存済み全例外 τ | normal-row factor |
|---:|---|---|
| 0 | 10/3,1/3,−2/3 | 1/3だけw、他はzw |
| 1 | 4/3,1/3,−1/6,−2/3 | 4/3だけz、他はzw |
| 2 | −2/3 | zw |

三つのgeneric完全核では全C_i,D_iが共に零、従ってzw因子である。
e0/τ1/3の追加K1方向はCを加えるが全Dは零。
e1/τ4/3の追加方向はDを加えるが全Cは零。
他の全壁も両三行が全列で零である。
これを全三＋八の11完全核の保存bytesから読み戻し、全行を短いcertificateに保存した。

generic kernelが持つminor/denominatorの全根は上記の全八literal wall入力に戻されている。
新しいrank条件、period壁、generic fiberのflatnessを仮定しない。
全τを以上の有限分割で覆う。各sourceに少なくともz又はwのproper因子があるため命題を得る。
元nonlog connectionのactual正規化・normal frameの全逆は既独立constant-source幾何から戻す。
δ零の補集合をC-exact源に無断で含めず、そちらは既EXACT_B/LOGの別入力のままである。

## 4. 有限batchと、証明に必要な範囲

許可されたONE full11 batchは全9inhomogeneous P2、e1の全homogeneous parameter、
同じ既存rank8/9minor、全350Jacobi、literalEuler/Ore、全Pf/Dを保持して21.440秒/RSS154736KiBで完了した。
全actualDはz²w因子、e0/1の上記critical壁ではw又はz²因子を持ち、純証明と整合した。

しかし本命題にはその新inhomogeneous消去や全残Jacobiの解法さえ不要である。
完全必要K1上のzero-row恒等式と、元P2がnormal columnsを持たないことだけで十分である。
このため追加ambient計算、higher jet、GB、new root/M1 profileには進まない。

`CLASSIFICATION_PLANE_M3_CONSTANT_CEXACT_INVARIANT_ROWS_20261108.json` が主となる短い証明書。
`CLASSIFICATION_PLANE_M3_CONSTANT_CEXACT_ORE_20261108.json` は全許可batchの追加検証である。
旧source/report/manifestは不変、旧五十三FunctionDef ASTとcurrent42契約はpreappendに固定し、
昔のG/E契約を新modeへ自動rebindしていない。

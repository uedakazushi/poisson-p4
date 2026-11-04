# Actual global P22と282-source basis

2026-11-04。FULLU=`D(a*d*F*B*C*J)`上で、凍結済み
`K:R22→R174` のglobal left inverse **P22** を有限線形回路として構成した。
元の22 basis vectors全てについて、**P22*K=I22の484成分をliteral QQ(a,e,d)で検証**。
分母は既知FULLU unitsだけである。minor/GB/新点/searchは用いない。

P22は174個のsource-cofactor coordinatesを入力する線形作用素であり、
2491 scalar nodes、1685 linear nodes、378 core atomsからなる。scalarはcoreだけ、
linear nodesはscalar倍と和だけであり、source変数同士を掛けたり割ったりしない。
従って単なるbasis上の値の表ではなく、R上で定義した実22×174線形写像である。
行列の全core係数は展開せず、回路がそれらを正確に定義する。

## 1. 二つのactual dual inputs

prime quadric basisを

`q'=(q0,q1,q3,q5,ellZ,ellW)=A*q`

とする。Aは既受理det−80の実変換である。lightは
`(x0,x1,x2,x3,x4,Z,W)=(v0,v1,v2,v4,v7,3v0−v1+v6,v3)`。

第一のinputはZ=W=L=x0=0上のactual q0/q5 coefficient pairである。
旧116-nonzero-term証明をarchiveから読み、実quad sourceへの12係数を再照合した。
旧native pairの第2行はC*q5である。第1行のcol2 pivotはp=5a、
remaining5 columnsのSchur rowはexactに `alpha*f`、alpha=−aCとなる。
保存cofactors bに対するliteral `sum b_i*f_i=1` も検査した。

remaining columnsでu_i=b_i/alpha、pivot columnで
`u_p=−sum(u_i*M0_i)/p` と置くと、uはq0を零、C*q5を1へ送る。
従ってq5 dualはC*u、q0 dualは `e_p/p−(M1_p/p)*u`。
**実unscaled q0/q5に対する全4 pairing=I2** をexactに検証した。
未知係数を求めるlinear solveや新しいBézout searchはしていない。

第二のinputはell=0四quadの既存21×4 right inverseである。
元データの第3core変数はP=ad。保存matrixにはa-clear shifts `(0,2,2,1)` がある。
P→adへ戻し、実unscaled四quadへのclearing multipliers m_iを照合した。
old right inverse L4に対して

`V4[:,i]=L4[:,i]*a^shift_i*m_i/g`, `g=−a^4*d²*J`

とする。使用したm_iは `(1,J²/25,J²/25,J/5)`。
actual ell=0 coordinate mapの全84係数を保存804-term sourceと照合しており、
P/a座標・row unit・変数順を推測で同一視していない。

## 2. 新22本目を取り出すactual extractor

Z=W=0へのrestrictionをbarとする。

`L_x1=−40C`, `ell_x1=−16CJ/243`,
`Delta=ellbar−(2J/1215)L`,
`g0=3Delta_x2−(7/2)Delta_x4=−80ad²/81`。

これらはFULLU unitsである。定係数微分作用素

```
D0=3∂x2−(7/2)∂x4
   −((3L_x2−(7/2)L_x4)/L_x1)∂x1
```

はLとx0を殺し、ellbarをg0へ送る。
quadratic pに対し `beta(p)=D0(∂x1 p)/(g0*L_x1)` と置く。
§1のq0/q5 dualを使って、pのL=x0=0 restrictionに含まれるq0/q5成分を引いてから
betaを適用する作用素lambdaを構成した。

実際に

`lambda(q0bar,q1bar,q3bar,q5bar)=(0,0,0,0)`、
`lambda(ellbar*L)=1`

の全5恒等式を検証した。prime cofactor c'のq3成分c'_2に対し

`n=−lambda(c'_2|Z=W=0)`

が、新22本目の係数を取り出す。旧21本を全て殺し、新22本目を1へ送ることを、
**元Kの全22列へ実際に適用して検証**した。StageNewは303 scalar/426 linear nodes。
これは他のdual段階から独立して完了したcheckpointとして保存する。

## 3. global six-quad dualを構成する

ell_x1がunitなので、lightをellと
`w=(x0,x2,x3,x4,Z,W)` に取り替えられる。§1のV4とell=0 restrictionから、
最初の四quadをpairingでI4へ、ellZ/ellWを零へ送る四つのfunctionals phi_iを得る。

psi_Z/psi_Wを、この座標でのellZ/ellW coefficient extractionとする。
残る二つを

`mu_Z=psi_Z−sum_i psi_Z(q'_i)*phi_i`、同様にmu_W

とすれば、これらは最初の四quadを殺し、ellZ/ellWにI2を与える。
このglobal6 dualについて、**actual六quadに対する全36 pairing=I6** を
literalに検証した。pointwise rank6だけからdualの存在を推論したものではない。

functionalsは6×28のscalar-root tableとして保存した。表の行はdual、列はquadratic
monomialsである。coefficient matrixをQ6×28とすれば、pairingは
`Q*(dual table)^t=I6`。各scalar atomの分母をnative側でも既知unitsで検査した。

## 4. 残り21係数とorientation

original cofactor cは `c'=A^-t*c` でprime cofactorへ移す。
新係数nとそのactual cofactor vectorを引いた後、最初の四cofactorにはKoszul部分
だけが残る。§3のdualを適用して、pair(i,j), i<j, i<4の14 Koszul係数を取り出す。
それらを引いた残りは

`(c'_4,c'_5)=(W*f,−Z*f)`, `f=l+a45*ell`。

選んだ6 light multiplesはx1を含まないので、Wx1 coefficientをell_x1で割るとa45、
残る6 coefficient extractionでlが得られる。**Wで変数除算はせず**、全てquadratic
coefficient extractionとして定義している。

Koszul係数を元basisへ戻すときの明示式は

`aOrig_kl=sum_(i<j) aPrime_ij*(A_ik*A_jl−A_il*A_jk)`。

すなわちskew matrixでは `S_original=A^t*S_prime*A`。
最後の6 linear multiplesと新22係数を加えた22 outputsの順序は、凍結Kの元順序に
厳密に一致する。**全22 outputs×全22 source basis vectors** のexact検査で
orientationも確認した。

schema上の注意: `global-P22.json` の `inverse_prime_cofactor_transform` は
A^-1のrow tableを保存する。実際のcofactor操作はそのtranspose、
`c'_j=sum_i (A^-1)_ij*c_i` である。field名だけを見てA^-1*cと読まないこと。

## 5. 同じcokernelを持つactual global282 basis

受理済み固定M304のsourceを `R174⊕R22⊕R108` と分け、map向きで
`M304=[A0 B0 C0]` と書く。B0は残る130 sourceの最初22個である。
PK=I22、A0K=0なので

`N282=[A0+B0*P22,C0]:R282→R239`

は同じ実cokernelを持つglobal表示となる。
具体的に、新source(z174,h108)を旧sourceへ `(z,Pz,h)` と持ち上げる。
任意の旧source(x,u,h)は

`k=Px−u`, `z=x−K*k`

により `kernel(Kk,0,0)+(z,Pz,h)` と一意に分解される。
従ってこれはglobalなsource basis変換であり、generic pointだけの操作ではない。
`global282-basis.json` は実304/319 source indicesとこの係数回路を結び付ける。

**282個の旧生成行を選んだsubmatrixとは主張しない。**
新174行の一部はquadratic sourceと22個のcubic sourceのP係数による和であり、
残り108行はそのまま保つ。全quartic239係数を展開せず、実global282 mapの
coefficient oracleを有限回路で定義した結果である。

このP22は既知22方向をsplitしたもので、新しいkernel方向数は増えない。
generic kernel rank43という残るgate、欠損曲線・点、fullquartic239、C122rank64は
なお未解決である。

## 6. 実行履歴と検証範囲

元120秒/2GB/32MiB枠の一回のconstructive jobで、StageNewは1.371秒で通過した。
続く保存right inverseのreaderが「1 polynomial=1 line」と仮定し、折返し行を含む
既存gzipで停止した。算術identityの失敗ではなく、失敗source/hashと完了checkpointを
保存した。readerをwhitespace streamへ修正したresumeでは、再構成prefixと保存済み
checked circuitのliteral一致を確認し、StageNew native検査を繰り返さず続行した。

resumeは3.581秒、full native45 scalar checks/484 matrix entriesは0.967秒。
旧20pivot/66minor jobの再実行・候補拡張・新しいBézout探索はしていない。
分母は全てFULLUの既知unit、全source/frame identitiesを保存している。

再現は新output先で

```
python code/s61_nonsquare_kernel22_extractor_20261104_v63.py --output <new-directory>
```

とする。`--resume-new22` は保存済みcheckpointのprefix一致を要求する再開用である。
公開packetは失敗・再開双方を保持し、成功したactual P22を独立監査に渡す。

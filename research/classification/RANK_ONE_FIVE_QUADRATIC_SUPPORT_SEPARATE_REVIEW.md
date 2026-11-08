# rank一long残五源：真正二次余核支持の独立監査

2026-11-08、PRIMARY E_multiple_hyperplane。**独立Fraction算術と形式切片の数学を受理。B1/B3/V1は25次元全成分。M1/M11のactual higher支持は未解決。smooth性・局所成分の唯一性は主張しない。**

## 1. 同じ五点と全actual方向

入力は凍結された元九点のうちM1/M11/B1/B3/V1だけである。前段の独立九点reviewが、元g,U,H,Fからの全150係数、全350 JacobiとPfaffian、全D／Euler／GL／scalar行列、既minor及びcomplete kernelを照合済みである。本検証では同じDとMのbytesをbindし、新しいD、点、normal profile又は高次jetを生成しない。

元六括弧から、全二次F六係数と次のcanonical源をparameter付きで再構成した。

| source | 全U | 全H |
| --- | --- | --- |
| M1 | −z∂z+(v0+vz z+vw w)∂w | t z |
| M11 | λz∂z+(−1−λ)w∂w | t zw |
| B1 | −z∂z−w∂w/2 | A z+B w² |
| B3 | −z∂z+w∂w | A z+B z²w |
| V1 | −z∂z | z(Aw²+Bw+C) |

全parameterで(U+1)H=0、全strict10／projective350 Jacobiが恒等零であることを独立sparse polynomial演算で確認した。V1のB方向、M1の全v方向とM11の全λ方向を保持する。元指定parameterへ評価した点は旧150-vectorと完全に一致する。

各全parameter方向を実微分して旧M51へ連結したFを再構成した。全DF=0、保存F50 minorの正確なnonzero determinant、全FがF50のspanへ入る恒等式を確認し、rank F=50を得た。単に群作用方向の数を50と見積もっていない。

## 2. 真の50cutと全150方向の直和

旧全D kernel Kのfree-coordinate identityから、指定されたkernel列を取り、true normal Nを再構成した。旧Dのindependent coordinate列Uと、F50／Nの本数は50+d+r=150である。

明示された50×150のcut行列について

```
cut F50=I50, cut N=0, cut U=0
```

を全entryで検証した。NはKの相異なる列なので独立、DUは旧nonzero minorでinjectiveである。従って[F50,N,U]は全150係数方向の直和となる。これは固定divisor源内だけの切片でも、都合のよいnormal部分空間でもない。

## 3. 全350二次式と全Dを消す余核

N_iの全Jacobi J(N_i)と、各i<jに対するJ(N_i+N_j)−J(N_i)−J(N_j)から、d(d+1)/2本の全350係数Qを独立に再構成した。旧selected equations I／coordinate columns Jのminor Aを使う二次補正Tについて

```
A T+Q[I,:]=0, C=Q+D[:,J]T
```

を検証した。Aは既知のnonzero定数なのでTは唯一であり、saved inverseだけを信用していない。Cの全entry、selected rowsの零、全350行のgenerator basis分解を確認した。

各generatorを取り出すfunctionalは、旧selected rowsへの係数Bを用いたφ=e_i−B e_Iである。独立に

```
D[i,:]−B D[I,:]=0  （全150列）,
Q[i,:]−B Q[I,:]=C[i,:]
```

を照合した。全Dの余核であり、chosen sourceの一部又は未解決の補正を捨てたquotientではない。

## 4. 固定degree d+1の証明書

全二次generator×全degree d−1 monomialからMacaulay行列を再構成し、元全entryとrankを照合した。成功時はdegree d+1の**全monomial**のmembership lift、失敗時は非零left annihilatorを元行列へ代入した。degree探索、適応的な延長又は高次項の計算はしていない。

| source | normal d | 固定degree | Macaulay rank／全monomial数 | 受理範囲 |
| --- | ---: | ---: | ---: | --- |
| M1 | 3 | 4 | 10/15 | quadratic radicalはplane t2=0、actual higher支持は未解決 |
| M11 | 2 | 3 | 2/4 | nondegenerate binary quadraticの相異なる二複素line、actual higher支持は未解決 |
| B1 | 2 | 3 | 4/4 | m³が真のquadratic initial idealへ入る |
| B3 | 1 | 2 | 1/1 | 真の式−t0²/2 |
| V1 | 1 | 2 | 1/1 | 真の式−12t0² |

M1の全generatorはt2を因子に持ち、残る三linear formが独立であることを別に確認した。従ってI2=t2(t0,t1,t2)。M11の全idealは `2(231t0²−110t0t1−70t1²)/1155` で、discriminant非零を確認した。

homogeneous quadratic idealがm-primaryなら、無限体上でdegree二generatorの一般線形結合からd本のregular sequenceを選べる。そのcomplete intersectionのsocle degreeはdなのでm^(d+1)は元idealへ入る。従って今回の固定degree testの失敗を単に「次数不足」と解釈しない。一方、初期二次idealが正dimensionでも、元Poisson切片にその枝が実在するとは結論しない。

ONE independent実行は3.389785秒、maximum RSS 18012KiB、全五caseのcertificateがPASS。標準ライブラリFractionと凍結済み独立engineのみを使い、producer／SymPy／CASはimportしていない。

## 5. 成分支持への形式的な矢印

実cutで残る100方向を(U,N)で書く。selected r equationsのU微分はAなので形式陰関数定理が適用でき、U=g(t)を消去出来る。DN=0よりgは二次以上で始まり、全Dを殺すφを適用した残余方程式の二次初期式は上記Cとなる。

B1/B3/V1ではそのinitial idealがm-primaryなので、真の形式切片のKrull dimensionは零。元150係数Jacobi局所環Rへ50本の実linear cutsを加えたquotientが零次元だから、Krullのdimension不等式によりdim R≤50。一方、別に受理済みのwhole既約25次元像のconeとEuler24方向は50次元の下界を与える。従ってその閉既約像は全Poisson schemeの25次元成分である。

接次元はB1/B3/V1で依然27/26/26であり、滑らかさは従わない。また `y(y−x²)=0` のx方向cutがy²を与える例のように、切片の原点支持だけでは接する別成分の不存在も従わない。局所唯一性、局所積分解、全横断代数又はnil指数は本結論に含めない。

## 6. 全source・相違・残余

全B1/B3/V1 sourceと係数境界、全short counterpartの帰属は、既whole source閉包とactual short-to-long曲線から従う。全prior／相互相違は独立pure入口のmarked factor／rank／incidenceを使う。新しいpointのfactorだけで比較しない。generic軌道23/22/21と商次元2/3/4は、前段の全源vertical Poisson場と同じ旧Mの下界を合わせた別の純補題を使う。

元有限十七sourceのうち十五源の全帰属が閉じ、残るのはM1/M11。二つのquadratic coneを元解集合の余剰枝として採録しない。全E又は全分類の完成を主張せず、台帳・Gitは本reviewでは変更しない。

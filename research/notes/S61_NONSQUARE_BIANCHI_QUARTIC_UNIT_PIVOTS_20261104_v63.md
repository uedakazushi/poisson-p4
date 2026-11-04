# Full Bianchi の11本の新しい quartic kernel directions

2026-11-04。ここで扱うのは元135式の full Schouten 恒等式から来る
**319生成行の間の関係**である。周囲239列の新しい行列、点、field、GBは
計算しない。既証の generic239 を再計算しない。

## 1. 結果と開集合

`S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63` の直接weight4 piecesは
15本である。degree4由来5本、degree5由来10本、degree6由来は零である。
33 genuine quartic columnsへのconstant coefficient blockは15×33、非零行11、
全体88 core termsである。

このblockに、既知unitだけを割る行操作を一回適用した。**rank11** のpivotが
得られ、残り4行はliteralに零となる。選んだquartic columnsでI11を持つよう
正規化された11本のactual319行関係を、有限DAGで保存した。

新しい11本は `a*F*B*C!=0` 全域で独立である。ここで

```
F=500a+90e²+3e−24,
B=(50a+18e²−9e−2)d+60e−40,
C=100a+18e²−9e−2.
```

旧quadratic22関係はquartic列が零であり、その独立性の既証開は
`FULLU: a*d*F*B*C*J!=0` である。従って両者を合わせて **FULLU上33本** が
点別に独立となる。80本のsource-kernel gateには、なお **47本不足** する。
旧22本をaFBC全域で独立と読み替えない。

## 2. 元135式とのliteral binding

公開operatorの全135 q/t斉次化式をq_fixed=1でlocal coordinatesへ戻し、
lazy oracleの全J/B atomsとliteralに照合した。零weight nuisance19/23を除く
各項は保存されたequation weightと一致する。nuisance19/23はTとUの両方で
零であり、全implicit graph上で恒等的に零である。

15本の対応する **full135** operator rowsについて、各係数と元斉次式を
掛け、cubic q/t polynomialの全係数が零であることを直接検査した。
続いて実際のweight4 pieceが元行のweight分解と一致することを照合した。
上流の175本全検査とtG6恒等式もhashで固定するが、tをunitにしてG6を
除去しない。tG6はweight6であり、今回の11方向に数えない。

## 3. T34/Tnと行単位

oracleの変数をy、旧quadratic/cubic sourceの変数をvとすると

`y=D*v, D=(a,a,a,1,a,B/50,a,a)`。

quarticは全33本について `Qhat_g(v)=F_g(D*v)/D_H²` とする。これにより
旧H² pivotの規約を保つ。block上の共通因子D_H²を除いた後、選択11列が
ちょうどI11になる。

既受理adapterの6 quadratic unitsと16 retained cubic unitsを使う。35行削除に
必要な残り5 genuine cubicの単位も固定した。元63行indexで
`48,55,57,59,61`、physical135行で `59,113,116,122,132` である。それぞれ

```
−81a³/200,
−567a/(10F),
−8100/(F*B),
140400/(F*B),
−1455300/(F*B).
```

whole-row covarianceの根拠は既受理のfull63 physical reconstructionである。
native graphはTn=T34D、Un=ud*Uを用い、抽出には共通unit ud³が掛かる。
保存した短いexact anchorsはこの既証whole-row scalarを決定するために使う。
anchorだけから未知の行恒等式を推測したものではない。

full63 archiveのheaderは`63 91 57`である。63/91は実際の行列サイズであり、
57はwriterのtarget metadataである。保存されたSchur parserは第3値を読み捨てる。
その57をrankの証拠に用いない。archive全体のfull member SHAも照合した。

## 4. 15関係の有限輸送

各weight4関係のmultiplierはlocal coordinateにlinearである。選択34式は
implicit graph上で零なので除く。残るweight2式のmultiplierにはweight2部分
`E+Z2`、weight3式にはweight1部分E、weight4式にはconstantだけが寄与する。
weight0 nuisanceは零である。正weightゆえ、Z3/Z4以後の項はこのmultiplierに
追加し得ない。元式のweight≤4 residualは既存oracleのexact finite graphである。

各係数をy=Dvへ移し、上のknown-unit row scalesで

* 6 quadrics ×29 weight2 monomials =174 rows;
* 21 normalized genuine cubics ×7 light variables =147 rows;
* 33 actual Qhat quartics

という354生成行へ移した。15本すべての係数ベクトルを保存した。既受理の
35個のunit恒等式を順に引き、deleted coefficientを零とし、全319 retained行へ
移した。35本はquartic columnを持たないので、11個のquartic pivotは保たれる。

この輸送は、元のliteral polynomial identityへの**有限代入・行単位変更・
既受理の行恒等式**である。大きな239列の全core coefficientを展開して
再加算する方式ではない。certificateは元piece番号、354係数DAG、35削除の
各step、319係数DAG、11個のunit正規化済み関係を保持する。初期1373 nodesから
最終10710 nodes、core atoms4228である。

## 5. 実行記録と限界

最初の実行は6.371秒で135式binding、15元恒等式、11pivotと残余零を通過した。
続くarchive headerの第3値を63と誤記したassertで停止した。失敗log・当時の
source hash・checkpointを保存した。許可された一回のresumeでは既存pivotを
読み、source検査/pivot計算を再実行しなかった。headerをwriter/parserで確認し、
輸送と35削除を32.653秒で完了した。合計39.024秒、2GB/60秒枠内である。

再現は新output先に

```
python code/s61_nonsquare_bianchi_quartic_unit_pivots_20261104_v63.py --output <new-directory>
```

とする。公開reportは実行済みのcheckpoint/resumeを記録する。新規実行なら
同じsource検査とpivotから始める。今回さらに再実行していない。

これは11個の**source関係の独立性**の証明である。full80 kernel、239点別下界、
actualC122のrank64、full quartic全射性、affine exhaustionを結論しない。
残余4本のquartic blockが零でも、他の列で新たな独立方向を持たないとは
主張しない。次はこの残余または異なるweightのsource関係を調べる必要がある。

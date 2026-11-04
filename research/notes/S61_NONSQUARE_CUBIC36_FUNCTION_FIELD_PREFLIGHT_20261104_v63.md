# 全 cubic defect の generic-a 有限 module 方針: 読出し preflight

2026-11-04。対象は frozen V6 の `after21`、SHA `c65b3d80b2207e75ba794821bb41ebc3ffc864959dab27307af77cdeab9d8fd1`。matrix の向きは **column map R^70→R^36**。QQ(a,e,d) 上の coker は零である。旧 raw42×70 は一般 coker に free rank6 を持つため、whole coker の Ann を使う対象ではない。

元63行の自由 row module から alternative ad57 行基底を選んだ段階では、triangular minor が a³d/10 なので d≠0 が必要。その後の21 unit Schur pivots/row cleaning が使う polynomial factors は a,F,B,C だけであり、J は使わない。保存36×70の entry は QQ[a,e,d] の polynomial で、parameter denominator はなく、scalar rational denominator だけである。全 polynomial record・degree・termcount・rational coefficient の32003 integrality を一回読み出す。旧 prefix 全体の scalar32003-unit 性を新たに再計算するものではない。

ここで F=500a+90e²+3e−24、B=(50a+18e²−9e−2)d+60e−40、C=100a+18e²−9e−2。比較対象の実 defect open は a*d*F*B*C≠0。F32003(a) 上で a は既に単元なので、次回計算が別途許可される場合の適切な presentation は

```
ring R=(32003,a),(e,d,z),dp;
// M は全36×70の source polynomial を変数順(a,e,d)から正確に移す。
poly U=d*(500*a+90*e^2+3*e-24)
 *(50*a*d+18*d*e^2-9*d*e-2*d+60*e-40)
 *(100*a+18*e^2-9*e-2);
// N=[M | (z*U-1)*I36], size36×106。
// module I=N; module G=std(I);
```

localizer relation を一列だけ付けるのは誤りであり、36個の target component 全部へ付加する。`z*a*F*B*C-1` だけの別 presentation は saved coker については正しいが d=0 を含み、alternative ad basis の元 cubic rank との比較には使えない。この余分な locus を避けるため上の d factor を明記した。

許可された場合の最初の成果候補は G 自体、source入力 SHA、term order、全 leading module、標準 module monomial basis、および e,d,z の coker multiplication matrices である。零次元ならその有限基底を保存し、各 multiplication の normal form と入力 generator の零 reduction を検査する。`primdec.lib` の既存 `Ann(module)` は whole coker の annihilator を返すが、内部に quotient1 計算を含むため最初の bounded GB と同時に無制限に実行しない。有限 module multiplication matrices からの characteristic polynomials は source support の必要方程式を得る別の有限手段である。

QQ で既知の codim≥2 により QQ(a) 上の実 support は零次元以下だが、mod32003 で同じ次元が保たれるとは仮定しない。GF の module が unit/有限になっても QQ へ直接昇格しない。QQ の exact certificate、または必要方程式から得た finite Z_(p)(a)-model と全 localizer の no-escape を含む別の proper/Nakayama 証明が要る。finite-field coefficient の a-denominators と leading coefficients は全て記録し、その消滅する a-special cases を generic-a probe が覆うと主張しない。

この packet は source format/size/provenance の **read-only preflightだけ**。GB、field arithmetic job、新しい点、元matrix再生成は一切実行しない。既存 Singular wrapper/runtime を再利用する計画であり、新規 install も行わない。

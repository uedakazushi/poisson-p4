# A25全G証明の座標符号とresidue pole位置の訂正

2026-11-07、producer plane_collision_review。
対象は凍結済 `CUBIC_SCROLL_A25_WHOLE_NONDOMINANCE.md`、
SHA256 `f08455bf52022f7cb3932f9428801c204d969189bd3903ed66a6a6f2c6c999a5`。
原bytesと実行記録は変更しない。結論、全source cover、parameter数又はimage上界に変更はない。
新CASは行っていない。

## 1. §6のincidence解の符号

全incidenceは

`a x3−b x2+c x0=0, a x4−b x3+c x1=0`。

a=x2=1のdense chartで、X=x2、Y1=x3、Z=x4とすれば、正しい式は

```
Y1=b−c x0,
Z=b²−bc x0−c x1,
XZ−Y1²=c(b x0−x1−c x0²).
```

原§6の`Y1=b+c x0`と、それに続く二つの式の対応する符号は誤記であった。
このchartはc=0でY1=b,Z=b²へcollapseし、conic equationのc-orderはgenericに正確に一である。
これが§5のc|Ω排除と§6のinverse-residue恒等式で用いた全情報である。
正しい式でもこれらは同じなので、

`Ω_L W=k*k1*k2`（node）、`Ω_L W=k*k0²`（cusp）、k非零

は変更されない。旧normal-Lie-bracket計算を結論の根拠にはしていない。

## 2. §8のresidue pole位置の適用範囲

Ωが二linesの和の場合、base inverse formの各component上のresidueは
Lとの交点に加え、他のΩ componentとの交点にもpoleを持ち得る。
従って原§8の「polesはL∩Ωだけ」は、**Ωがsmooth irreducible conicと判明した後**に限る。
reducible Ωの排除で用いるのは、η_Bが各componentのgeneric pointで非零であることと、
少なくとも一componentでrow ratio f0がconstantなのでdlog f0=0となることだけ。
pole位置の誤った限定をその排除に用いていない。

smooth irreducible Ωでは、residueはΩ¹_Ω(L|Ω)の非零sectionであり、
pole位置は正確にL∩Ω（intersection multiplicitiesを含む）である。
row-ratio同型f0:Ω→P1のzero/pole二点をそこへ制限する後半の議論は、この正しい範囲で成立する。

## 3. 結論と独立確認

nodeの全22/31分類、Ω reducedを先に証明する順序、31排除、cuspの全五block、
真のpresentation fibers、全PGL像≤24は変わらない。
第二rowの非primitive性、全full-P2 fibers、Rの全非被約境界を引き続き保持する。
第2項は `CUBIC_SCROLL_A25_BOUNDARY_SEPARATE_REVIEW.md` でも明示されたscope補正である。
第1項はproducerの後続再読による符号訂正であり、独立確認の受理と混同しない。

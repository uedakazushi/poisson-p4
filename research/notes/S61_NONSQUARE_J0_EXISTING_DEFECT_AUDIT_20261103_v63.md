# J=0 の有限性を既存の全三次欠損定理から読めるか

## 結論

**三core全開で full cubic rank57 欠損が次元1以下、という定理は既に証明・独立監査されている。しかし、既存の heavy 欠損曲線除去と合わせても、J=5ad−6e−4=0 上の full cubic 欠損の有限性はまだ従わない。**

今回の作業は既存ノート・報告・handoff の読解と適用範囲の監査であり、大行列solverや既存証明書の再実行は行っていない。J0用18×52 matrixの追加消去も行っていない。

## 1. 既存の全三次定理は product の定理ではない

`S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102.md` §1–4 と同日の独立 REVIEW は

\[
 U_0=\{aFBC\ne0\}\subset\mathbf A^3_{a,e,d},\quad
 \Sigma=\{\operatorname{rank}N<57\}\subset U_0,
 \qquad \dim\Sigma\le1
\]

を証明する。Nは **63×91の全重み三次行列** で、六二次式のlight一次倍42本と、元135式由来の三次式21本をすべて含む。指定単位のみの21 pivotにより42×70のMを得て、全U0で `rank N=21+rank M` である。従ってこれは product41だけ、あるいはheavy7だけの次元上界ではない。

主report `reports/S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102.json` の結論は `cubic_rank_drop_codimension_lower=2`、`cubic_rank_drop_dimension_upper=1` である。証明は次の三段を明示的に組み合わせている。

- MODULAR：有限体上で実際の五つの36 minorの共通未知因子を排除する。最初の二つのminorのfactor bagに共通する三つの因子H1,H2,H3は、**全欠損集合の曲線成分のリストではない**。残る三minorがそれらの各超曲面上で非零になることを使う。
- HORIZONTAL：指定minor PのQQ上のd次数79、p-goodな先頭係数、既知B因子の正確な重複度12により、未知のd依存共通因子の還元が非定数のまま残ることを保証する。
- VERTICAL：P(a,e,3)のe次数213、lex先頭単項式a^39 e^213のp単元係数、a,F,Cの正確な重複度37,10,26により、残るa,e因子を排除する。

有限体のgcdだけからQQへ推測した結果ではない。独立REVIEWはこれらの同一行列・先頭係数・Gaussの補題の接続を監査し、垂直証明を別の作業先で再実行している。今回その大計算は再実行していない。

## 2. 現在のfull cubic matrixとの同一性

上の主reportが使うSchur matrix SHA256は

`0946ea9d6b5ce952d12c6ac7095c91b4f7ef60f20c05fd2b4f558f5c00c0a0a2`。

Cloud側の `S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json` も同じSHAを持つ。さらに同packetは、旧archive member `native_data/hyperplane-full-cubic-schur36-input.txt` と現在のSchurの **全2940 entryが同順序で一致** することを記録している。旧member SHA `b24dc826…73107` が異なるのは、旧headerが目標rank36を、新headerが既消去pivot数21とmetadataを記すためであり、別の行列を扱ったという意味ではない。

V6の27×61 packetは、この同じ全三次表示をalternative ad57-row basisへ移した後、さらに指定単位pivotを行ったものである。ただしその30 pivotのうち後半はJを使う。V6 READMEとSTATUSは明記している通り、**aFBCdJ≠0上のrank-preserving reduction** であり、all-point rank57 theoremではなく、J=0へそのまま代入できない。J0-safeなのは最初の21 pivotまでである。

## 3. Heavy欠損Yの除去はΣの全曲線の分類ではない

claim ID `S61_NONSQUARE_HEAVY_DEFECT_REMOVAL` のproof sourceは、`S61_BINARY_NONSQUARE_HEAVY_CURVE_CUBIC_20261103_v62.md` と同QUARTIC、およびCUBIC_REVIEWである。claim自体のscopeも “Specified heavy-defect locus, not all other cubic or quartic defect curves.” と限定されている。

15×7 heavy-light matrixの欠損集合をYとすると、既存結果は

\[
 \dim Y\le1,\qquad \dim(Y\cap\Sigma)\le0
\]

である。さらにYの各曲線一般点で四次全射と実横断長65が成立し、同じH11入射、nuisance収縮、次元評価からYに由来するambient bad部分を一般β上で除去する。

従って **Σの曲線成分はYに含まれない**。もし含まれれば `dim(Y∩Σ)≤0` と矛盾する。未処理のfull cubic欠損曲線は、むしろ一般点でheavy rank7である。Yのexplicit annihilator25式が、Σのannihilatorであるという包含・等号は証明されておらず、ノートも別の三次・四次欠損曲線の可能性を明示している。

したがって「heavy曲線の全成分でJが恒等零でない」と確認しても、Σの曲線成分がJ=0内に無いことの代用にはならない。Y上の300 minorやその固有模型による移送を、Yの外のΣへ適用することもできない。

H11相対自由65の定理も、適格な三次・四次minorが非零の開で成立する局所結果である。そこに現れる「不適格補集合のcore次元1以下」は、同じ残存曲線の有限性を既に証明したという意味ではない。係数chart U0と、長65の適格開を同じ記号Uで混同しないことが必要である。

## 4. J=0の有限性に本当に必要な追加命題

既にΣは高々曲線なので、求める命題は

\[
 \text{Σのどの1次元既約成分Γについても }J|_\Gamma\not\equiv0
\]

である。これが証明されれば、各ΓとJ=0の交叉は零次元であり、有限個の成分と既存の孤立点を合わせてJ0の欠損は有限となる。

**J0の一般点でrank57であることだけでは足りない。** それはdimΣ≤1とdimJ0=2から既に従う。今回のexplicitなQQ点(1,2)の非零minorは有用な整合検証であるが、この既存一般性を強める有限性証明ではない。一般に、A3内の曲線V(x,y)は次元1であり、超曲面V(x)の一般点を避けつつその超曲面に全体が含まれる。

また旧B-free horizontalのQQ(a,e)[d]上の全射は、Σの曲線が(a,e)平面を優越しないことに依存する。それをJ0へ戻しても、一つの未知曲線p(a,e)=0上でJ0と一致する可能性は残る。QQ(a)[e]とQQ(e)[a]に必要な、一座標だけのannihilator証明とは区別する。

既存資料からΓの全リストや、Σ∩V(J)の次元0証明書は見つからなかった。次の対象は、現在の18×52 cokernelの支持、または同じものを記述する実際のfull minor idealに **Jを加えた支持** である。既存の五minorをJ0へ制限して追加のgcd・良い還元を証明する方法は候補になるが、既存の三変数gcd=1からその制限gcd=1を自動的に推論してはならない。

## 参照

- 全cubic欠損次元上界：`research/notes/S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102.md` と同REVIEW、主reportおよびMODULAR/HORIZONTAL/VERTICAL各report。
- 現在のmatrixへの接続：`research/notes/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62.md`、その`report.json`と`source-report.json`、`research/handoff/20261103_v62/v6_cubic27_coker/{README.md,STATUS.json}`。
- Heavy Yの範囲：`research/notes/S61_BINARY_NONSQUARE_HEAVY_CURVE_{CUBIC,QUARTIC}_20261103_v62.md`、CUBIC_REVIEW、claim `S61_NONSQUARE_HEAVY_DEFECT_REMOVAL`。
- 残存全体の完了条件：`research/notes/S61_NONSQUARE_MARKED_REMAINING_20261103.md` §5。

本監査は既存定理の適用範囲を整理したものであり、新しい全三次欠損・四次・ambient除去の定理を追加しない。

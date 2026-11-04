# 残るbad curve：明示sourceと欠ける有限certificate

2026-11-04。新しい固定cut transferにより、actual full319×239四次欠損はU0={aFBC!=0}で次元≤1。さらに六quad定数rank6が全U0で証明され、旧Z2は空である。したがって必要coverは

`Z_any65 union Z3`

へ縮む。ここでZ_any65は適格cutの存在する旧開U_anyの補集合、Z3はfull cubic rank57の欠損。これは集合のcontainmentで、Z3がfullquartic欠損を尽くすという主張ではない。

## 1. Z3の明示的な五determinant-circuit cover

CUBIC_MODULARのsafe21 Schurは42×70、source SHA `0946ea9d6b5ce952d12c6ac7095c91b4f7ef60f20c05fd2b4f558f5c00c0a0a2`。最初の36minorと隣接minorはchoices[0,1]、残り三minorはfull_rank_witnessesのrow/column indicesで完全に指定される。五つとも同じQQ sourceの実際のminorなので、Z3はその共同零点に入る。

旧horizontal/vertical good-leading/Gauss argumentは、この五つのminorに共通height1 factorがないことをaFBC局所化上で証明する。従って五つの共同零点自体が必要な次元≤1のQQ coverとなる。単独のmodp因子分解をQQ factorizationへ読み替えたものではない。

ただし、このcoverは小さい式にまだなっていない。

- first minorのQQ d次数は79、B^12によるglobal整除が確認済み。Bがunitの開では残りのd次数は67。d=3ではQQ e次数213とlex-leading p-unitも確認済み。
- first全次数310・個別次数(a,e,d)=(161,213,79)、adjacent全次数307・個別次数(160,214,77)という列挙は**GF32003 factorbagsの値**。adjacentのQQ d次数77は別に証明されていない。
- 完全な二つのGF factorbagは旧handoff archive内に保存される。一方に123396項・totaldegree144、他方に104388項・totaldegree135のterminal factorがある。全pairwise factor intersectionsを小さな有限algebraと仮定できない。
- 二bag共通の未知GF hypersurface factorsは21,92,6項、個別次数は(2,4,2),(6,6,5),(1,2,1)。追加三minorの非消滅はこれらが全minor共通因子でないことだけを示す。これら三因子の和が全bad curveを覆うわけではない。

raw42×70はgeneric rank36なので、そのcolumn coker全体のannihilatorはfree rank6のため零。ann計算の候補には、既存global ad57 basisから21safe pivotsを施した**36×70**を使う必要がある。frozen V6/J0 packetに元入力とunit21検証があり、現在のfileは `work/nonsquare-v62/cubic27-v6/after21`、SHA `c65b3d80b2207e75ba794821bb41ebc3ffc864959dab27307af77cdeab9d8fd1`、約13MB。Jはまだ逆にしておらず、ad basisを使うのでd!=0を仮定する。後の27×61巨大cacheと混同しない。

この36×70を既存primeのF32003(a)[e,d]で一度profileする案は、degree/lengthの候補を得るための**未実行案**に留まる。GF profileからQQ finite modelへ自動移送してはならない。

## 2. Z_any65には明示QQ ann polynomialがまだない

BFREE_HORIZONTALはQQ(a,e)[d]上の全横断代数free65を証明する。証明はdegree259のmonic minorとfinite-local Nakayamaを用い、必要なscalar denominator f(a,e)の存在を与えるが、そのQQ式やunit membership liftは保存されていない。

旧TRANSVERSE_CODIMENSIONのa-degree575 polynomialはd=3へのsliceであり、e-only degree942 polynomialはa=1,d=3のsliceである。それらは、不適格**divisor因子**がdに依存しない、またはeだけに依存するという仮定の下で使われる。一般のbad curve上でこれらslice polynomialsが零になるとは言えない。必要なcurve equationsとして使用しない。

既知D0面自身のgeneric pointには元raw quartic非零minor、C面にはalternate cutがあるが、これも全残余curveの明示方程式を供給しない。

## 3. 既にproper modelのある例外との区別

heavy7 defect Yには小さいann/finitefree4 modelとleading-exception fibersがあり、その全curveの一般点でfull cubic/quarticを検証済み。しかしYは全Z3や全Z_any65のcoverではない。

finite gcd-exception Sigmaにはmonic e4/P28、全localizer-unit、28次元special algebraのproper modelがあり、raw24およびrawquad78を全Sigmaで証明済み。Sigmaもここで残るfullquartic bad curvesのcoverではない。

## 4. 次に必要な最小のcertificate

実quarticをcurve上で有限QQ(parameter) algebraへ移すには、少なくとも次の実dataが必要である。

1. fullquartic coker supportを含む**明示的なQQ必要curve cover**。source minor/ann membershipを保存し、単なる特殊sliceやmodp factorを使わない。Z3については上の五minor circuitがあるが、Z_any65側が不足する。
2. そのcoverのone-parameter有限射影を証明するmonic equations/GBまたは有限basis。generic leading chartだけでなくleading coefficient零の場合も別に保持する。適切なparameterへの有限射影が全curve componentを覆えば、genericparameter上の余核消滅から残りは有限点となる。
3. Z_(p)(parameter)上でも同じfinite modelを作り、数値leading coefficientsがp-unitであることと全localizerのunit性をENTIREspecialmodel上で確認する。inverse変数を入れただけでfiniteと仮定せず、nonproper escapeを閉じる。
4. 同じsourceの実quartic mapをそのfinite algebraで検証する。例えばC122/consistent-frame quotientを使う場合はDeltaH/DeltaQもmodel上のunitとして扱うか、その零caseを別に覆う。GF rankだけをQQへ昇格しない。

現在、この四段を満たす小さい全badcurve coverは旧報告に存在しない。ここで新しいGB、determinant expansion、field、点計算は実行していない。inventory.jsonは五minorの具体indices、factorbag統計、source/archive hashと未完了gateを保存する。全三次・全四次・ambientの未解決性を維持する。

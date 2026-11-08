# rank二conductor plane・全位数の点別被覆：独立合成査読

2026-11-08、別作者 `/root/plane_collision_review`。**純証明の合成査読 PASS**。
対象は rank 四 Poisson 構造の actual 被約既約 quintic D と、一般 restriction rank 二の conductor plane S である。
root の `PLANE_RANKTWO_M2_WHOLE_POINTWISE_COVERAGE.md` と
`PLANE_RANKTWO_CONDUCTOR_WHOLE_POINTWISE_COVERAGE.md` の全分割・仮定・帰結の接続を独立に確認した。
既受理の全源・全逆写像・完全核証明を入力として再使用し、旧行列を新たに再計算したとは主張しない。
新 CAS、operator、kernel、minor、point、profile、GB、jet の実行は零である。

## 1. 最初と最後の補集合

S 上の restriction は非零 section of Λ²T_P²=O(3) なので、零 cubic R は非零。
C 上の全 cubic は integral、line+irreducible conic、三つの相異直線、2L+M、3L に尽くされる。
line-conic は横断／接触、三相異線は triangle／concurrent を区別する。
追加の smoothness、ordinary normal crossing、vertical divisor 不存在は仮定しない。
Integral R は `PLANE_IRREDUCIBLE_CUBIC_RANK_TWO_REVIEW.md` により全 normal order で空。

S は Sing(D) の plane 支持なので generic normal order m は二以上。
全 degree 五式を base/normal に分け、m を最初の非零 normal degree とする。
ある leading amplitude が零となることを元 raw source 全体の空性とせず、次の m へ戻す。
m>=5 なら D は二つの normal variables だけの homogeneous binary quintic であり、C 上で分解する。
したがって対象の既約 D には m=2,3,4 だけが残る。D 恒等零は rank 四でない。

## 2. 位数二：全根型と全 R

G2 は非零 binary quadratic。generic squarefree と generic double が全補集合である。
SF では UFD から G2=h ell1 ell2、h detM=f、primitive row degrees は
00/01/02/03/11/12 の六つ。nonprimitive な極限は共通因子を h に戻した他 flag で扱う。
Double では G2=h ell²、primitive degree は零／一。root basepoints、unused R components、
非被約 R の高位極、非log exact connection を含む元の全入口を保つ。

| 全 m二 source | 全点の十分な帰結 | 使用する強い証明 |
|---|---|---|
| SF、reduced R、全六 flags | MPS 閉和又は空 | `PLANE_RANK_TWO_M2_SQUAREFREE_REDUCED_CUBIC_COVERAGE.md` |
| SF、R=2L+M | 全空 | `PLANE_DOUBLE_LINE_M2_SQUAREFREE_SEPARATE_REVIEW.md` の全 I/II・ratio・Gauss 証明 |
| SF、R=3L、03/11/12 | 全空 | `PLANE_NONREDUCED_M2_SPLIT_ROOT_ENTRANCE.md` |
| 同 02、01、00 | 空、空、essential三八成分又は空 | 02 の curl、01 の全五 source、00 の両 bad-vector と translation の各 receipt |
| Double、moving degree一、nonred R | 全 P0 空 | `PLANE_NONREDUCED_M2_MOVING_DOUBLE_SEPARATE_REVIEW.md` |
| Double、constant、R=3L | essential三八成分又は空 | `PLANE_TRIPLE_LINE_M2_CONSTANT_DOUBLE_ORE_SEPARATE_REVIEW.md` |
| Double、constant、R=2L+M | T36／T38／essential三八成分又は空 | 全34源／51壁の receipt と全八 HF の後続 receipt |
| Double、reduced transverse line-conic | 全空 | 全 connection／bad-vector の double-root receipt |
| Double、reduced tangent line-conic | T37／T40／essential三八成分又は空 | `PLANE_M2_TANGENT_POINTWISE_SEPARATE_REVIEW.md` |
| Double、triangle | 全空 | `PLANE_TRIANGLE_M2_ROOT_PURE_REVIEW.md` と元全 source 逆 |
| Double、concurrent three-lines | 全 T37 又は空 | 全17源／84壁の入口と `PLANE_M2_CONCURRENT_T37_SEPARATE_REVIEW.md` |

SF/reduced の文章中の MPS 所属は SNC 点への actual 全曲線と「SNC と交わる全 component」の theorem を使う。
SNC endpoint の所属だけを元全 family へ外挿していない。

R2L+M constant-double の旧 receipt は八 HF を像23以下としか述べなかった。
後続の八 HF receipt が六全源の proper factor／二つの全 A10 の T36/T38 閉包所属を与えるため、
ここで初めて全点被覆となる。唯一 HA3 nil whole p は、全 p の rank十一 bundle、FULL9 P2、
T零の essential三、非critical の u因子、critical p=-2/7 の元 row34=-7m²/9 を戻す。
p=4 の leading零は本 m二 domain 外へ、chosen-basis p=5 壁等は全 kernel bundle 内へ保持されている。

Reduced double の旧 G 上界も単独では使用しない。接触の QQ(t)／全七壁は有限 shear の全逆と
FULL9 affine P2 を戻し、t=-9/2 の全二係数を T37、t=-3 の全二係数を新しい全 X122/T40 theorem へ送る。
三角形の全十壁では元 row279 の m² を使い、concurrent の h21 は元 m²、二 A8 は実 normalizer と
全 polynomial curve／completed T37 の実恒等式を使う。旧 raw Puv だけの rank 判定は不使用。

## 3. 位数三：旧 G のみの source を残さない

`PLANE_RANKTWO_M3_WHOLE_POINTWISE_COVERAGE.md` を、各後続の強い E 結論と照合した。

| 全 m三 source | 全点の帰結 |
|---|---|
| reduced nonconcurrent の全 root covers／全九 source | 全空。最後の二全 fiber にも六 shear の実逆と FULL9 P2 の因子を使う |
| reduced concurrent の全 SF／repeated／Ore／nonOre | nonOre 空、K0 essential三、K22 nil22-long23、K4a T37、他 nonzero K 空 |
| nonred SF | 既 essential三又は空。最後の c-dependent g2 も全 normal-row 因子で空 |
| nonred generic21、constant00 | 全十四 source・二十二 module・C-exact 全三と全八壁から essential三又は空 |
| 同 moving-double10 | 全十一 source・全壁、critical -8/3 の元二乗係数と full normal-row から空 |
| 同 moving-simple01/02 | 全37 source／22壁／全 actual inverses を戻し、HB は Ga9-distinct23、HA は moving-simple2/15 成分、他は空 |
| nonred generictriple | 正しい minus curvature の全200 routing、全 special/literal 逆から essential三八成分／新二23成分又は空 |

Generictriple の旧 plus-curvature 90分類は撤回されたままで、本査読の入力ではない。
正しい表の二 weak walls (4/21,1/12) は空でなく、新二23成分の全 W六へ戻す。
後続の二 Y walls は元二乗係数と全 actual factor で空となったため、dimension-only source は残らない。
Concurrent では未実行330 root候補の有限列挙を完全性の根拠にせず、受理済み全 strict/Ore theorem を使う。
全 central/cohom/磁気と全 unit 零境界は元全源のままである。

## 4. 位数四：最後の六源以外も明示的に戻す

G4∈O(1)⊗Sym⁴C² の coefficient image rank一では D=L g4+g5 が essential高々三。
既約性から exactly 三となり、全 essential-three 八成分 theorem を適用できる。
Generic partitions 22/4 は固定 multiple-root divisor の Bertini 議論で image rank一へ戻る。
高 rank の残 partition は SF／211／31 の三つである。

* SF の非特別 invariant 枝は image rank一へ戻る。special I0/J0 の全 Witt/Hankel pencils は保持される。
  reduced concurrent の高 rank は全共通外積因子で空、nonred は全六 actualD source の全 Pfaffian 零で空。
  reduced nonconcurrent に残る三 special source は最後の六 shear theorem の前三行で全空。
* 211 は全 fixed-double-root／symmetric-pencil cover を戻す。reduced concurrent は共通外積因子で空。
  reduced nonconcurrent の三源は最後の六 theorem の後三行で空。
  nonred の両 ordering・全 Jordan parameter は `PLANE_M4_TWO_LINE_CASIMIR_SEPARATE_REVIEW.md` により空。
* 31 は R の種類を問わず F0/F1 の二 actualD 正常形へ戻る。F0 は gradient distribution の偶数 rank から空。
  **F1 に旧「像<=21」だけを使ってはならない。** 後続の
  `PLANE_TRIPLE_ROOT_F1_SEPARATE_REVIEW.md` は全必要四 kernel が Zi∧∂w と明示し、
  全 Pfaffian が零となることを独立確認している。これが全31を点別に空へ送る必要な橋である。

最後の六源は全六 shear の定数可逆写像を持ち、任意 P1 と raw P2九を同時に戻した
r_new=r+q(alpha) の全域 triangular inverse を使う。P1 除去後に初めて forcing が零となり、
残全九 P2 の完全核から actualD の proper normal factor 又は零を得る。
homogeneous kernel零を、非斉次 P2 そのもの零と誤って置く議論ではない。
これで m四全 R／全 partition は essential三八成分又は空となる。

## 5. 独立読解の判定と正確な範囲

真の未被覆 source flag は検出しなかった。短い producer manifest の十五入力だけでは
31/F1、R3L SF01/02/00 等への直接の追跡が省略されるため、本 review manifest には対応する
後続の強い receipt を追加固定した。これは旧証明の欠陥を補う新計算ではなく、必要な既証明の接続である。
旧 G-only reports、旧 failure／false／撤回 records は全て不変であり、membership の根拠として代用しない。

各 actual 点は、全 normalizer の逆と全零補集合を保持した有限分割により、現129既知閉成分の有限和又は矛盾へ送られる。
したがって、対象の全 rank二 conductor-plane stratum の pointwise coverage を受理する。新成分増分は零。
これは一般点の非優越だけの命題ではなく、全点の包含である。
rank零 plane、他 factor、全 E、全非被約 scheme 又は総成分数の確定を本 receipt 単独からは主張しない。

有限 pure-proof job の gate と全 input hash は PASS。gate／hash確認は数学の証明ではない。
本査読は合成の独立読解であり、全依存の数千行列を今回再構成したという報告ではない。
Producer 二 notes と旧証明書、global files、主稿、Git を変更していない。

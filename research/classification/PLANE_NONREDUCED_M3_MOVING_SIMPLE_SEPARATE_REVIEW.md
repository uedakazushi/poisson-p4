# 非被約 R・m三 moving-simple 全 degrees(0,1)/(0,2) の独立検証

2026-11-08。producer と別作者が、全正常形・全逆写像・全 coefficient source を独立に監査した。
**本 stratum 全体の G 像次元は二十四以下**を受理する。
弱上界しか得ていない源の E 帰属、generic triple、m二、全分類完了は主張しない。

独立算術は **23.307768 秒、RSS 99,880 KiB、全 PASS**。
`classification_plane_quintic_casimir_review.py --moving-simple-m3` は既存十四 FunctionDef の
AST を全て保存した isolated append で、元 producer／current fullrank runner／CAS を
import・実行していない。Fraction、QQ(τ)、有限 Laurent 恒等式だけを使った。
実行契約、全 input SHA、checker の全 bytes、全 structured 結果は
`CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_SEPARATE_REVIEW_20261108.json`
（SHA256 `c5481bc3fc4ca2697be6438fb75a5eeaa13485ba37536c173c2eca429c9a2e81`）に固定した。

## 1. 全 primitive flag と first connection

actual D は被約既約 quintic、rank二 plane に沿う order は三、zero cubic の支持は
2L+M 又は 3L、normal leading は double+simple とする。
constant double row を z、primitive moving simple row を Az+Cw とし、
`G3=h z²(Az+Cw)`、d=degA=degC は一又は二、deg h=2−d である。
gcd(A,C)=1、h と C の支持は R に含まれる。
全下三角 connection の parallel／flat 条件は

`p=dlog(f/(hC)), r=−d(A/C)+(A/C)(q−p), dq=0`。

q を logarithmic と仮定しない。

f=a²b では全 closed q は
`dlog(f)/3+αdlog(a/b)+γd(b/a)+δd(c/a)`。
d=1,C=a なら r の最高 a 極で γ=δ=0、h二選択×A=b,c で四 whole α line。
C=b なら最高 b 極で h=a の α=4/3、h=b の α=7/3 が強制され、
全 γ,δ を保持する。moving-double first-connection §5 と同じ actual root-row stabilizer は
`δ'=tδ/s, γ'=t(γ−δu/s)/r` を実現するので、δ非零／δ零γ非零／両零を尽くす。
f=a³,d=1 は h=C=a,A=b、global 性により
`q=dloga+γd(b/a)+δd(c/a)`、同じ三場合を持つ。
全逆操作は shear の符号反転と非零 diagonal の逆であり、parameter を chosen slice へ捨てていない。

R2L,C=b の四 EXACT_C whole stratum は、δ非零のまま constant-root 必要式で空となる。
`q(Cbad)=0` の全係数を掛け戻すと

`Ca(ρa ab−γb²−δbc)+Cb(ρb a²+γab)+δCc ab=0`。

ρb は −1 又は −2。bc、a²、ab の順で Ca=Cb=Cc=0、
しかし既約 D は Cbad非零を必要とするので矛盾する。
R の被約性や接続の単極性を流用せず、非log 二方向も保持した直接比較である。
この四源を除く degree(0,1) は四 whole α line＋十一 fixed となる。

## 2. primitive quadratic row の全正常形と零補集合

d=2,f=a²b では C=a²,ab,b²。
元 pole 係数の global 性はそれぞれ次を強制する。

- C=a²：q=−2dloga+3dlogb、A は a 非可除の任意 quadratic。
- C=ab：q=2dloga−dlogb、A は a,b のいずれでも可除でない任意 quadratic。
- C=b²：q=4dloga−3dlogb+γd(b/a)、
  A=k(a²+γab)+m b²、k非零。実 relative row scale と normal Borel で
  k=1,m=0、残 actual a/b ratio で γ=0,1。

最後の全必要式は b=0 上で
`a∂c A=δA, a∂a A=2A, ∂b A=kγa` である。
非零有限次数 polynomial の ∂c 固有値は零なので δ=0、A|b0=ka²、
残 db 係数で上の全 A となる。δ≠0 又は k=0 を黙って単位化していない。

A は modulo C で考えられ、c² 係数、次いで c-linear 各係数の零／非零で場合分けする。
C=a² では、c² 非零の complete square 後
`c²+εab+ηb²` の四型、c²零では
`bc,c(a+b),ac+b²,b²,ab+b²` の五型、計九。
C=ab は complete square 後 `c²+εa²+ηb²` の四型、
c²零で `c(a+b),ac+b²,bc+a²,a²+b²` の四型、計八。
ε,η は零又は一。各非零係数の戻しは非零 complex diagonal による。
例えば c-linear が a+b の時、c→c+ra+sb は modulo C の残 quadratic 二方向を全て除ける。
c-linear が a のみなら、primitive 性が残 b² の非零を強制する。
primitive 性を失う補集合は元の UFD factorization で小さい degree flag へ戻る。
これらに C=b² の二型を加えて R2L は十九 fixed。

f=a³,d=2 は C=a²,p=q=dloga,r=−d(A/a²)。
全 q の高位 exact 項を戻して global 性を比較すると全て零となる。
A/a² は affine quadratic で、primitive 性により quadratic part は非零。
affine GL2 と translations、定数差を用い、rank二、rank一で kernel-linear 零／非零の
全三場合 `A=bc,b²,b²+ac` へ戻る。
全逆 map は上の actual base affine change と normal Borel、非零 scale で存在する。
従って d=2 は全二十二 fixed、前節と合わせ **四 line＋三十三 fixed＝三十七**を尽くす。

## 3. 全 global lift、全 original kernels、全 α 壁

三つの producer report FIRST／ORE／MODULES の埋込 bytes を hash 確認だけし、
全入力式を独立 sparse arithmetic へ読んだ。
各 ω=3fβ−df の全三成分、Euler horizontality、curl/3 trace-free inverse、
全三未知 K の D2 三十係数、P0 actualD と全350 projective Jacobi を再構成した。
二 fixed は global K の augmented rank 増加で空、残りは一意 global K。

全30 P1／全9 P2 の original350 operators、全24 actualD3 map、
leading line の全二十三 annihilator 行、actual six shear を再構成した。
全 constrained kernel は元全行積零、free-coordinate identity、正確な非零 rank minor により完全。
constant-root character `3f q(Cbad)` の全六 base-quadratic coefficient 行も元 ω22+df から作り、
全 kernel を殺すことを再確認した。
四 line の全 generic minor／kernel 分母の根を exact linear-factor 除算で取り尽くし、
合計二十二の元 α 壁を復活した。全三十七と壁を合わせ五十九 source の独立 receipt である。
未知の複素根を有理標本で置き換えていない。

必要 source の全像上界は B+k1+k2+5、固定 α 壁では B を一減する。
B は全 marked root／closed-form parameter を戻した安全な上界で、
一般 orbit の次元を特殊点へ外挿しない。
各 P1 に対する P2 fiber は、元二次 compatibility を全て課す前でも空又は kerJ2 の affine 部分集合である。
三十 whole label はこの段で空又は二十四以下、残七 label を次へ送る。

## 4. 旧七 vertex 診断と全三 bad-column 排除

全七 label の weight零 Euler gauge と fullP1 kernel＋P2全九の positive-weight gauge を再構成した。
六 label と二 line の空でない旧壁では vertex 非依存、
全60 Q-coordinate map／全 W十八 map と Q-space complete rank certificate が一致する。
DEG2_CA2_0 は本当に b 依存が残り、旧 source_covered=false を変更しない。
この全七診断は、空でない旧壁六を合わせ十三 source で確認した。
chosen independent Q/W が元 source に等しいと仮定するものではない。

より弱い条件で、HA_CA_OTHER(t=−8/3)、HB_CA_OTHER(t=4/3)、CA2_0 の全 P1 kernel 末三行が零。
元 P0 の全 z-row は z 倍、P2九は base-base のみなので、全 actual tensor の z-row は z 倍となる。
Euler/Pfaffian の各項から D は z で割れる。
従って三 source は既約 D に矛盾し、CA2 の vertex 障害を解消する必要はない。
全 higher source が rank二だと呼んだものでもない。

## 5. 残四 whole source の全逆橋

残二 line HA/HB と fixed TRIPLE/CB2 の全 P1十三＋P2九に二段 Euler gauge を行う。
全 source が `Q(a,b,z,w)+∂c∧W`、Q四係数＋W十八係数へ入ることを literal に確認した。
四 joint22 coordinate matrix の determinant は全て **−1**。
旧 generic kernel の全分母と source-Q の分母は復活済み壁に支持され、
空でない旧壁は full kernel dimension 十三、同じ full cohom＋common-source が全列を張る。
したがって全元 source の双方向被覆であり、Q の chosen slice からの推論ではない。

全 four Q の strict Jacobi、元 leading target への正確な amplitude、
actual central shear の full polynomial identity を独立に確認した。
HA は k=−(3t+5)/(3t+8)、HB は R=(3t+2)/(3t−4),S=−(6t+1)/(3t−4)。
central change は

- HA：`u[(kδ−α)z²−2βzw]`。
- HB：`u[−γz²−(S+R)βw²]`。
- TRIPLE：`u[−(δ+α+β)z²−2βzw]`。
- CB2：`u[(γ+β)z²+2δzw]`。

全源で exact D3 非零が u非零を強制する。
HA は t=−5/3,−8/3、HB は t=4/3,−2/3 が元 fullleading 像零で空。
残 source の u は actual diagonal/tensor scale で一へ戻せる。
sa=s,sb=sℓ/k,sc=r,sz=k,sw=ℓ,σ=s/r とすれば
HA の u'=s²ℓu/(rk²)、HB は u'=s²u/(rk)。
TRIPLE/CB2 は sb=sa,sw=sz で u'=s²u/(rk)。
この全 Q と W0 の pushforward 恒等式も独立検証した。

HA/TRIPLE/CB2 の central quotient は w² 一方向、HB(t≠1/3) は zw 一方向である。
非零 residual を一へ戻す自由 diagonal character が残る。
HB(t=1/3) では image が z² のみに落ち、F=εzw+ηw² の全四場合を保持する。
u=1 を保つ r=s²/k の下で zw と w² はそれぞれ s/k と s/ℓ の character を持ち、
二つを独立に零／非零正常化出来る。
従って η零を whole line に含め、η一の ε零／一の二 fixed を追加する。
全 inverse は minus shear と非零 diagonal の逆であり、全 W十八へ regauge する。
これで四 whole line＋六 fixed の **全十 module**へ入る。

## 6. 全十 module の元 affine 条件と二つの非空弱上界

全350 affine W十八 matrix／constant、augmented matrix を定義から再構成した。
full kernel／particular の元全行積と rank minors、全 actualD affine 係数を検証した。
全 generic minors、augmented minors、kernel／particular／factor 分母の全根を元 matrices へ戻す。
proper factor は particular だけでなく全 affine kernel 方向を正確に割る。
全十と空でない壁七、計十七 module receipt が PASS。
leading零で既に排除された壁も明示的に保持する。

producer NONDOMINANCE の「全 residual一は不整合、全壁でも同じ」は、
全 report に記録された次の二例外を省略した文章である。
**全結果を空性へ強めてはならない。**

| whole source | full W 核 | proper factor／矛盾／最弱上界 |
|---|---:|---|
| HA_0 generic | 7 | z² |
| HB_0 generic | 5 | z² |
| HB_EXTRA_0 | 5 | z² |
| TRIPLE_0、CB2_0 | 各7 | z² |
| HA_1 generic、HB_1 generic、HB_EXTRA_1、TRIPLE_1、CB2_1 | ― | augmented rank 増加 |
| **HB_0,t=0** | **7** | factor 無し、全像 **23以下** |
| **HA_1,t=2/15** | **6** | consistent、factor 無し、全像 **23以下** |
| HB_0,t=1/3 | 5 | z² |
| HB_0,t=−1/6 | 5 | z |

HB_0,t=0 の上界は B8＋Q-source3＋kernel7＋5=23。
HA_1,t=2/15 は B8＋Q-source4＋kernel6＋5=23。
Q-source は、元 actual normalizer による exact central image と u／残 residual を全て戻した
有限 parameter source の上界である。
他の各 label でも report の image_offset がこの B+Q-source+5 と一致することを確認した。
一般 W 核五を壁の核七へ外挿したり、正規化後一点の kernel だけで whole bound を宣言していない。
二例外の pointwise E 帰属は未解決のままである。

以上で degrees(0,1)/(0,2) の全 G 像は二十四以下となる。
既に別作者検証済みの constant(0,0) と moving-double(1,0) の whole G と合わせれば、
`2 d_double+d_simple≤2` が許す四 flag を尽くし、**非被約 R・m三 generic21 の全 G が非優越**となる。
この aggregate の constant(0,0) 入力は `NONREDUCED_M3_CONSTANT_SEPARATE_REVIEW.md`、
moving-double 入力は `PLANE_NONREDUCED_M3_MOVING_DOUBLE_SEPARATE_REVIEW.md` に明示的に依存する。
generic triple、他 normal order、reduced R、未確認 E 帰属を消していない。

# tangent m二・t=−9/2の全源は既知T37へ：第二のplaneを使う

2026-11-08、`/root/quartic_review`。新しいlocal model／normal jet／profileは不要。
whole二parameter源の実逆と、同じcanonical点が持つ**別のPoisson plane**を使う。
結論は全源の既知T37への正の包含であり、新成分は加えない。
本新しい幾何的適用の別作者読解は作成時点では未了。

## 1. whole源と一度だけのambient

`PLANE_TANGENT_M2_MINUS_NINE_HALVES_WHOLE_ORBIT.md`のP(p,q)について、
q!=0は全てcanonical P*=P(0,1)とactual PGL同値で、q0は同じ軌道閉包に入る。
全inverse、finite coverη²=q、全zero補集合及び元全9P2はその純証明のまま保持する。
P*のDは非零でC絶対既約、generic rank四。

許可された一canonical fullambient jobはgate PASS後、
元全150variation／350Jacobi、25Euler、25GL5＋scalarから

    rank D=101, rank Euler=24, rank M=48

を得た。完全kernelは49／1／3、全maxminorと元全rowでのkernel零を保存した。
射影actual orbit23、接空間24、normal quotient一。
**このrankだけでは成分性又は新成分を宣言しない**。
reportは `CLASSIFICATION_PLANE_TANGENT_M2_MINUS_NINE_HALVES_AMBIENT_20261108.json`。
整数scale2、変数名は(a,b,c,z,w)でz=u,w=v。
2.507617秒、RSS18,472KiB、全150／350／Pf／D*Mを検証した。
ここで追加normal計算を止め、次の全体的なplane幾何へ移る。

## 2. 第二のPoisson planeはconcurrent零cubic・normal位数三

元のtangent planeはS=(u=v=0)だが、P*は

    S'=(c=u=0)

もPoisson planeとして持つ。
実際、Pcj及びPujの全entryは(c,u)で割れるidealに入り、
P*|S'を(a,b,v)で書くと

    Pab=−3b²−9v², Pav=−6bv, Pbv=0.

従ってrestrictionの実zero cubicは

    a Pbv−b Pav+v Pab=3v(b²−3v²).

C上で三つの相異なる直線で、共通点は[1:0:0]。
restriction rankは一般点で二であり、三線の非被約性を仮定する場面ではない。

actual Dの−3/2倍を除いたFは

    −a u⁴+9c(bu+cv)²+8bu³v+9cu²v².

I(S')=(c,u)のorderは正確に三で、leading cubicは

    G3=9v²c³+18bv c²u+9(b²+v²)cu²+8bv u³.

そのbinary discriminantは

    −2916 v⁴(b²−3v²)²

で非零。従ってgeneric normal rootsはsquarefree。
この点で確認したのはrestriction／normal coneの全多項式であり、
元P*を変えてchosen pointのoperatorを新計算したものではない。

## 3. 既受理のwhole concurrent m三定理を適用する

既受理claim `E_PLANE_CONCURRENT_M3_K4A_T37_AND_WHOLE_KNOWN_COVERAGE` は、
rank四・actual D被約既約quintic・rank二Poisson plane・reduced concurrent zero cubic・
normal order三の**全点**を、T37／既nil22-long23／essential-three八成分へ被覆する。
normal root type、全connection、全higher／magnetic complementを含む命題である。
新P*は§2により全仮定を満たすので、そのまま適用出来る。
whole normalizer theoremは `PLANE_CONCURRENT_M3_WHOLE_NONDOMINANCE.md`、
K22の正の認識は `K22_KNOWN_QUINTIC_ORBIT_RECOGNITION.md`、
K4aの既知T37正の認識は `K4A_T37_ACTUAL_PGL_IDENTIFICATION.md` とその別作者full150 reviewを用いる。

この既知閉和への包含だけでpointwise Eは閉じる。
さらにexact T37を選ぶには、残二候補を次で排除出来る。

* Fはessential5。五偏微分の関係にau³,c³v,b²u²,bcu²,u⁴を順に比較すれば独立。
  essential≤3は閉determinantal条件なので八成分に入らない。
* Fの被約特異支持は正確に二plane、{u=v=0}と{u=c=0}。
  全K22はknown nil22-long dense23の一actual軌道そのものへ正規化され、
  そのDの被約特異支持は三plane。したがってP*は全K22 actual軌道でない。

残るのは全K4a actual軌道であり、既T37の閉包に入る。
特異plane数が一般の退化で保存されるという仮定は使わない。
ここで比べたK22は**全caseが一つのactual orbit**という既whole定理である。

## 4. 全原sourceへ戻す

全q!=0はP*の同じactual軌道、T37は閉PGL安定成分なので全unit枝はT37内。
全(p,q)はaffine polynomial Poisson族で、q!=0は稠密。
よって全q0もT37内である。
元の六shear／全source正常化のactual逆を戻してもT37内への包含は保たれる。

従ってt=−9/2の**全moving-double source**は既知T37で被覆される。
同じambientを使った新成分count、normal obstruction/model又は追加profileは必要ない。
t=−3の残二parameter源は別の未解決義務として保持する。

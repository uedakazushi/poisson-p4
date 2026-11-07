# H+既約四次曲面：非零半単純cofactorの全源支持

2026-11-07。課題 `E_reduced_nonnormal`。
本稿の型番号は `CLASSIFICATION_HYPERPLANE_QUARTIC_SEMISIMPLE_SOURCES_20261107.json`
の24行に対応する。これは全支持を覆う表であり、24型が相異なるとは主張しない。
型18と19は全scalar倍・座標交換で同値である。

入口の有限化は `HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE.md`、その独立検証は
`CLASSIFICATION_HYPERPLANE_QUARTIC_ENTRANCE_REVIEW_20261107.json` を参照する。
全90係数の各固定cofactorに対する `div P=Z` と `[Z,P]=0` の全線形源を構成し、
全四Jacobi、多項式Pfaffian、最高次Pfaffianを保持した。
別のFraction実装による全24源の独立再構成は
`CLASSIFICATION_HYPERPLANE_QUARTIC_SEMISIMPLE_SOURCES_REVIEW_20261107.json` にある。
以下の積恒等式・重みは reusable runner の `--support` mode、
`CLASSIFICATION_HYPERPLANE_QUARTIC_SEMISIMPLE_SUPPORT_V2_20261107.json` に保存した。
計画checkerのPASSは数学の証明ではない。

## 1. 結論と射程

型1,2,3,5,6,8,9,10,11,13,14,15,16,18,19は、全線形源の実Pfaffianが
固定座標因子を持つ。従ってJacobiを課す前から既約quarticを持たない。
型4,17,21,22,24も、全Jacobiの被約支持を以下の分岐で調べると可約又は
最高四次部分零となる。
型7,12の全点は、既約quarticの条件の下でMPS40閉包の和へ入る。
型20は全定数crossing階数分岐・七つの正重みSNC曲線でMPSへ、型23は全支持可約性で排除した。
その全証明は `HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md` を参照する。
H+既約quartic全体の被覆や全分類の完了は主張しない。

全変数は複素数体上で扱う。Pfのaffine表示を次数4に斉次化する。
従ってaffine Pfの次数が4未満なら斉次化にh因子が生じ、単純Hと異なる
既約quarticという仮定に反する。

## 2. 型4と17

型4の全Jacobiは、非零定数倍の重複を除き

    n0(n4−1)=0,  n3(n2−n5)=0,
    n0 n3−n1(n2−1)=0

となる。n0=0ではPfはz因子、n3=0ではy因子を持つ。
n0 n3≠0ならn4=1,n5=n2,n0=n1(n2−1)/n3であり、実Pfは

    (n1 y²−2n3 xz) (z−(n2−1)wy/n3)

と等しい。全開分岐と両閉補集合を保持した。

型17の全Jacobiには `(2n5−1)²` と `n3(3n5−1)` が含まれる。
被約支持上でn5=1/2,n3=0となり、実Pfはwyで割れる。
この議論は根基上の結論であり、元の非被約理想の等式ではない。

## 3. 型21、22、24

以下v=wx、f=n0+n1z+n2z²、g=n3+n4z+n6z²と置く。
全分岐の正確な置換と二因子は証明書の `factorization_certificates` に記録した。

型21ではn5=0,n9=0ならPf=(f+v)g。
n5=0,n9≠0なら全Jacobiからn3=n4=0となり
Pf=z²[n6 f+(n6+n9)v]。
n5≠0ならn9=1/4,n2=0,n1=n4/(4n5),n0=n3/(2n5)であり、残式は

    n4 n3=0,  n4(3n6+1)=0,
    n3(2n6+1)+n4²=0

を含む。最終式にn4を掛けるとn4³=0なので被約支持上n4=0。
n3=0ならv因子、n3≠0ならn6=−1/2で

    Pf=(n3/(2n5)+v/2)(n3−z²/2+n5v)

となる。n4³の理想所属も小さい多項式恒等式で独立に検査できる。

型22では全Jacobiから

    n3=0, n4 n5=0, n1(4n4+1)=0,
    (2n4+1)n0+2n1n2=0

を得る。t=z²−2yと置くと

    Pf=(n0+n1t+v)(n2+n4t)−n5(2n4+1)v.

n4≠0ならn5=0で積。n4=0ならn1=n0=0で(n2−n5)vとなる。

型24のn5≠0部分ではn10=1/2,n6=−2n2n5であり、n2n6=0からn2=0、
次いでn6=n1=0となる。さらにn0(n4+1)=0。
n0=0ならv因子、n0≠0ならn4=−1で

    Pf=(v+2n0)(n5v+n3−z)/2.

n5=0,n6≠0部分ではn10=n9=n2=n1=n0=0となりPf=vg。
n5=n6=0では、全Jacobiからn10n3−n4n9=0であり

    Pf=f(n3+n4z)+v[n3+n9+(n4+n10)z].

二つのz一次式は比例する。n4≠0なら

    Pf=(n3+n4z)[f+(n4+n10)v/n4].

n4=0,n3≠0ならn10=0、Pfの次数は2以下。
n4=n3=0ならPf=v(n9+n10z)。これで全閉補集合を尽くした。

## 4. 正則paired型7の全点をMPSへ接続する

型7のcofactor固有値は(−3,−1,1,3)。全Jacobiは次の七式と同値である。

    n4(2n6+2n8−3)=0,  n2n4=0,
    n1n4−n3(n6+n8)=0,  n0(n6−1)=0,
    n7(n5−n8)=0,
    2n0n7−2n2n5+n2=0,
    n1(1−n5)−n2n3=0.                         (7)

κ=3n5−n6−n8と置く。n2n4=0により二源に分け、表の各係数を
n_j→t^(weight_j)n_jとする。他の対角係数n5,n6,n8は固定する。

| 部分 | 正の重み |
|---|---|
| n4=0 | (n0,n7,n2,n3,n1)=(1,1,2,1,3) |
| n2=0 | (n0,n7,n4,n1,n3)=(1,1,1,1,2) |

(7)の各式はこの重みで斉次なので、全tに対して全Jacobiが成立する。
t=1は元の全源点、t=0はlogcanonical tensorでありPf=κxyzw。
κ≠0なら実射影退化因子はhxyzwで、五つの独立超平面のSNCである。

κ=0の補集合を捨てない。a=n0n7、E=a−n2(n5−1/2)とすると

    a²=aE+(n2/2)[aκ+n7 n0(n6−1)−n0 n7(n5−n8)].       (8)

従って(7),κ=0の被約支持上でn0n7=0。
n3≠0ならn4≠0は不可能である。実際n4≠0ならn6+n8=3/2、κ=0よりn5=1/2、
n2=0、最終式からn1=0となり、第三式がn3=0を強制する。
よってn4=0、続いてn5=n2=n1=0。
この全枝はn0=0ならx因子、n7=0ならw因子を持つ。
n3=0,n4≠0ならn5=1/2,n1=n2=0となり、同じくx又はw因子を持つ。
最後にn3=n4=0ではn0=0ならz因子、n7=0ならy因子を持つ。
これら六分岐の積恒等式を全て記録した。
従って既約quarticを持つ全点ではκ≠0であり、上の真のSNC曲線を持つ。

## 5. 正則paired型12の全点をMPSへ接続する

型12の固有値は(−2,−1,1,2)。全Jacobiは

    n4(n7+n8−1)=0,  n2n4=0,
    n1n4−n3(n7+n8)=0,  n0(n7−1)=0,
    n6(n5−n8)=0,  n2(2n5−1)=0,
    n0n6+n1(1−n5)−n2n3=0.                    (12)

κ=2n5−n7−n8と置く。全点に対し次の二重みが真のPoisson曲線を与える。

| 部分 | 正の重み |
|---|---|
| n4=0 | (n0,n6,n2,n3,n1)=(1,1,1,1,2) |
| n2=0 | (n0,n6,n4,n1,n3)=(1,1,1,2,3) |

各残Jacobiは斉次であり、t=0のPfはκxyzw。
κ≠0は型7と同じSNC接続を与える。

κ=0,n4≠0ならn2=0,n7+n8=1,n5=1/2。
もしn0n6≠0ならn7=1,n8=n5=1/2となり和が1に反する。
従ってn0n6=0、最終式と第三式からn1=n3=0となる。
n0=0ならx因子、n6=0ならw因子を持つ。
κ=0,n4=0,n2≠0でもn5=1/2,n7+n8=1であり、同じ議論から
n0n6=0,n1=n3=0。今度はn0=0ならz因子、n6=0ならy因子を持つ。
最後のn2=n4=0では最高四次Pfはκxyzw=0なので、単純Hの仮定に反する。

## 6. SNCからの帰属に使う定理と実例

上の曲線の像閉包は既約であり、一つの全Poisson成分Cに含まれる。
κ≠0ではCは実際にSNC log-symplectic点を持つ。
ここで用いる外部入力はMPS Theorem1.4の**SNC locusと交わる全成分が既知40個**
という尽くしである。従ってCはMPS40の一つであり、元のt=1点もその閉包へ入る。
これはSNC点が既知の閉和に属するだけという推論ではない。
T01の一意性も要求しない。
既存の同じ適用の詳細は `DOUBLE_QUADRIC_SUPPORT.md` のSNC成分定理適用箇所と
`KNOWN_COMPONENT_AUDIT.md` のMPS一次文献監査を参照する。

型7には実際に既約quarticを持つ点がある。

    n0=1,n2=n7=0,n1=3/2,n3=n4=n5=n6=1,n8=1/2

では全四Jacobiが零で、実Pfは

    [3(wx+1)(wx+yz+1)−2wy³]/2.

これはzの一次式で、係数(3/2)y(wx+1)と
定数項(3/2)(wx+1)²−wy³はC[x,y,w]上で互いに素である。
前者の既約因子y,wx+1のどちらも後者を割らないので、Gauss補題により絶対既約。
斉次化にもh因子はない。
従って「全半単純型が既約quarticを持たない」という強い主張は誤りである。
この点はκ=3/2で上の全点SNC曲線を持ち、MPSに被覆される。
新成分の存在を意味しない。

## 7. 純冪零入口と残る義務

同じ全90係数runnerの `--nil-sources` modeは純冪零10型も構成した。
`CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json` に全源・全Jacobi・
実Pfを保存した。型ごとの全線形源次元は

    (4,0):13, (4,4):10;
    (31,0):19, (31,1):18, (31,3):15;
    (22,0):27, (22,2):20;
    (211,0):36, (211,1):28, (211,2):25.

括弧の第二成分はtranslationのcoker-Nでのtop block長で、0はtranslation零。
最初の試行は未制約Pfの因数分解で180秒上限へ達したため、全支持に不要な
未制約因数分解を省略し、全rawPfを保持する同じrunnerで13.4秒で完了した。
これは純冪零のJacobi支持を解いたことではない。
半単純型20/23の全支持は上記別稿で閉じた。純冪零型の残帰属は別途全支持を扱う。

## 型20/23の独立読解の範囲

別担当の全支持証明を、入口・全90係数源を構成した担当が独立に読解した。
C-rank0/1/2、非零平方rのunipotent正規化と有限矛盾、q=r=0の対称rank1行列と
Cartan射影、r=tzwの九式と全generic/special値、七正重みSNC曲線を照合した。
generic tでE3から直接e(1−h)=0とする記述には、先にE5/E7からbc=0を導く
説明が必要であると指摘し、別稿で修正された。具体的にはbc非零なら
t=1−2e=1−2hなのでe=h=1/2,t=0となりgeneric仮定に反する。
この確認は独立な証明読解であり、第三の算術実装・全計算再実行ではない。
全90係数源の既存独立Fraction再構成とは区別する。読解時点のsnapshot SHA256は以下。

* code/classification_hyperplane_quartic_semisimple_reduction.py: `d8a46fbeeff7a327835674f7f9650b731c7caccc46fc4dacceaa9e829b220719`
* reports/CLASSIFICATION_HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_COMPLETE_20261107.json: `6834df5e8a4274f37b48f3473811270cd13c39dabb7b6063468fbdc35ffd2043`
* research/classification/HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md: `c6ba0d9a0ae1ba30ced7082b104d4692f665316e1283fc561abcd823fb969e15`

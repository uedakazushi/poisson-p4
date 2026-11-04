# actual D の第四の embedded associated closure：B18

2026-11-04。D は全射影被約成分 CA と CB の実 scheme 交叉であり、全 Poisson
Jacobi scheme 自体ではない。独立 source 監査と独立初期関係・深さ監査により
元 P0 での cone depth19、projective depth18が確定した。本稿は、その全 scheme
群積と実軌道を使う global associated-prime corollary である。新しい係数・点・
jet・GBを計算しない。

## 1. P0 の軌道閉包は既知 B18 そのものである

点は既存の

    P0=partial_y wedge D0,
    D0(h)=x^2, D0(x)=0, D0(z)=hx, D0(w)=h^2+z^2。

全25 GL5微分の独立計算により、GL5とtensor scalarの軌道 O0 の cone 次元は
19である。特に E_yy P0=P0 で scalar方向は実軌道に入る。作用群は連結既約
なので O0 は smooth irreducible locally closed orbit、閉包も既約19次元。

既受理の全被約境界 B18=(R intersect N)_red は既約 projective18次元であり、
cone次元19である。その全 chain source は

    D(h)=a*x^2, D(x)=0, D(z)=b*hx+c*x^2, D(w)=G2(h,x,z)。

P0 は literal に a=b=1、c=0、G2=h^2+z^2 の source点である。B18 は作用群で
安定なので closure(O0) subset B18。既約性と同じ次元19により

    closure(O0)=B18。

従って O0 は B18 の稠密開軌道である。B18 の名前を接空間の似た別の formal
層へ付けたのではなく、全 proper boundary theorem の同じ閉像と同定している。

## 2. depth0 横断環の maximal associated prime

受理済みの actual 全群積は、nilを含めて

    completed O_(D,P0)=A_perp[[g0,...,g18]]

を与える。A_perp は9次元、embedding dimension20の complete Noetherian local
ring。19群方向は正則列であり、新しい depth19 theorem により depth(A_perp)=0。
ゆえにその最大 ideal m_perp は Ass(A_perp) に属する。

実際 H0_m(A_perp)!=0 の有限長 module には非零socleがあり、その非零元 z は
Ann(z)=m_perp を満たす。具体的な z の source式や H0の長さは不要であり、
ここでは求めていない。同じ z を群積へ定数係数で入れると

    Ann(z in A_perp[[g19]])=m_perp[[g19]]。

群変数の非零冪級数はこの class を殺さない。したがって completion に
19次元 quotient k[[g19]] を持つ associated prime が得られる。

## 3. これは actual orbit ideal の completion であり、globalへ降りる

上の群積は任意に選んだ抽象 factorization ではない。実 canonical126空間で
19群列の minor -32 と107横断補空間を用いる actual作用写像

    Gamma19 times (P0+H107) -> canonical126

の formal IFT同型である。Gammaの最初の変数は literal tensor scalar、他は
実GL方向。CA、CB、その和ideal Dを保つことは全Artin基底上で証明済みで、
同じ積が N,R にも使える。

横断点を P0 に固定した像はまさに軌道 O0 の形式近傍。従って
m_perp[[g19]] は algebraic orbit ideal I_(B18,P0) の completion である。
この明示的同定により、任意の formal-fiber prime の contraction の次元を
無根拠に同じとする問題を避ける。

faithfully flat completion の associated-prime contraction から

    I_(B18,P0) belongs to Ass O_(D,P0)。

localization の associated-prime 定理により、global cone D の associated
prime も I_(B18) そのものである。B18 は minimal N28 と R20 の両方より
次元が小さく、その proper 境界なので、この associated prime は embedded。

## 4. 第四の相異なる閉包と projective版

受理済み Pa、Qb、CA-critical V の三つの embedded associated closure は全て
cone20/projective19次元であり、互いに異なる。今回の B18 は
cone19/projective18なので、これら三つのいずれとも異なる。

全閉包は scalar-stable cone で、P0では q106=1。非零 homogeneous chart では
cone ring は projective chart と Laurent scalar の積となる。同じ associated
prime が scalar factor を除いて降りるため、projective actual D にも B18 が
embedded associated closure を与える。

従って既知の global embedded associated prime closure は少なくとも四つ：
Pa、Qb、CA-critical V（各 projective19次元）と B18（projective18次元）。
これは全 associated-prime list、primary ideal、各長さ・重複度の尽くしではない。
元 D0 の具体的socle代表、H0長、nil三乗の非零性や全nil指数も決めていない。
特に H0!=0 から Nil^3!=0 を推論しない。

Poisson Jacobi scheme の irreducible component を追加する主張ではなく、
全Poisson成分下界58は変更しない。既存の frozen manuscript・台帳は編集せず、
本稿と同名 report だけに証明とexact依存SHA256を保存する。

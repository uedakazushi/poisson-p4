# tangent m二・t=−9/2：全残二源を一実軌道へ

2026-11-08、`/root/quartic_review`。純構造証明、新しいM1／profile計算無し。
全source逆は `PLANE_M2_TANGENT_MOVING_WHOLE_SOURCE_REDUCTION.md` §4に固定する。
本稿は同節の**全(p,q)**を扱い、q!=0の全像を一actual軌道へ、q=0をその閉包へ送る。
canonical軌道の次元下界、成分性又は既成分名は次の別義務である。

## 1. fullsourceと実Borel

座標(a,b,c,u,v)、全P(p,q)は

    Pab=3ac−3b²−p uv−9qv²,
    Pac=u(pu+8qv), Pau=−3(2bu+cv), Pav=−3(au+2bv),
    Pbc=q u², Pbu=−3cu/2, Pbv=−3(2bu+3cv)/2,
    Pcu=Pcv=0, Puv=−3u²/2.

全(p,q)で元350projective Jacobiは恒等零。
全9P2を戻した後の完全解であり、全affine particularを省略していない。
実座標変更

    a'=a+2λb+λ²c, b'=b+λc, c'=c,
    u'=u, v'=v−λu                                       (1)

はdet1、逆はλ→−λ。
P0(t)は全tでこの(1)にliteralに不変である。
conic ac−b²とroot bu+cvもそのまま保たれる。
全残P2に対しても、十括弧のLeibniz則から厳密に

    gλ_*P(p,q)=P(p+10λq,q).                             (2)

例えばPacへの追加λ項は8λq u²+2λq u²、Pbc=q u²は不変。
Pabへのu²項は−9λ²q+8λ²q+λ²qで相殺され、uv項は−10λq。
他の追加括弧は零又はP0不変なので、(2)はfulltensorの恒等式であってDだけの比較ではない。

## 2. 全unit枝と全inverse

q!=0ではλ=−p/(10q)によりp'=0。
次いでη²=qを選び、actual diagonal

    (a,b,c,u,v)→(a,b,c,ηu,ηv)

を使うとP0は不変、(p,q)→(p/η²,q/η²)である。
したがって全unit枝はP(0,1)へ戻る。
inverseはη^-1、λの符号反転で明示的。平方根はC*の有限surjective double cover上の選択だけで、
新しいsource parameterやdiscarded componentではない。

q=0ではactual Dは

    −(3/2)(bu+cv)[9c(bu+cv)+p u³]

と分解する。既約quinticの条件下では空である。
さらにP(p,q)はaffine A²上の全Poisson polynomial族、q!=0は稠密で、
全unit像は上の一軌道に含まれる。従って全q0もこの**同じ軌道閉包**に含まれる。
canonical軌道に帰属する既知／新成分が分かれば、全sourceへ閉性で戻せる。
boundary所属からgenericを推測する逆向き論法ではない。

## 3. canonical点のactual既約性・rank四

P*=P(0,1)の実Dは−3/2倍を除いて

    F=−a u⁴+9c(bu+cv)²+8bu³v+9cu²v².

aの一次式で先頭係数−u⁴、残項のu=0制限は9c³v²!=0。
したがってC[b,c,u,v]でも両係数は互いに素で、Gaussの補題からFはC上絶対既約。
全Pfは非零なのでprojective一般rank四。

∂aF=−u⁴からSing Fはu=0に入り、そこで∂vF=18c³v。
逆包含も全偏微分から確認出来るので、被約特異支持は二つの平面

    {u=v=0} ∪ {u=c=0}.

{u=v=0}の一般点に沿うnormal quadraticは9c(bu+cv)²で**double root**。
これらの幾何を先に確認してから、一canonical点の全ambientによる成分性／親義務へ進む。
rankや成分数を本純証明だけから宣言しない。

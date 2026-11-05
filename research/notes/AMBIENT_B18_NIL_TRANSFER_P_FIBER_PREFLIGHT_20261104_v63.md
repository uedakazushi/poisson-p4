# transfer W の P 被覆：正しい九次元 honest fiber と不足 source

2026-11-04。元 D0、N、K、P=N/K、D=k[[ell]] は既存 packet と同じである。ここでは W→P の実値を計算せず、利用可能な512 sourceと必要な exact fiberを指定する。

## 1. 被覆した場合の正確な帰結

P=psi(W) なら N=K+W。既知 K³=0 と WN²=0 により全ての三次積が零なので、これは **N³=0** を直接与える十分条件である。gamma=0からの上界四より強いが、被覆自体は未証明である。

P/psi(W) は D-finiteなので、Nakayama による正しい有限 test は

    P/(psi(W)+ellP)=0。

P は D-free7 のため、0→K/ellK→N/ellN→P/ellP→0 は exact。honest Tnil の既受理 sequence と合わせると

    S9 := Tnil / image(K/ellK),
    0→P/ellP → S9 → Cprime[ell]→0,
    length(S9)=9, length(P/ellP)=7, length(Cprime[ell])=2。

従って S9 を512 actual W像の **全 R0/(ell)-submodule** で割った長さが2なら被覆し、nilcubezeroが従う。originの定数 rankだけではこの長さを証明しない。S9 の actual表示・K/ellK注入はまだ構成されていない。

## 2. 正規化 ell-fiber にある不可視の二〜四方向

NE=Nil(E4) は D-free9、Q=NE/P=ker(C→Cprime)。既知 C のD-rank/fiberは4/8、Cprimeは2/4なので、tors_D C の minimal数は4、tors_D Cprimeは2である。

0→Q→C→Cprime→0 から Q[ell]→C[ell]→Cprime[ell] が exact。QのtorsionはCのtorsionの部分なので、そのgenerator数 t≤4。他方 dim C[ell]=4、dim Cprime[ell]=2より t≥2。従って

    2≤t=dim Q[ell]≤4。

0→P→NE→Q→0 の tensor sequence は

    0→Q[ell]→P/ellP→NE/ellNE→Q/ellQ→0。

ゆえに P/ellP→NE/ellNE は必ず非単射で、その kernel は二〜四次元である。可視の像は三〜五次元に過ぎない。NE/ellNEの中で七次元の W像を求める test は不可能であり、七 generator存在・rankだけをもって primitive lattice としない。

また whole old specialized comparison T_K の長さは L、K/ellK は L−4。T_K 全体で Tnil を割ると長さ五となり、必要な四つの C[ell] edge方向を落とす別の quotient である。K/ellK の代わりに用いない。

## 3. 512 actual source の利用可能性

transfer packet は全512の actual44-cycle root、元 C0 free9代表、正規化像 f(n_alpha e_gamma) の C-four-component sourceを保存している。NEは二 nil generator×四 C-basisによりR0-moduleとして張られるので、これら512は WのR0-module像も生成する。

各 actual cycle z を ell=0 に制限し、v=mzへ original canonical three-sector Koszul homotopyを適用して d2 b=v を取れば、honest128代表(z,b)を得る。Phi(z,b)=z。その後、既受理 minimal5×128 target projectionを適用できる。homotopyはS-linearでないため、512個の original formal functionそれぞれへ適用してから線形に延長する。この source recipe は利用可能だが、新しい512 homotopy operator/valueはここで作らない。

不足しているのは、(a)Pのactual七D-basisと512像の正確な regular coordinate/lattice証拠、または(b)S9を作るhonest K/ellK injectionと、512像をそこへ移す完全な比較である。既存のTnil5×849や全849一次値だけでこの注入を補わない。新source/jetを自動で要求せず、all-order edge-map/source gateを先に閉じる。

本 packet は被覆、実 Qの t、beta/gammaの消滅、nil三乗零のいずれも決定しない。

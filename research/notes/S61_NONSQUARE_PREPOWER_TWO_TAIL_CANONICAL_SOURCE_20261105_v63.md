# 二 tail の共通分母・局所 source 表示と最小の既知単元 case

2026-11-05。係数値の新しい計算・factorization・GCD を含まない SOURCE。
元の T0/T1 係数と独立受理済み商を結び、同じ共通分母を使える点を固定する。

\[
T_0=a^4C^4J_0/D,\qquad T_1=a^3C^4J_1/D.
\]

両 D は文字通り同じ197 bit整数。J0は4500項、J1は4811項。
通常次数はそれぞれ49/50、boxは (16,32,18)/(16,33,18)。
\(R=\mathbf Q[a,e,d,1/(adFBCJ)]\) では scheme として

\[
(T_0,T_1)R=(J_0,J_1)R.
\]

従って既知の十分開は \(D(D_{12})\cap(D(J_0)\cup D(J_1))\) と書ける。
残る \(V(D_{12})\cup V(J_0,J_1)\) はこの十分開の補集合であって、
実 bad Fitting scheme と同一と主張しない。

## 既存 full source での二列表示

pair 四行は current15 の [2,4,5,8] で、同一の保存済み κ を使う。
\(T_j=B\kappa M_j\)（j=1,2、T2の名称をT0とした）から

\[
\Phi=\frac{DB}{a^3C^4}\kappa,\quad
v_1=M_1,\quad v_0=a^{-1}M_2,\qquad
\Phi v_1=J_1,\quad\Phi v_0=J_0.
\]

これは既知 FULLU unit による二列の可逆変更だけである。
任意の parameter-only α,β に対し、source 全行で同じ列変更を施せば
\(\Phi(\alpha v_1+\beta v_0)=\alpha J_1+\beta J_0\)。
新しい混合列の値は計算していないし、他の tail をこの二本が生成するとも言わない。

受理済み共通十五角の source/target 順序・倍率を保ち、
\(E_1=D D15_1/(a^3C^4)\)、\(E_0=D D15_0/(a^4C^4)\) とすると
\[
A E_i=-J_iD_{12}
\]
が source domain の恒等式である。A、κ3、旧Ω因子の反転は導入しない。
全330 source roots、H3、他の七 tail と full57 lift を保持する。

## 次に選べる最小の既知単元 case

保存 support だけで両 J の a,d 付値は0。C の付値も両原 tail の正確な付値4から0。
独立 ROOT Horner により \(v_F(J_1)=0\) で、共有する既知 F 因子もない。
未処理の既知 FULLU 因子は B,J。小さな原式から
\(B+10J-Cd=-80\) を全係数で確認するが、これは補集合消滅を意味しない。

J=5ad−6e−4 の d 根は \((6e+4)/(5a)\) で、先頭係数5aは既知 FULLU unit。
従って **一回だけ J1 の J 付値** を調べる proposal は未知の分母 chart を増やさない。
固定 \(D5^{18}a^{18}\) を払い、k=0..18の最初の非零まで、二変数層だけを保持できる。
層支持上限1667、root cache190、whole owned29336（提案34376）、
L1上限372 bit（hard518）、全 product上限120万、RAW8 MB。
これは SOURCE proposal のみであり、J の値・新 worker・商・GCD は未実行。

この既知単元を除く目的は将来の二 tail の補集合計算を小さくすることである。
普通の多項式因子を除外できても局所化後の unit ideal、radical coprimality や全 FULLU57は従わない。

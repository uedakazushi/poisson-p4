# 保存matrix-vectorによるmixed16-minorのg13判定

2026-11-04。一回の保存入力Nm1*Vnormal mod g13と、そのgcd/Bézoutだけを計算する。field逆元、e-Euclid、新しいnorm、点、素数、GBは使わない。

旧V=adj(Nk)e0は受理済みのdetとNV証拠から定まる。全16成分共通のQQ contentとa24だけを除いたVnormalを使い、成分ごとの独立scalingはしない。保存mod g13入力をnumeric共通分母DA,DVで整数化し、全16個のraw整数多項式行列-vector積を取る。

整数g13の最高係数Lに対し、各行をpseudo-divisionしてL^d*raw_i=g13*Q_i+R_iを得る。全商・剰余・指数・積恒等式を保存し、別のPython整数convolutionで再検査する。実際のW_iはR_i/(DA*DV*L^d)で、逆にするのはQQ定数だけ。全剰余次数は13未満。

Nm1とNkは同じmonic h1代数での乗法なので可換し、adjugateとも可換する。従ってunitを戻したW_iはNkの列iをNm1の列0で置き換えたmixed16-minorである。Fitt0(T)はg13と全W_iを含む。gcdがQQ単元ならT=0で旧D(G)badは空。全Wが零やproper gcdなら候補を残し、非空・e有理・length13を自動的に結論しない。

gcdの各変化には両quotientと整数Bézout恒等式を保存する。G0のQ/L側や全Bbadへの組合せは親の別gateであり、このprimaryが全域空性を単独で主張するものではない。60秒/2GB/32MiB/150k項を上限とし、別methodへ自動切替しない。

## 確定結果

承認された一回はnative0.0539秒、全体1.0958秒で完了した。**全16成分のWがg13を法として厳密に零**だった。16個のraw整数積と全pseudo-division恒等式を保存し、別Python整数convolutionで全係数一致を確認した。実際に必要だったLの冪は全16本で0、従って整数のままraw_i=g13*Q_iが成立した。raw積は合計391項、全商は183項である。

従ってgcdはg13のままであり、単位criterionは成立しない。13次candidateを保持し、actual Bbadの実在、空性、e有理性、exact length13は今回結論しない。追加のfield inverse、e-Euclid、別normや方法への切替は行っていない。
